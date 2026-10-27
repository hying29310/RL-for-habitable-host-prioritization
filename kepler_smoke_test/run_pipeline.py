"""One-star engineering smoke test, not trained RL or planet confirmation.
Run: python run_pipeline.py (downloads missing data; cached reruns are offline).
"""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplconfig'))
import json
import time
import hashlib
import platform
from datetime import datetime, timezone
import requests
import numpy as np
import scipy
from scipy.ndimage import median_filter
import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.timeseries import BoxLeastSquares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA = ROOT / 'data'
URLS = {
    'target_pixels.fits.gz': 'https://archive.stsci.edu/missions/kepler/target_pixel_files/0069/006922244/kplr006922244-2009166043257_lpd-targ.fits.gz',
    'pipeline_lightcurve.fits': 'https://archive.stsci.edu/missions/kepler/lightcurves/0069/006922244/kplr006922244-2009166043257_llc.fits',
}
QUERY = "select pl_name,hostname,ra,dec,st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,st_metratio,st_rad,st_mass,pl_orbper,pl_tranmid,pl_rade,pl_orbsmax,st_refname,pl_refname from ps where pl_name='Kepler-8 b' and default_flag=1"

def save_json(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')

def download():
    DATA.mkdir(exist_ok=True)
    for name, url in URLS.items():
        dest = DATA / name
        if not dest.exists():
            r = requests.get(url, timeout=90)
            r.raise_for_status()
            with fits.open(__import__('io').BytesIO(r.content)):
                pass
            dest.write_bytes(r.content)
    if not (DATA / 'catalog.json').exists():
        r = requests.get('https://exoplanetarchive.ipac.caltech.edu/TAP/sync',
                         params={'query': QUERY, 'format': 'json'}, timeout=90)
        r.raise_for_status()
        rows = r.json()
        assert len(rows) == 1
        save_json('data/catalog.json', {'query': QUERY, 'request_url': r.url, 'rows': rows})

def best_result(pg):
    j = int(np.argmax(pg.power))
    return {k: float(pg[k][j]) for k in
            ['period', 'duration', 'transit_time', 'depth', 'depth_err', 'depth_snr', 'power']}

class SearchEnvironment:
    """Minimal reset/step environment with hidden period-and-epoch evaluator.
    Fixed abstract costs test accounting only; wall time is recorded separately.
    No star selection, learning, HZ reward, or full false-positive vetting yet.
    """
    costs = {'coarse': 2, 'refine': 3, 'submit': 1, 'stop': 0}

    def __init__(self, t, y, dy, features, evaluator, budget=6):
        self.model = BoxLeastSquares(t, y, dy)
        self.features = features
        self._evaluator = evaluator
        self.initial_budget = budget
        self.reset()

    def reset(self):
        self.budget = self.initial_budget
        self.candidate = None
        self.history = []
        self.done = False
        self.refined = False
        return self.observe()

    def observe(self):
        actions = ['stop']
        if self.candidate is None:
            actions += ['coarse']
        else:
            if not self.refined:
                actions += ['refine']
            actions += ['submit']
        return {'stellar_features': self.features, 'budget_remaining': self.budget,
                'candidate': self.candidate,
                'valid_actions': [] if self.done else [a for a in actions if self.costs[a] <= self.budget]}

    def step(self, action):
        if action not in self.observe()['valid_actions']:
            raise ValueError('Action unavailable, unaffordable, or episode finished')
        start = time.perf_counter()
        reward = -0.05 * self.costs[action]
        evaluation = None
        if action == 'coarse':
            # Fixed search grid; no known period or epoch enters the search.
            self.periodogram = self.model.autopower(
                [0.08, 0.12, 0.16, 0.20], minimum_period=0.5,
                maximum_period=10.0, frequency_factor=2.0)
            self.candidate = best_result(self.periodogram)
        elif action == 'refine':
            p = self.candidate['period']
            self.candidate = best_result(self.model.power(
                np.linspace(p * 0.99, p * 1.01, 2001),
                np.linspace(0.06, 0.24, 19), oversample=20))
            self.refined = True
        elif action == 'submit':
            evaluation = self._evaluator(self.candidate)
            reward += 1.0 if evaluation['period_and_epoch_match'] else 0.0
            self.done = True
        elif action == 'stop':
            self.done = True
        self.budget -= self.costs[action]
        self.done = self.done or self.budget == 0
        entry = {'action': action, 'abstract_cost': self.costs[action],
                 'wall_seconds': time.perf_counter() - start, 'reward': reward,
                 'budget_remaining': self.budget, 'candidate': self.candidate,
                 'evaluation': evaluation}
        self.history.append(entry)
        return self.observe(), reward, self.done, entry

def main():
    download()
    cat = json.loads((DATA / 'catalog.json').read_text())['rows'][0]
    with fits.open(DATA / 'target_pixels.fits.gz') as hd:
        hdr = hd[0].header.copy()
        table = hd[1].data
        aperture = (hd[2].data & 2) != 0
        cube = np.array(table['FLUX'], dtype=float)
        errcube = np.array(table['FLUX_ERR'], dtype=float)
        t = np.array(table['TIME'], dtype=float)
        quality = np.array(table['QUALITY'])
        cadences = np.array(table['CADENCENO'])
        # FLUX is already background-subtracted. Do not subtract FLUX_BKG again.
        flux = np.sum(cube[:, aperture], axis=1)
        error = np.sqrt(np.sum(errcube[:, aperture] ** 2, axis=1))
    good = (quality == 0) & np.isfinite(t) & np.isfinite(flux) & np.isfinite(error) & (error > 0) & (flux > 0)
    image = np.full(cube.shape[1:], np.nan)
    observed = np.any(np.isfinite(cube[good]), axis=0)
    image[observed] = np.nanmedian(cube[good][:, observed], axis=0)
    t, f, e = t[good], flux[good], error[good]
    order = np.argsort(t)
    t, f, e = t[order], f[order], e[order]
    # Robust ~2-day baseline in each continuous segment, preserving deep dips.
    breaks = np.r_[0, np.flatnonzero(np.diff(t) > 0.3) + 1, len(t)]
    trend = np.empty_like(f)
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        window = min(101, hi-lo if (hi-lo) % 2 else hi-lo-1)
        trend[lo:hi] = median_filter(f[lo:hi], size=max(1, window), mode='nearest')
    y, dy = f / trend, e / trend
    np.savez_compressed(ROOT / 'extracted_lightcurve.npz', time_bkjd=t,
                        aperture_flux=f, flux_error=e, trend=trend,
                        normalized_flux=y, normalized_error=dy)
    fits.writeto(ROOT / 'median_image.fits', image, overwrite=True)
    with fits.open(DATA / 'pipeline_lightcurve.fits') as hd:
        lc = hd[1].data
        common, itpf, ilc = np.intersect1d(cadences, lc['CADENCENO'], return_indices=True)
        ok = good[itpf] & np.isfinite(lc['SAP_FLUX'][ilc])
        relative = np.abs(flux[itpf][ok] - lc['SAP_FLUX'][ilc][ok]) / np.abs(lc['SAP_FLUX'][ilc][ok])
    separation = float(SkyCoord(hdr['RA_OBJ'], hdr['DEC_OBJ'], unit='deg').separation(
        SkyCoord(cat['ra'], cat['dec'], unit='deg')).arcsec)
    features = {'temperature_K': cat['st_teff'], 'temperature_error_K': cat['st_tefferr1'],
                'metallicity_dex': cat['st_met'], 'metallicity_error_dex': cat['st_meterr1'],
                'metallicity_definition': cat['st_metratio'], 'radius_solar': cat['st_rad'],
                'source': 'NASA Exoplanet Archive ps default row; not a LAMOST retrieval',
                'reference_html': cat['st_refname']}
    def evaluator(candidate):
        p = cat['pl_orbper']
        perr = abs(candidate['period'] - p) / p
        dt = abs((candidate['transit_time'] - (cat['pl_tranmid'] - 2454833) + p/2) % p - p/2)
        return {'relative_period_error': perr, 'epoch_error_days': dt,
                'period_tolerance_fraction': 0.001, 'epoch_tolerance_days': 0.1,
                'period_and_epoch_match': bool(perr < 0.001 and dt < 0.1)}
    env = SearchEnvironment(t, y, dy, features, evaluator)
    initial_observation = env.observe()
    for action in ['coarse', 'refine', 'submit']:
        env.step(action)
    candidate = env.candidate
    evaluation = env.history[-1]['evaluation']
    # Actual state-machine checks: reject overspending and post-terminal actions.
    rejected_terminal = False
    try:
        env.step('submit')
    except ValueError:
        rejected_terminal = True
    low_budget = SearchEnvironment(t, y, dy, features, evaluator, budget=1)
    rejected_budget = False
    try:
        low_budget.step('coarse')
    except ValueError:
        rejected_budget = True
    checks = {'usable_cadences_over_1000': len(t) > 1000,
              'identity_separation_under_1_arcsec': separation < 1,
              'aperture_matches_official_SAP_within_1e_minus_5': float(np.max(relative)) < 1e-5,
              'period_and_epoch_recovered': evaluation['period_and_epoch_match'],
              'budget_exhausted_exactly': env.budget == 0,
              'reject_post_terminal_action': rejected_terminal,
              'reject_unaffordable_action': rejected_budget}
    result = {'run_utc': datetime.now(timezone.utc).isoformat(),
              'target': 'Kepler-8 / KIC 6922244', 'quarter': int(hdr['QUARTER']),
              'status': 'PASS' if all(checks.values()) else 'FAIL', 'checks': checks,
              'image_shape': list(image.shape), 'total_cadences': len(cube),
              'retained_cadences': len(t), 'aperture_pixels': int(aperture.sum()),
              'baseline_days': float(np.ptp(t)), 'coordinate_separation_arcsec': separation,
              'max_aperture_SAP_relative_difference': float(np.max(relative)),
              'stellar_parameters': features, 'candidate': candidate,
              'evaluation': evaluation, 'catalog_period_days': cat['pl_orbper'],
              'catalog_radius_earth': cat['pl_rade'],
              'initial_agent_observation': initial_observation,
              'episode': env.history, 'total_reward': sum(s['reward'] for s in env.history),
              'rl_training_performed': False, 'habitable_zone_assessment_performed': False,
              'limitations': ['Single known giant planet selected as a positive control.',
                  'Scripted actions, abstract cost units, and period/epoch-only reward.',
                  'No LAMOST/Gaia cross-match, Hubble ingestion, or TESS validation.',
                  'No false-positive vetting, injection recovery, learned ranking, or GRPO training.',
                  'Box depth and formal SNR are not physical planet fits or detection significance.'],
              'software': {m: v for m,v in [('python',platform.python_version()),('numpy',np.__version__),
                  ('scipy',scipy.__version__),('astropy',astropy.__version__),('matplotlib',matplotlib.__version__)]},
              'data_manifest': [{'file': p.name, 'bytes': p.stat().st_size,
                  'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                  'url': URLS.get(p.name, json.loads((DATA/'catalog.json').read_text())['request_url'])}
                  for p in sorted(DATA.iterdir()) if p.is_file()]}
    save_json('results.json', result)
    phase = (t - candidate['transit_time'] + candidate['period']/2) % candidate['period'] - candidate['period']/2
    fig, axes = plt.subplots(2,2,figsize=(12,8),layout='constrained')
    ax = axes[0,0]
    im = ax.imshow(image,origin='lower',cmap='magma')
    for row,col in np.argwhere(aperture):
        ax.add_patch(plt.Rectangle((col-.5,row-.5),1,1,fill=False,edgecolor='cyan',lw=1.5))
    fig.colorbar(im,ax=ax,label='Median flux (electrons / second)')
    ax.set(title='Downloaded Kepler pixels; cyan = aperture',xlabel='Cutout column',ylabel='Cutout row')
    axes[0,1].plot(t,y,'.',ms=2,color='#156082')
    axes[0,1].set(title='Light curve extracted from the images',xlabel='Time (BJD - 2454833)',ylabel='Normalized flux')
    pg = env.periodogram
    axes[1,0].plot(pg.period,pg.power,color='#156082',lw=.8)
    axes[1,0].axvline(cat['pl_orbper'],color='#d35a26',ls='--',label='Catalog period (evaluation only)')
    axes[1,0].set(title='Blind box-least-squares search: 0.5–10 days',xlabel='Trial period (days)',ylabel='BLS power')
    axes[1,0].legend(fontsize=8)
    axes[1,1].plot(phase*24,y,'.',ms=2,alpha=.5,color='#156082')
    axes[1,1].set(xlim=(-8,8),title=f"Recovered period: {candidate['period']:.6f} days",xlabel='Hours from fitted mid-transit',ylabel='Normalized flux')
    fig.suptitle('Kepler-8 pipeline smoke test — known giant planet, no RL training',fontsize=14)
    fig.savefig(ROOT/'diagnostic.png',dpi=160)
    plt.close(fig)
    print(json.dumps({'status':result['status'],'checks':checks,'candidate':candidate,
                      'evaluation':evaluation,'stellar_parameters':features},indent=2))
    assert all(checks.values()), 'One or more smoke-test checks failed; inspect results.json'

if __name__ == '__main__':
    main()
