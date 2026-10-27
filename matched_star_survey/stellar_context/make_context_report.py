"""Summarize actual Gaia availability and downloaded light-curve segments."""
from pathlib import Path
from collections import Counter
import json,os,math
import numpy as np
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    stars=json.loads((ROOT/'availability.json').read_text());manifest=json.loads((ROOT/'lightcurve_manifest.json').read_text());byid={r['gaia_dr3_source_id']:r for r in manifest}
    out=[];lines=['RADIUS, LUMINOSITY AND LIGHT-CURVE AVAILABILITY','',
        '31 Hubble--LAMOST--Gaia matched stars checked.',
        '24 have Gaia DR3 FLAME radius AND luminosity estimates.',
        'One additional star has a Gaia GSP-Phot radius estimate (25 radii total).',
        '23 have a downloaded, readable light-curve segment: 21 Kepler + 2 TESS QLP.',
        '21 stars have both FLAME properties and a downloaded light curve.',
        '', 'FIELD       GAIA DR3 SOURCE ID    RADIUS(Rsun)  LUMINOSITY(Lsun)  LIGHT CURVE']
    for s in stars:
        gid=s['gaia_dr3_source_id'];g=s['gaia_parameters'];lc=byid[gid];r=g.get('radius_flame');l=g.get('lum_flame');source='Gaia DR3 FLAME'
        if r is None:
            r=g.get('radius_gspphot');source='Gaia DR3 GSP-Phot' if r is not None else None
        low=g.get('radius_flame_lower') if source=='Gaia DR3 FLAME' else g.get('radius_gspphot_lower')
        high=g.get('radius_flame_upper') if source=='Gaia DR3 FLAME' else g.get('radius_gspphot_upper')
        status=lc.get('observation',{}).get('provenance_name','not obtained') if lc['status']=='downloaded' else 'not obtained'
        out.append({'gaia_dr3_source_id':gid,'field':s['field'],'radius_solar':r,'radius_source':source,'radius_16th_percentile_solar':low,'radius_84th_percentile_solar':high,'luminosity_solar':l,'luminosity_source':'Gaia DR3 FLAME' if l is not None else None,'luminosity_16th_percentile_solar':g.get('lum_flame_lower'),'luminosity_84th_percentile_solar':g.get('lum_flame_upper'),'gaia_raw_parameters':g,'lightcurve':lc})
        rr='missing' if r is None else f'{r:.3f}';ll='missing' if l is None else f'{l:.3f}'
        lines.append(f"{s['field']:11s} {gid:19s} {rr:>12s} {ll:>17s}  {status}")
    complete=sum(x['radius_solar'] is not None and x['luminosity_solar'] is not None and x['lightcurve']['status']=='downloaded' for x in out)
    assert len(out)==31 and complete==21
    (ROOT/'stellar_properties_and_lightcurves.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    lines+=['','HOW TO READ THIS','Rsun and Lsun mean times the Sun\'s radius and luminosity.',
        'Gaia values are model-dependent catalog estimates. Lower/upper entries are',
        '16th/84th percentiles; systematic uncertainty is not fully represented.',
        'GSP-Phot and FLAME are distinct estimates; do not mix their errors.',
        'No luminosity was invented for the star with only a GSP-Phot radius.',
        '', 'LIGHT-CURVE LIMITS',
        'The initial MAST search found nearby time-series records for 24 stars.',
        '23 passed the selection for a standard downloadable Kepler/QLP product.',
        'Selection required a target coordinate within 1 arcsec of Gaia, with our',
        'Gaia star the nearest of the cached neighboring Gaia sources. The FITS',
        'Kepler/TIC target ID was checked against the selected archive record.',
        'This is a positional association, not a complete cross-catalog identity audit.',
        'Positions were compared at the Gaia reference epoch; proper-motion differences',
        'between archive input catalogs were not fully corrected in this extra step.',
        'Only one segment per star was downloaded. Kepler segments are individual',
        'quarters; TESS segments are individual sectors. This is not the full baseline.',
        'Time, flux, errors where available, and quality flags are saved in NPZ files.',
        'Usable preview samples require finite time/flux and quality==0.',
        'One QLP file uses SAP_FLUX with no corresponding error column; its exported',
        'errors are NaN. No uncertainty was substituted from a different flux column.',
        'No transit search, aperture-blending assessment or planet validation was done.',
        'No nearby result does not prove no telescope coverage. A broader catalog',
        'search or custom full-frame-image extraction may recover more light curves.',
        '', 'FILES AND RUNNING',
        'check_availability.py: queries Gaia AP and MAST time-series observations.',
        'download_lightcurves.py: selects and downloads one light-curve segment per star.',
        'make_context_report.py: creates this report, unified JSON and preview.',
        'Run these scripts sequentially with the project work/venv/bin/python.',
        'lightcurve_manifest.json: file paths, target IDs, data columns, times and warnings.',
        'gaia_properties.json: unmodified mapped Gaia parameters, including FLAME flags.',
        'data/: cached archive responses and original FITS light curves.',
        'lightcurves/: exported NPZ arrays.',
        '', 'REFERENCES',
        'Gaia DR3 stellar-parameter methods: https://www.aanda.org/articles/aa/pdf/2023/06/aa43919-22.pdf',
        'MAST API: https://mast.stsci.edu/api/v0/MastApiTutorial.html',
        'Gaia TAP schema descriptions are saved in data/gaia_schema.json.']
    (ROOT/'RESULTS.txt').write_text('\n'.join(lines)+'\n')
    good=[x for x in out if x['lightcurve']['status']=='downloaded'];nc=4;nr=math.ceil(len(good)/nc)
    fig,axs=plt.subplots(nr,nc,figsize=(14,nr*2.4),layout='constrained')
    for ax in axs.flat:ax.axis('off')
    for ax,x in zip(axs.flat,good):
        a=np.load(x['lightcurve']['arrays_file']);mask=a['usable'];t=a['time_days'][mask];f=a['flux'][mask];median=np.median(f)
        ax.axis('on');ax.scatter(t-t.min(),f/median,s=.7,rasterized=True);ax.set_title(f"{x['field']} | Gaia ...{x['gaia_dr3_source_id'][-6:]}\n{x['lightcurve']['observation']['provenance_name']} / {x['lightcurve']['flux_column']}",fontsize=8);ax.tick_params(labelsize=7);ax.set_xlabel('Days since first usable sample',fontsize=7);ax.set_ylabel('Flux / median',fontsize=7)
    fig.suptitle('23 downloaded light-curve segments — preview only\nQuality-zero points; no detrending, blending checks or transit search applied.',fontsize=13)
    fig.savefig(ROOT/'lightcurve_preview.png',dpi=130);plt.close(fig)
    final={'stars_checked':31,'radius_estimates':sum(x['radius_solar'] is not None for x in out),'luminosity_estimates':sum(x['luminosity_solar'] is not None for x in out),'downloaded_lightcurves':len(good),'all_three_available':complete,'lightcurve_missions':dict(Counter(x['lightcurve']['observation']['obs_collection'] for x in good))}
    (ROOT/'final_summary.json').write_text(json.dumps(final,indent=2)+'\n');print(final)
if __name__=='__main__':main()
