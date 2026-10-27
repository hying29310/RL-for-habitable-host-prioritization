"""Verify final counts and make a plain-text report and Hubble/spectrum previews."""
from pathlib import Path
import os,json,math
from collections import Counter
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import simple_norm

def write(p,o):p.write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')

def main():
    matches=json.loads((ROOT/'matches.json').read_text());props=json.loads((ROOT/'retrieved/stellar_properties.json').read_text())
    byid={x['gaia_dr3_source_id']:x for x in props};cuts={x['gaia_dr3_source_id']:x for x in json.loads((ROOT/'cutout_manifest.json').read_text())}
    gids={m['gaia_dr3_source_id'] for m in matches}
    assert len(gids)==len(matches) and len({m['hsc_match_id'] for m in matches})==len(matches)
    assert len(matches)>=30,'Do not report success below the requested sample size'
    assert gids==set(byid),'Spectral retrieval and final accepted matches must agree'
    assert gids<=set(cuts),'Missing image cutouts'
    records=[]
    for i,m in enumerate(matches,1):
        prop=byid[m['gaia_dr3_source_id']]
        assert prop['lamost_obsid']==m['lamost_obsid']
        assert prop['quality']['stellar_parameters_available']
        assert prop['spectrum']['sample_count']>0
        assert m['neighbor_check']['passed'] and m['lamost_gaia_separation_arcsec']<1.5
        records.append({'number':i,**m,'stellar_properties':prop,'hst_cutout':cuts[m['gaia_dr3_source_id']]})
    write(ROOT/'star_catalog.json',records)
    oldp=ROOT/'previous_M67_matches.json'
    old=json.loads(oldp.read_text());assert not(gids & {m['gaia_dr3_source_id'] for m in old})
    combined=[dict(m,field='M67',method='previous direct-image pipeline') for m in old]+matches
    write(ROOT/'all_fields_including_M67.json',combined)
    counts=dict(sorted(Counter(m['field'] for m in matches).items()))
    summary={'outside_M67':len(matches),'M67_previous':len(old),'total_distinct_stars':len(combined),'outside_M67_fields':counts,
        'spectra_downloaded':len(props),'hst_cutouts_downloaded':len(gids & set(cuts)),
        'lamost_zero_fiber_mask':sum(m['fibermask']==0 for m in matches),
        'lamost_flagged':sum(m['fibermask']!=0 for m in matches),
        'has_HSC_zero_flag_detection':sum(m['hst_has_zero_flag_detection'] for m in matches),
        'both_lamost_and_at_least_one_HSC_detection_zero_flag':sum(m['fibermask']==0 and m['hst_has_zero_flag_detection'] for m in matches),
        'multiple_independent_white_light_detections':sum(m['distinct_white_light_detections']>=2 for m in matches),
        'method':'HSC v3 image detections crossmatched to Gaia DR3 with epoch correction, then checked against LAMOST DR11 v2.0 positions and observation IDs.'}
    write(ROOT/'summary.json',summary)
    # Plot real downloaded FITS pixels. A circle marks the requested cutout center,
    # not a newly fitted centroid; the catalog measurement is recorded in JSON.
    ncols=5;nrows=math.ceil(len(matches)/ncols)
    fig,axs=plt.subplots(nrows,ncols,figsize=(15,nrows*3),layout='constrained')
    for ax in axs.flat:ax.axis('off')
    for ax,r in zip(axs.flat,records):
        with fits.open(r['hst_cutout']['path']) as h:a=next(x.data for x in h if x.data is not None and x.data.ndim==2).astype(float)
        ax.imshow(a,origin='lower',cmap='gray',norm=simple_norm(a,'asinh',percent=99.5));ax.add_patch(plt.Circle(((a.shape[1]-1)/2,(a.shape[0]-1)/2),6,edgecolor='cyan',fill=False,lw=1))
        ax.set_title(f"{r['number']}. {r['field']} | Gaia ...{r['gaia_dr3_source_id'][-6:]}\nLAMOST flag {r['fibermask']} | HSC flag {r['best_hst_detection']['Flags']}",fontsize=9)
    fig.suptitle(f'{len(matches)} stars outside M67: actual Hubble image cutouts\nCyan circles mark requested positions; flags are retained for review.',fontsize=16)
    fig.savefig(ROOT/'hubble_contact_sheet.png',dpi=120);plt.close(fig)
    fig,axs=plt.subplots(nrows,ncols,figsize=(15,nrows*2.2),layout='constrained')
    for ax in axs.flat:ax.axis('off')
    for ax,r in zip(axs.flat,records):
        ax.axis('on');s=r['stellar_properties'];npz=np.load(s['spectrum']['arrays_file'])
        wave=npz['wavelength_angstrom'];flux=npz['catalog_flux'];good=np.isfinite(flux)
        ax.plot(wave[good],flux[good],lw=.35);lo,hi=np.percentile(flux[good],[1,99]);ax.set_ylim(lo,hi)
        ax.set_title(f"{r['number']}. {r['field']} | {s['temperature_K']:.0f} K\n[Fe/H] {s['Fe_H_dex']:+.2f}, fiber flag {s['fibermask']}",fontsize=9);ax.tick_params(labelsize=6)
    fig.suptitle('Downloaded LAMOST spectra (original flux scale; wavelength in Angstrom)',fontsize=14)
    fig.savefig(ROOT/'spectra_contact_sheet.png',dpi=120);plt.close(fig)
    lines=['MATCHED STAR SURVEY — RESULT',f"{len(matches)} distinct matched stars outside M67; {len(combined)} including the previous M67 sample.",
           '','Field counts:']+[f'{f}: {n}' for f,n in counts.items()]
    lines += ['',f"Retrieved {len(props)} LAMOST spectra and {len(gids & set(cuts))} Hubble FITS cutouts.",
              f"LAMOST fiber mask zero: {summary['lamost_zero_fiber_mask']}; flagged: {summary['lamost_flagged']}.",
              f"At least one zero-flag Hubble detection: {summary['has_HSC_zero_flag_detection']}.",
              f"Zero LAMOST mask AND at least one zero-flag HSC detection: {summary['both_lamost_and_at_least_one_HSC_detection_zero_flag']}.",
              'Zero flags do not prove perfect measurements or an unblended spectrum.',
              '', 'One row per unique Gaia star. Repeated spectra/images do not increase the count.',
              'These are catalog associations, not detected planets, confirmed habitability,',
              'or confirmed cluster membership. See README.txt for methods and limitations.','',
              'NO. FIELD      GAIA DR3 SOURCE ID    LAMOST OBSID   TEFF(K)  [Fe/H]  FIBER MASK']
    for r in records:
        s=r['stellar_properties'];lines.append(f"{r['number']:3d} {r['field']:10s} {r['gaia_dr3_source_id']:19s} {r['lamost_obsid']:12d} {s['temperature_K']:8.1f} {s['Fe_H_dex']:7.3f} {s['fibermask']:5d}")
    (ROOT/'RESULTS.txt').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
