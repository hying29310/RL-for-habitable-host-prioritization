"""Verify existing HST detections against all Gaia-linked stars in the cached LAMOST cone.
Select one observation per star: prefer zero fiber mask, then highest g-band S/N.
Preserve flagged observations when no clean one exists. No new detections are made.
"""
from pathlib import Path
import json
import numpy as np
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord,Distance
from astropy.time import Time
from astropy.io import fits
import astropy.units as u
from retrieve_matched_stars import retrieve
ROOT=Path(__file__).resolve().parent

def main():
    sources=json.loads((ROOT/'sources.json').read_text())
    sky=SkyCoord([r['ra_deg'] for r in sources]*u.deg,[r['dec_deg'] for r in sources]*u.deg)
    cat=parse_single_table(ROOT/'data/lamost_cone.xml').to_table(use_names_over_ids=True)
    gj=json.loads((ROOT/'data/gaia.json').read_text())
    gaia={str(r[0]):dict(zip([m['name'] for m in gj['metadata']],r)) for r in gj['data']}
    with fits.open(ROOT/'data/icmp09011_drz.fits') as h:epoch=Time(h[0].header['EXPSTART'],format='mjd')
    matches=[];audit=[]
    for gid in sorted(set(str(v) for v in cat['gaia_source_id'])):
        if gid not in gaia:continue
        g=gaia[gid]
        c=SkyCoord(ra=g['ra']*u.deg,dec=g['dec']*u.deg,pm_ra_cosdec=g['pmra']*u.mas/u.yr,
            pm_dec=g['pmdec']*u.mas/u.yr,distance=Distance(parallax=g['parallax']*u.mas),
            radial_velocity=0*u.km/u.s,obstime=Time(g['ref_epoch'],format='jyear'))
        predicted=c.apply_space_motion(new_obstime=epoch)
        sep=sky.separation(predicted).arcsec;order=np.argsort(sep);i=int(order[0])
        rows=cat[np.asarray([str(v)==gid for v in cat['gaia_source_id']])]
        # Finite parameters and S/N are required for this particular extraction.
        usable=[r for r in rows if not any(np.ma.is_masked(r[k]) for k in ['teff','feh','snrg'])
            and np.isfinite([r['teff'],r['feh'],r['snrg']]).all() and r['teff']>0 and r['feh']>-5]
        usable.sort(key=lambda r:(int(r['fibermask'])!=0,-float(r['snrg'])))
        accepted=bool(sep[i]<0.3 and np.sum(sep<0.3)==1 and len(usable)>0)
        a={'gaia_source_id':gid,'lamost_observation_count':len(rows),
           'closest_hst_source_id':sources[i]['source_id'],'separation_arcsec':float(sep[i]),
           'second_nearest_arcsec':float(sep[order[1]]),'positional_match_accepted':accepted}
        audit.append(a)
        if not accepted:continue
        row=usable[0];mask=int(row['fibermask'])
        matches.append({'source_id':sources[i]['source_id'],'lamost_obsid':int(row['obsid']),
            'lamost_release':'DR11 v2.0 LRS','gaia_dr3_source_id':gid,
            'match_status':'positionally_validated','ra_deg':sources[i]['ra_deg'],'dec_deg':sources[i]['dec_deg'],
            'epoch_corrected_match_arcsec':float(sep[i]),'fibermask':mask,
            'spectrum_quality_status':'zero_fiber_mask' if mask==0 else 'flagged_requires_review',
            'fiber_flag_description':'WHOPPER (bit value 64)' if mask==64 else ('none' if mask==0 else 'see LAMOST documentation'),
            'available_observations':len(rows),'selection_rule':'prefer zero fibermask, then highest snrg'})
    assert len({r['source_id'] for r in matches})==len(matches),'Multiple catalog stars matched to one HST detection'
    (ROOT/'all_three_way_matches.json').write_text(json.dumps(matches,indent=2)+'\n')
    (ROOT/'three_way_match_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    manifest=retrieve(ROOT/'all_three_way_matches.json',ROOT/'data/lamost_cone.xml',ROOT/'data',ROOT/'retrieved_three_way')
    print(json.dumps({'distinct_stars':len(matches),'matches':matches,'retrieval':manifest},indent=2))
    assert not manifest['failed']

if __name__=='__main__':main()
