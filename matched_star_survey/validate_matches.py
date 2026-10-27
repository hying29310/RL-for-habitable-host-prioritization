"""Validate LAMOST--Gaia--Hubble Source Catalog associations with epoch corrections.
Uses published HSC detections rather than rerunning DAOStarFinder on full images.
Every candidate and rejection is saved; repeated observations count as one star.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,warnings,hashlib
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io.votable import parse_single_table
from astropy.table import vstack,unique
import astropy.units as u
from search_fields import ROOT,DATA,fetch

def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def position_at(g,mjd):
    """Propagate Gaia proper motion to each Hubble observation date (no radial velocity)."""
    co=SkyCoord(ra=g['ra']*u.deg,dec=g['dec']*u.deg,pm_ra_cosdec=g['pmra']*u.mas/u.yr,
        pm_dec=g['pmdec']*u.mas/u.yr,obstime=Time(g['ref_epoch'],format='jyear'))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return co.apply_space_motion(new_obstime=Time(mjd,format='mjd'))

def gaia_sources(ids):
    """Fetch exact Gaia identifiers without converting 64-bit IDs to floats."""
    ids=sorted(ids);out={}
    for i in range(0,len(ids),250):
        chunk=ids[i:i+250];key=hashlib.sha256(','.join(chunk).encode()).hexdigest()[:16]
        q='SELECT source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,ruwe,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE source_id IN ('+','.join(chunk)+')'
        j=json.loads(fetch('https://gea.esac.esa.int/tap-server/tap/sync',DATA/f'gaia_{key}.json',params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':q}))
        cols=[x['name'] for x in j['metadata']]
        out.update({str(row[0]):dict(zip(cols,row)) for row in j['data']})
    return out

def details(mid):
    """Read individual HSC image detections, preserving image IDs, flags and dates."""
    p=DATA/'hsc_details';p.mkdir(exist_ok=True)
    j=json.loads(fetch('https://catalogs.mast.stsci.edu/api/v0.1/hsc/v3/detailed.json',p/f'{mid}.json',params={'MatchID':mid,'flatten_response':'false','pagesize':50000}))
    if len(j.get('data',[]))>=50000:raise ValueError('HSC response may be truncated')
    return j['data']

def main():
    paths=list(DATA.glob('*/hsc_candidates.json'))
    # Prefer the larger cone when both pilot and expanded results are available.
    paths=[p for p in paths if not(p.parent.name.endswith('_pilot') and (DATA/p.parent.name.removesuffix('_pilot')/'hsc_candidates.json').exists())]
    candidate_by_star={};tables={}
    for p in paths:
        targets=json.loads((p.parent/'targets.json').read_text());hs=json.loads(p.read_text())
        tables[p.parent.name]=parse_single_table(p.parent/'lamost.xml').to_table(use_names_over_ids=True)
        for h in hs:
            t=targets[h['_target_index']];gid=t['gaia_dr3_source_id']
            entry=candidate_by_star.setdefault(gid,{'target':t,'field':p.parent.name.removesuffix('_pilot').removesuffix('_core'),'catalog_field':p.parent.name,'hsc':{}})
            entry['hsc'][str(h['MatchID'])]=h
    gaia=gaia_sources(candidate_by_star);save(ROOT/'gaia_sources.json',gaia)
    jobs={};audit=[]
    for gid,e in candidate_by_star.items():
        g=gaia.get(gid)
        if not g or g['pmra'] is None or g['pmdec'] is None:continue
        for mid,h in e['hsc'].items():
            if h['NumImages']<2 or h['AbsCorr']!='Y':continue
            pred=position_at(g,(h['StartMJD']+h['StopMJD'])/2)
            sep=pred.separation(SkyCoord(h['MatchRA']*u.deg,h['MatchDec']*u.deg)).arcsec
            if sep<0.6:jobs[mid]=None
    print('Candidates:',len(candidate_by_star),'HSC detail queries:',len(jobs),flush=True)
    with ThreadPoolExecutor(4) as pool:
        futures={pool.submit(details,mid):mid for mid in jobs}
        for f in as_completed(futures):jobs[futures[f]]=f.result()
    matches=[];selected_rows=[]
    for gid,e in candidate_by_star.items():
        g=gaia.get(gid);a={'gaia_dr3_source_id':gid,'field':e['field'],'candidate_hsc_matches':len(e['hsc']),'evaluations':[]}
        if not g or g['pmra'] is None or g['pmdec'] is None:
            a['status']='missing_gaia_astrometry';audit.append(a);continue
        accepted=[]
        for mid,h in e['hsc'].items():
            if mid not in jobs:continue
            rows=[r for r in jobs[mid] if r['Det']=='Y' and r['SourceRA'] is not None and r['StartMJD'] is not None]
            if not rows:continue
            times=[(r['StartMJD']+r['StopMJD'])/2 for r in rows]
            pred=position_at(g,times);co=SkyCoord([r['SourceRA'] for r in rows]*u.deg,[r['SourceDec'] for r in rows]*u.deg)
            seps=co.separation(pred).arcsec
            good=[]
            for r,s,pr in zip(rows,seps,pred):
                if s<0.3:
                    d={k:r[k] for k in ['CatID','SourceID','ImageID','ImageName','SourceRA','SourceDec','XImage','YImage','StartMJD','StopMJD','Filter','Detector','Flags','ExposureTime','MagAper2']}
                    d.update(separation_arcsec=float(s),gaia_ra_at_hst_epoch=float(pr.ra.deg),gaia_dec_at_hst_epoch=float(pr.dec.deg));good.append(d)
            nimages=len({r['ImageName'] for r in good});nsource=len({r['SourceID'] for r in good});med=float(np.median(seps))
            a['evaluations'].append({'hsc_match_id':mid,'images_within_0p3_arcsec':nimages,'distinct_white_light_detections':nsource,'median_separation_arcsec':med})
            # At least two filter images and a consistent median; preserve independent-detection count.
            if nimages>=2 and med<0.3:accepted.append((h,good,med,nsource))
        if len(accepted)!=1:
            a['status']='ambiguous_hsc_association' if len(accepted)>1 else 'insufficient_consistent_hst_detections'
            audit.append(a);continue
        h,good,med,nsource=accepted[0];t=e['target'];flags=sorted({r['Flags'] for r in good})
        best=min(good,key=lambda r:(r['Flags']!=0, bool(r['Flags']&6),r['separation_arcsec']))
        m={'source_id':'Gaia_DR3_'+gid,'field':e['field'],'gaia_dr3_source_id':gid,'lamost_obsid':t['lamost_obsid'],
           'lamost_release':'DR11 v2.0 LRS','catalog_field':e['catalog_field'],
           'match_status':'epoch_checked_HSC_catalog_match','hsc_release':'v3','hsc_match_id':str(h['MatchID']),
           'median_epoch_corrected_separation_arcsec':med,'hst_images_within_0p3_arcsec':len({r['ImageName'] for r in good}),
           'distinct_white_light_detections':nsource,'hst_detection_flags':flags,'hst_has_zero_flag_detection':0 in flags,
           'fibermask':t['fibermask'],'gaia_ruwe':g['ruwe'],'gaia_g_mag':g['phot_g_mean_mag'],
           'ra_deg':g['ra'],'dec_deg':g['dec'],'best_hst_detection':best,'hst_detections':good,
           'caution':'Catalog positional association; HST photometry may be saturated/extended. Cluster membership and planets not established.'}
        matches.append(m);a['status']='accepted_pre_neighbor_check';a['hsc_match_id']=str(h['MatchID']);audit.append(a)
        table=tables[e['catalog_field']];selected_rows.append(table[np.asarray(table['obsid'])==t['lamost_obsid']])
    # An HSC source cannot be assigned to two different Gaia stars.
    mid_counts={m['hsc_match_id']:sum(x['hsc_match_id']==m['hsc_match_id'] for x in matches) for m in matches}
    matches=[m for m in matches if mid_counts[m['hsc_match_id']]==1]
    save(ROOT/'preliminary_matches.json',matches);save(ROOT/'matching_audit.json',audit)
    if selected_rows:
        unique(vstack(selected_rows,metadata_conflicts='silent'),keys='obsid').write(ROOT/'selected_lamost.xml',format='votable',overwrite=True)
    print('PRELIMINARY UNIQUE MATCHES',len(matches),flush=True)
    print({f:sum(m['field']==f for m in matches) for f in sorted({m['field'] for m in matches})},flush=True)
if __name__=='__main__':main()
