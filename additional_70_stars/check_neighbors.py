"""Reject associations if another Gaia star fits the selected Hubble detection.
Save local-neighbor query results so the check can be audited and repeated offline.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.io.votable import parse_single_table
from validate_matches import ROOT,DATA,fetch,position_at,save

def check(m, lamost):
    gid=m['gaia_dr3_source_id'];d=m['best_hst_detection'];p=DATA/'gaia_neighbors';p.mkdir(exist_ok=True)
    # 5-arcsec catalog-epoch neighborhood: candidate selection, not matching tolerance.
    q=f"SELECT source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,ruwe,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{m['ra_deg']},{m['dec_deg']},{5/3600}))=1"
    j=json.loads(fetch('https://gea.esac.esa.int/tap-server/tap/sync',p/f'{gid}.json',params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':q}))
    cols=[c['name'] for c in j['metadata']];rows=[dict(zip(cols,r)) for r in j['data']]
    ref=SkyCoord(d['SourceRA']*u.deg,d['SourceDec']*u.deg);dist=[]
    for r in rows:
        has_pm=r['pmra'] is not None and r['pmdec'] is not None
        co=position_at(r,(d['StartMJD']+d['StopMJD'])/2) if has_pm else SkyCoord(r['ra']*u.deg,r['dec']*u.deg)
        dist.append({'gaia_source_id':str(r['source_id']),'separation_arcsec':float(ref.separation(co).arcsec),'proper_motion_available':has_pm})
    dist.sort(key=lambda r:r['separation_arcsec'])
    others=[x for x in dist if x['gaia_source_id']!=gid]
    lcoord=SkyCoord(float(lamost['ra'])*u.deg,float(lamost['dec'])*u.deg)
    ldist=[]
    for r in rows:
        co=position_at(r,Time(str(lamost['obsdate'])).mjd) if r['pmra'] is not None and r['pmdec'] is not None else SkyCoord(r['ra']*u.deg,r['dec']*u.deg)
        ldist.append((float(lcoord.separation(co).arcsec),str(r['source_id'])))
    ldist.sort()
    lsep=next((sep for sep,rid in ldist if rid==gid),None)
    lpass=bool(ldist) and ldist[0][1]==gid and lsep is not None and lsep<1.5
    passed=bool(dist) and dist[0]['gaia_source_id']==gid and (not others or others[0]['separation_arcsec']>0.6) and lpass
    m['lamost_gaia_separation_arcsec']=lsep
    return m,{'gaia_dr3_source_id':gid,'passed':passed,'lamost_gaia_separation_arcsec':lsep,'lamost_nearest_gaia_source_id':ldist[0][1] if ldist else None,'lamost_position_check_passed':lpass,'neighbors':dist,'rule':'Selected Gaia source must be nearest; no other Gaia source within 0.6 arcsec of selected HSC detection. Missing-PM neighbors use catalog position. Gaia must also be nearest to LAMOST coordinates at its observing date and within 1.5 arcsec.'}

def main():
    ms=json.loads((ROOT/'preliminary_matches.json').read_text());out=[];audit=[]
    table=parse_single_table(ROOT/'selected_lamost.xml').to_table(use_names_over_ids=True)
    byobs={int(r['obsid']):r for r in table}
    with ThreadPoolExecutor(4) as pool:
        fs={pool.submit(check,m,byobs[m['lamost_obsid']]):m for m in ms}
        for f in as_completed(fs):
            try:
                m,a=f.result();audit.append(a)
                if a['passed']:
                    m['neighbor_check']=a;out.append(m)
            except Exception as e:audit.append({'gaia_dr3_source_id':fs[f]['gaia_dr3_source_id'],'passed':False,'error':str(e)})
    out.sort(key=lambda m:(m['field'],m['gaia_dr3_source_id']))
    save(ROOT/'matches.json',out);save(ROOT/'neighbor_audit.json',audit)
    print('ACCEPTED',len(out),'REJECTED/FAILED',len(audit)-len(out),flush=True)
if __name__=='__main__':main()
