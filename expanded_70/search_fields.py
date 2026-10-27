"""Search public LAMOST catalogs and Hubble Source Catalog v3; save raw evidence.
These are preliminary candidates, not validated matches. Run validate_matches.py next.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,time,requests,numpy as np
from astropy.io.votable import parse_single_table
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data';DATA.mkdir(exist_ok=True)
# Approximate search centers in ICRS degrees. Labels do not assert membership.
FIELDS={'NGC6791':(290.22,37.77),'NGC6819':(295.33,40.187),'NGC2420':(114.535,21.575),
'NGC2158':(91.86,24.10),'M35':(92.28,24.33),'NGC7789':(359.34,56.72),
'NGC6633':(276.85,6.61),'NGC752':(29.42,37.78),'M13':(250.423,36.461),
'M3':(205.548,28.377),'M92':(259.281,43.136),'M15':(322.493,12.167),
'M2':(323.363,-0.823),'M53':(198.230,18.169),'NGC5466':(211.364,28.534),
'NGC5053':(199.113,17.698),'NGC4147':(182.526,18.543),'NGC2419':(114.535,38.882),
'NGC7006':(315.372,16.187),'NGC6229':(251.745,47.527),'NGC2264':(100.27,9.88),
'Orion':(83.82,-5.39),'M31':(10.684,41.269),'Pleiades':(56.75,24.117)}

def fetch(url,path,params=None,data=None):
    """Cache successful archive responses atomically; retry transient failures."""
    if path.exists():return path.read_text()
    for attempt in range(3):
        try:
            r=requests.post(url,json=data,params=params,timeout=(30,180)) if data is not None else requests.get(url,params=params,timeout=(30,180))
            r.raise_for_status()
            tmp=path.with_suffix('.partial');tmp.write_text(r.text);tmp.replace(path)
            return r.text
        except Exception:
            if attempt==2:raise
            time.sleep(2*(attempt+1))

def eligible_rows(table):
    """Keep distinct Gaia-linked stars with usable catalog parameters and g SNR>20."""
    by_id={}
    for row in table:
        if any(np.ma.is_masked(row[k]) for k in ['gaia_source_id','teff','feh','logg','snrg','fibermask']):continue
        gid=str(row['gaia_source_id'])
        if not gid.isdigit() or int(gid)<=0:continue
        if not np.isfinite([row['teff'],row['feh'],row['logg'],row['snrg']]).all():continue
        if row['teff']<=0 or row['feh']<=-5 or row['snrg']<=20 or str(row['class'])!='STAR':continue
        by_id.setdefault(gid,[]).append(row)
    return [sorted(rows,key=lambda r:(int(r['fibermask'])!=0,-float(r['snrg'])))[0] for rows in by_id.values()]

def search_field(item, radius=0.25):
    """Crossmatch a field's eligible LAMOST targets with HSC within 3 arcsec for discovery."""
    name,(ra,dec)=item;p=DATA/name;p.mkdir(exist_ok=True)
    fetch('https://www.lamost.org/dr11/v2.0/voservice/conesearch',p/'lamost.xml',params={'RA':ra,'DEC':dec,'SR':radius})
    t=parse_single_table(p/'lamost.xml').to_table(use_names_over_ids=True);rows=eligible_rows(t)
    targets=[{'ra':float(r['ra']),'dec':float(r['dec'])} for r in rows]
    candidates=[{'gaia_dr3_source_id':str(r['gaia_source_id']),'lamost_obsid':int(r['obsid']),'ra':float(r['ra']),'dec':float(r['dec']),'fibermask':int(r['fibermask'])} for r in rows]
    (p/'targets.json').write_text(json.dumps(candidates,indent=2))
    hits=[]
    if targets:
        for i in range(0,len(targets),500):
            raw=fetch('https://catalogs.mast.stsci.edu/api/v0.1/hsc/v3/detailed/crossmatch/json',p/f'hsc_candidates_{i}.json',params={'flatten_response':'false'},data={'targets':targets[i:i+500],'radius':0.000833,'resolve':False})
            j=json.loads(raw)
            for h in j.get('data',[]):
                h['_target_index']=int(h['_searchID_'])+i;hits.append(h)
    (p/'hsc_candidates.json').write_text(json.dumps(hits,indent=2))
    stats={'field':name,'lamost_observations':len(t),'eligible_stars':len(rows),'hsc_candidate_records':len(hits),'stars_with_hsc_candidates':len(set(h['_target_index'] for h in hits))}
    (p/'search_summary.json').write_text(json.dumps(stats,indent=2));return stats

def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wide',action='store_true',help='Also query 0.25-degree cones for all 24 exploratory fields')
    args=parser.parse_args()
    jobs=[(n,FIELDS[n],0.25) for n in ['NGC6791','NGC6819','NGC2420']]
    jobs += [(n+'_core',FIELDS[n],0.08) for n in ['M13','M3','M92','M15','NGC2158','NGC2264','Orion','NGC5466','M53','NGC6229']]
    if args.wide:jobs += [(n,c,0.25) for n,c in FIELDS.items() if n not in ['NGC6791','NGC6819','NGC2420']]
    (ROOT/'search_plan.json').write_text(json.dumps([{'field':n,'ra':c[0],'dec':c[1],'radius_deg':r} for n,c,r in jobs],indent=2))
    results=[]
    with ThreadPoolExecutor(max_workers=3) as pool:
        fs={pool.submit(search_field,(n,c),r):n for n,c,r in jobs}
        for f in as_completed(fs):
            try:s=f.result()
            except Exception as e:s={'field':fs[f],'error':str(e)}
            results.append(s);print(json.dumps(s),flush=True)
    (ROOT/'search_summary.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
