from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import requests,json
from astropy.io.votable import parse_single_table
ROOT=Path(__file__).resolve().parent/'data'
ROOT.mkdir(parents=True,exist_ok=True)
fields={'NGC6791':(290.22,37.77),'NGC6819':(295.33,40.187),'NGC2420':(114.535,21.575),'NGC188':(11.80,85.25)}
def job(name,coord,kind):
 p=ROOT/name;p.mkdir(exist_ok=True)
 ra,dec=coord
 if kind=='lamost':
  url='https://www.lamost.org/dr11/v2.0/voservice/conesearch';params={'RA':ra,'DEC':dec,'SR':0.08}
  r=requests.get(url,params=params,timeout=90);r.raise_for_status();(p/'lamost.xml').write_text(r.text)
  t=parse_single_table(p/'lamost.xml').to_table(use_names_over_ids=True)
  good=[a for a in t if not any(__import__('numpy').ma.is_masked(a[k]) for k in ['teff','feh','gaia_source_id']) and a['teff']>0]
  return name,kind,{'rows':len(t),'parameter_rows':len(good),'unique_gaia':len(set(str(a['gaia_source_id']) for a in good))}
 req={'service':'Mast.Caom.Cone','params':{'ra':ra,'dec':dec,'radius':0.08},'format':'json','pagesize':5000}
 r=requests.post('https://mast.stsci.edu/api/v0/invoke',data={'request':json.dumps(req)},timeout=90);r.raise_for_status();j=r.json();(p/'mast.json').write_text(r.text)
 h=[a for a in j.get('data',[]) if a.get('obs_collection')=='HST' and a.get('dataproduct_type')=='image' and a.get('instrument_name') in ['ACS/WFC','WFC3/UVIS']]
 return name,kind,{'rows':len(j.get('data',[])),'hst_imaging':len(h),'examples':[{k:a.get(k) for k in ['obsid','obs_id','filters','t_exptime']} for a in sorted(h,key=lambda a:a.get('t_exptime') or 1e10)[:5]]}
with ThreadPoolExecutor(max_workers=4) as ex:
 fs=[ex.submit(job,n,c,k) for n,c in fields.items() for k in ['lamost','mast']]
 for f in as_completed(fs):
  try:print(json.dumps(f.result()),flush=True)
  except Exception as e:print('ERROR',repr(e),flush=True)
(ROOT/'search_centers.json').write_text(json.dumps(fields,indent=2))
