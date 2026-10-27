"""Download the selected HST images, LAMOST catalogs and Gaia astrometry.
Existing cache files are reused. Images total approximately 434 MB.
"""
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from astropy.io.votable import parse_single_table
import requests,json,numpy as np
ROOT=Path(__file__).resolve().parent/'data'
ROOT.mkdir(parents=True,exist_ok=True)
centers={'NGC6791':(290.22,37.77),'NGC6819':(295.33,40.187)}
for name,(ra,dec) in centers.items():
 p=ROOT/name;p.mkdir(exist_ok=True)
 if not (p/'lamost.xml').exists():
  r=requests.get('https://www.lamost.org/dr11/v2.0/voservice/conesearch',params={'RA':ra,'DEC':dec,'SR':0.08},timeout=120)
  r.raise_for_status();parse_single_table(__import__('io').BytesIO(r.content))
  (p/'lamost.xml').write_bytes(r.content)
images={'NGC6791':'ibc603zoq_drc.fits','NGC6819':'ib2o06020_drc.fits'}
def job(name,kind):
 p=ROOT/name
 if kind=='image':
  filename=images[name];url='https://mast.stsci.edu/api/v0.1/Download/file';params={'uri':'mast:HST/product/'+filename};dest=p/filename
 else:
  t=parse_single_table(p/'lamost.xml').to_table(use_names_over_ids=True);ids=sorted(set(str(r['gaia_source_id']) for r in t if not np.ma.is_masked(r['gaia_source_id']) and str(r['gaia_source_id']).isdigit() and int(r['gaia_source_id'])>0))
  q='select source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,ruwe from gaiadr3.gaia_source where source_id in ('+','.join(ids)+')'
  url='https://gea.esac.esa.int/tap-server/tap/sync';params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':q};dest=p/'gaia.json'
 if not dest.exists():
  with requests.get(url,params=params,stream=True,timeout=(30,180)) as r:
   r.raise_for_status()
   with dest.with_suffix('.partial').open('wb') as f:
    for chunk in r.iter_content(1024*1024):f.write(chunk)
  dest.with_suffix('.partial').replace(dest)
 print(name,kind,dest.stat().st_size,flush=True)
with ThreadPoolExecutor(max_workers=4) as e:
 fs=[e.submit(job,n,k) for n in images for k in ['image','gaia']]
 for f in as_completed(fs):f.result()
(ROOT/'selected_images.json').write_text(json.dumps(images,indent=2))
