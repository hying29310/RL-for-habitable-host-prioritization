"""Cache LAMOST spectra for new matches with known Gaia FLAME properties."""
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from search_fields import ROOT,DATA
from retrieve_matched_stars import fetch_spectrum
ms=json.loads((ROOT/'matches.json').read_text());g=json.loads((ROOT/'gaia_stellar_parameters.json').read_text())
ms=[m for m in ms if g.get(m['gaia_dr3_source_id'],{}).get('radius_flame') is not None and g.get(m['gaia_dr3_source_id'],{}).get('lum_flame') is not None]
print('Downloading',len(ms),'spectra',flush=True);errors=[]
with ThreadPoolExecutor(4) as pool:
 fs={pool.submit(fetch_spectrum,m['lamost_obsid'],DATA/'spectra'):m for m in ms}
 for i,f in enumerate(as_completed(fs),1):
  try:f.result()
  except Exception as e:errors.append({'obsid':fs[f]['lamost_obsid'],'error':str(e)})
  if i%10==0:print('Done',i,flush=True)
(ROOT/'spectrum_prefetch_errors.json').write_text(json.dumps(errors,indent=2));print('Failures',len(errors),flush=True)
