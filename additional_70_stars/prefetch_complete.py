"""Cache evidence for complete new matches while further field searches run."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from validate_matches import ROOT,DATA
from retrieve_matched_stars import fetch_spectrum
from download_evidence import cutout
if __name__=='__main__':
 ms=json.loads((ROOT/'matches.json').read_text());g=json.loads((ROOT/'gaia_sources.json').read_text())
 ms=[m for m in ms if all(g[m['gaia_dr3_source_id']].get(k) is not None and g[m['gaia_dr3_source_id']][k]>0 for k in ['radius_flame','lum_flame'])]
 with ThreadPoolExecutor(4) as pool:
  fs={}
  for m in ms:
   fs[pool.submit(fetch_spectrum,m['lamost_obsid'],DATA/'spectra')]=(m['gaia_dr3_source_id'],'spectrum')
   fs[pool.submit(cutout,m)]=(m['gaia_dr3_source_id'],'image')
  for f in as_completed(fs):
   try:f.result()
   except Exception as e:print('FAILED',fs[f],repr(e),flush=True)
 print('Prefetch finished for',len(ms),'stars',flush=True)
