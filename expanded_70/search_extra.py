from concurrent.futures import ThreadPoolExecutor,as_completed
import json
from search_fields import search_field,ROOT
JOBS=[('M35_west',(92.0,24.33),.25),('M35_east',(92.55,24.33),.25),('M31_PHAT1',(11.3,41.7),.25),('M31_PHAT2',(11.7,42.0),.25),('M31_PHAT3',(11.8,41.8),.25),('M31_SW2',(10.1,40.8),.25),('M33_north',(23.45,30.9),.25)]
with ThreadPoolExecutor(3) as pool:
 fs={pool.submit(search_field,(n,c),r):n for n,c,r in JOBS};out=[]
 for f in as_completed(fs):
  try:r=f.result()
  except Exception as e:r={'field':fs[f],'error':str(e)}
  out.append(r);print(json.dumps(r),flush=True)
  (ROOT/'search_results_extra.json').write_text(json.dumps(out,indent=2))
