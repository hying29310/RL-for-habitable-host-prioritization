"""Search larger mosaic footprints and neighboring tiles; distinct stars are deduplicated."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from search_fields import search_field,ROOT
import json
JOBS=[('M35_wide',(92.28,24.33),.5),('M31_NE1',(11.1,41.5),.3),('M31_NE2',(11.6,41.75),.3),('M31_NE3',(12.1,42.0),.3),('M31_SW1',(10.3,41.0),.3),('M31_SW2',(9.8,40.7),.3),('M33',(23.46,30.66),.3),('COSMOS',(150.12,2.21),.5),('COSMOS_N',(150.12,2.8),.4),('COSMOS_S',(150.12,1.65),.4),('GOODS_N',(189.23,62.24),.18),('EGS',(214.83,52.82),.3),('NGC6791_wide',(290.22,37.77),.5),('Orion_wide',(83.82,-5.39),.5)]
if __name__=='__main__':
 (ROOT/'mosaic_search_plan.json').write_text(json.dumps(JOBS,indent=2))
 with ThreadPoolExecutor(3) as pool:
  fs={pool.submit(search_field,(n,c),r):n for n,c,r in JOBS}
  for f in as_completed(fs):
   try:print(json.dumps(f.result()),flush=True)
   except Exception as e:print(json.dumps({'field':fs[f],'error':str(e)}),flush=True)
