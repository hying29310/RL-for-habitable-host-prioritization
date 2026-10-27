"""Search additional northern Hubble target regions without changing match thresholds."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from search_fields import search_field,ROOT
import json
FIELDS={'M5':(229.638,2.082),'M10':(254.288,-4.100),'M12':(251.809,-1.948),'M14':(264.400,-3.245),'M71':(298.443,18.779),'M56':(289.148,30.184),'NGC6712':(283.268,-8.706),'NGC6535':(270.960,-0.298),'M37':(88.074,32.553),'M38':(82.183,35.824),'M11':(282.77,-6.27),'M52':(351.2,61.6),'NGC663':(26.55,61.25),'NGC457':(19.9,58.28),'NGC869':(34.75,57.15),'NGC884':(35.60,57.13),'NGC188':(11.80,85.25),'M67_extended':(132.825,11.80),'NGC7142':(326.3,65.77),'NGC6793':(290.95,22.13),'NGC6811':(294.34,46.38),'NGC6866':(300.98,44.16),'NGC7160':(328.42,62.60),'NGC2168':(92.28,24.33),'NGC2169':(92.13,13.96),'NGC2244':(98.0,4.95),'NGC1980':(83.85,-5.90),'NGC1977':(83.85,-4.84),'IC348':(56.14,32.17),'NGC1333':(52.3,31.3),'NGC1501':(61.75,60.92),'NGC2301':(102.93,0.46)}
if __name__=='__main__':
 (ROOT/'extra_search_centers.json').write_text(json.dumps(FIELDS,indent=2))
 with ThreadPoolExecutor(3) as pool:
  fs={pool.submit(search_field,(n,c),.12):n for n,c in FIELDS.items()}
  for f in as_completed(fs):
   try:print(json.dumps(f.result()),flush=True)
   except Exception as e:print(json.dumps({'field':fs[f],'error':str(e)}),flush=True)
