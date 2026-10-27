"""Broad HST/LAMOST overlap discovery; centers are search regions, not memberships."""
import json,argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from search_fields import ROOT,DATA,FIELDS,search_field
MORE={'M67':(132.825,11.8),'M11':(282.77,-6.27),'M71':(298.444,18.779),'M56':(289.148,30.184),
'M29':(305.98,38.52),'NGC6811':(294.34,46.38),'NGC6866':(300.98,44.16),
'NGC663':(26.5,61.25),'NGC654':(26.0,61.88),'NGC457':(19.9,58.28),
'NGC869':(34.75,57.13),'NGC884':(35.6,57.15),'NGC1501':(61.75,60.92),
'NGC1502':(61.93,62.32),'NGC1907':(82.0,35.33),'NGC1912':(82.18,35.85),
'NGC2099':(88.08,32.55),'NGC2168':(92.28,24.33),'NGC2323':(105.68,-8.36),
'NGC2345':(107.08,-13.2),'NGC2506':(120.01,-10.77),'NGC2539':(122.67,-12.83),
'NGC2682_east':(132.95,11.8),'NGC6709':(282.83,10.32),'NGC6793':(290.62,22.14),
'NGC6940':(308.62,28.28),'NGC6939':(307.9,60.66),'NGC7142':(326.29,65.77),
'NGC7160':(327.71,62.60),'NGC7243':(333.78,49.88),'NGC7380':(341.75,58.1),
'NGC7510':(347.79,60.57),'NGC7654':(351.2,61.6),'NGC7788':(359.42,61.4),
'NGC7790':(359.6,61.2),'NGC225':(10.85,61.77),'NGC559':(22.4,63.3),
'NGC1664':(72.2,43.68),'NGC1778':(77.0,37.03),'NGC1960':(84.07,34.14),
'NGC2281':(102.26,41.08),'NGC2437':(115.44,-14.81),'NGC2395':(111.78,13.61),
'NGC2244':(98.0,4.95),'NGC2232':(96.99,-4.75),'NGC2169':(92.12,13.97),
'NGC2024':(85.42,-1.9),'NGC2068':(86.67,0.05),'IC348':(56.13,32.16),
'NGC1333':(52.3,31.3),'NGC7129':(325.76,66.1),'NGC7023':(315.4,68.15),
'NGC6535':(270.96,-0.3),'NGC6712':(283.27,-8.71),'Pal5':(229.02,-0.11),
'Pal4':(172.32,28.97),'Pal3':(151.38,0.07),'NGC6934':(308.55,7.4),
'NGC7492':(347.11,-15.61),'M31_NE':(11.1,41.6),'M31_NE2':(11.5,41.9),
'M31_SW':(10.3,40.9),'M33':(23.46,30.66),'NGC6791_outer':(290.35,37.9)}

def main():
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['cached','new','all'],default='all');args=p.parse_args()
 jobs=[]
 if args.phase in ['cached','all']:jobs += [(n,c,.25) for n,c in FIELDS.items()]
 if args.phase in ['new','all']:jobs += [(n,c,.12 if not n.startswith('M31') else .25) for n,c in MORE.items()]
 (ROOT/f'search_plan_{args.phase}.json').write_text(json.dumps([{'field':n,'ra':c[0],'dec':c[1],'radius_deg':r} for n,c,r in jobs],indent=2))
 results=[]
 with ThreadPoolExecutor(4) as pool:
  fs={pool.submit(search_field,(n,c),r):n for n,c,r in jobs}
  for f in as_completed(fs):
   try:r=f.result()
   except Exception as e:r={'field':fs[f],'error':str(e)}
   results.append(r);print(json.dumps(r),flush=True)
   (ROOT/f'search_results_{args.phase}.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
