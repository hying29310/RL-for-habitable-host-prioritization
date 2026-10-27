"""Fetch Gaia DR3 radius/luminosity estimates for new matched stars only."""
from pathlib import Path
import json,requests,hashlib
ROOT=Path(__file__).resolve().parent
COLS=['source_id','radius_flame','radius_flame_lower','radius_flame_upper','lum_flame','lum_flame_lower','lum_flame_upper','flags_flame','radius_gspphot','radius_gspphot_lower','radius_gspphot_upper','teff_gspphot','teff_gspphot_lower','teff_gspphot_upper','libname_gspphot']
def retrieve(ids):
    out={};p=ROOT/'data/gaia_parameters';p.mkdir(exist_ok=True)
    for i in range(0,len(ids),200):
        chunk=sorted(ids)[i:i+200];key=hashlib.sha256(','.join(chunk).encode()).hexdigest()[:16];dest=p/f'{key}.json'
        if not dest.exists():
            q='SELECT '+','.join(COLS)+' FROM gaiadr3.astrophysical_parameters WHERE source_id IN ('+','.join(chunk)+')'
            r=requests.get('https://gea.esac.esa.int/tap-server/tap/sync',params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':q},timeout=120);r.raise_for_status();dest.write_text(r.text)
        j=json.loads(dest.read_text());cols=[x['name'] for x in j['metadata']]
        out.update({str(row[0]):dict(zip(cols,row)) for row in j['data']})
    (ROOT/'gaia_stellar_parameters.json').write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--preliminary',action='store_true');a=p.parse_args()
    ms=json.loads((ROOT/('preliminary_matches.json' if a.preliminary else 'matches.json')).read_text());g=retrieve([m['gaia_dr3_source_id'] for m in ms])
    print('matched',len(ms),'Gaia rows',len(g),'FLAME complete',sum(v['radius_flame'] is not None and v['lum_flame'] is not None for v in g.values()),'GSP-Phot radii',sum(v['radius_gspphot'] is not None for v in g.values()))
