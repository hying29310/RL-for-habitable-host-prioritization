"""Retrieve Gaia radius/luminosity estimates and search MAST time-series records.
Archive position matches are candidates, not proof that light curves are unblended.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,requests
from astropy.coordinates import SkyCoord
import astropy.units as u
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data';DATA.mkdir(exist_ok=True)
MATCHES=json.loads((ROOT.parent/'matches.json').read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def gaia():
    dest=DATA/'gaia_astrophysical_parameters.json'
    if not dest.exists():
        cols=['source_id','radius_flame','radius_flame_lower','radius_flame_upper','lum_flame','lum_flame_lower','lum_flame_upper','flags_flame','radius_gspphot','radius_gspphot_lower','radius_gspphot_upper','teff_gspphot','teff_gspphot_lower','teff_gspphot_upper','distance_gspphot','ag_gspphot','libname_gspphot','logposterior_gspphot']
        q='SELECT '+','.join(cols)+' FROM gaiadr3.astrophysical_parameters WHERE source_id IN ('+','.join(m['gaia_dr3_source_id'] for m in MATCHES)+')'
        r=requests.get('https://gea.esac.esa.int/tap-server/tap/sync',params={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':q},timeout=120);r.raise_for_status();j=r.json();write(dest,j)
    j=json.loads(dest.read_text());cols=[x['name'] for x in j['metadata']]
    return {str(row[0]):dict(zip(cols,row)) for row in j['data']}
def mast(m):
    dest=DATA/f"mast_{m['gaia_dr3_source_id']}.json"
    if not dest.exists():
        req={'service':'Mast.Caom.Filtered.Position','params':{'columns':'*','filters':[{'paramName':'dataproduct_type','values':['timeseries']}],'position':f"{m['ra_deg']}, {m['dec_deg']}, {5/3600}"},'format':'json','pagesize':50000}
        r=requests.post('https://mast.stsci.edu/api/v0/invoke',data={'request':json.dumps(req)},timeout=120);r.raise_for_status();j=r.json()
        if j.get('status')!='COMPLETE':raise ValueError(str(j)[:500])
        if j.get('paging',{}).get('pagesFiltered',1)>1:raise ValueError('MAST result needs pagination')
        write(dest,j)
    j=json.loads(dest.read_text());out=[];pos=SkyCoord(m['ra_deg']*u.deg,m['dec_deg']*u.deg)
    for r in j.get('data',[]):
        if r.get('s_ra') is None or r.get('s_dec') is None:continue
        sep=pos.separation(SkyCoord(r['s_ra']*u.deg,r['s_dec']*u.deg)).arcsec
        if sep<=3:
            d={k:r.get(k) for k in ['obsid','obs_id','obs_collection','provenance_name','target_name','s_ra','s_dec','t_min','t_max','t_exptime','dataURL','sequence_number','instrument_name']}
            d['separation_from_gaia_2016_arcsec']=float(sep);out.append(d)
    return m['gaia_dr3_source_id'],out

def main():
    gs=gaia();write(ROOT/'gaia_properties.json',gs);results={};errors=[]
    with ThreadPoolExecutor(4) as pool:
        fs={pool.submit(mast,m):m for m in MATCHES}
        for f in as_completed(fs):
            try:
                gid,rows=f.result();results[gid]=rows;print(gid,len(rows),sorted({r['obs_collection'] for r in rows}),flush=True)
            except Exception as e:errors.append({'gaia_dr3_source_id':fs[f]['gaia_dr3_source_id'],'error':str(e)})
    stars=[]
    for m in MATCHES:
        gid=m['gaia_dr3_source_id'];g=gs.get(gid,{});lc=results.get(gid,[])
        stars.append({'gaia_dr3_source_id':gid,'field':m['field'],'gaia_parameters':g,'light_curve_candidate_records':lc,'light_curve_query_completed':gid in results})
    write(ROOT/'availability.json',stars);write(ROOT/'query_errors.json',errors)
    s={'stars':len(stars),'Gaia_radius_FLAME':sum(x['gaia_parameters'].get('radius_flame') is not None for x in stars),'Gaia_luminosity_FLAME':sum(x['gaia_parameters'].get('lum_flame') is not None for x in stars),'Gaia_radius_GSPPhot':sum(x['gaia_parameters'].get('radius_gspphot') is not None for x in stars),'with_nearby_archive_timeseries':sum(bool(x['light_curve_candidate_records']) for x in stars),'errors':errors}
    write(ROOT/'summary.json',s);print(json.dumps(s),flush=True)
if __name__=='__main__':main()
