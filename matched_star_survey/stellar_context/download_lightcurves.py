"""Download one archive light-curve segment per position-checked target.
This is an availability/sample download, not a complete transit-search dataset.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,requests,hashlib
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
ROOT=Path(__file__).resolve().parent;DATA=ROOT/'data'
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def choose(star):
    gid=star['gaia_dr3_source_id'];nj=json.loads((ROOT.parent/'data/gaia_neighbors'/f'{gid}.json').read_text())
    cols=[x['name'] for x in nj['metadata']];neighbors=[dict(zip(cols,r)) for r in nj['data']];good=[]
    for r in star['light_curve_candidate_records']:
        if r['separation_from_gaia_2016_arcsec']>1:continue
        co=SkyCoord(r['s_ra']*u.deg,r['s_dec']*u.deg)
        ds=[(float(co.separation(SkyCoord(n['ra']*u.deg,n['dec']*u.deg)).arcsec),str(n['source_id'])) for n in neighbors]
        ds.sort()
        if not ds or ds[0][1]!=gid:continue
        if r['obs_collection']=='Kepler' or r['provenance_name'] in ['SPOC','TESS-SPOC','QLP']:good.append(r)
    good.sort(key=lambda r:(0 if r['obs_collection']=='Kepler' else 1,r['separation_from_gaia_2016_arcsec'],r['obs_id']))
    return good[0] if good else None

def download(star):
    gid=star['gaia_dr3_source_id'];r=choose(star)
    if r is None:return {'gaia_dr3_source_id':gid,'status':'no_selected_standard_product'}
    uri=r['dataURL'];product=None
    if r['obs_collection']=='Kepler':
        p=DATA/f"products_{r['obsid']}.json"
        if not p.exists():
            req={'service':'Mast.Caom.Products','params':{'obsid':str(r['obsid'])},'format':'json','pagesize':50000}
            rr=requests.post('https://mast.stsci.edu/api/v0/invoke',data={'request':json.dumps(req)},timeout=90);rr.raise_for_status();j=rr.json()
            if j.get('status')!='COMPLETE':raise ValueError('Incomplete MAST product response')
            write(p,j)
        prods=json.loads(p.read_text())['data'];prods=[x for x in prods if x.get('productFilename','').endswith('_llc.fits')]
        if not prods:raise ValueError('No individual long-cadence FITS product')
        prods.sort(key=lambda x:x['productFilename']);product=prods[len(prods)//2];uri=product['dataURI']
    if not uri:raise ValueError('Missing download URI')
    p=DATA/'lightcurves'/f'{gid}.fits';p.parent.mkdir(exist_ok=True)
    if not p.exists():
        url='https://mast.stsci.edu/api/v0.1/Download/file'
        rr=requests.get(url,params={'uri':uri},timeout=(30,120));rr.raise_for_status()
        tmp=p.with_suffix('.partial');tmp.write_bytes(rr.content)
        with fits.open(tmp) as h:assert len(h)>1 and h[1].data is not None
        tmp.replace(p)
    with fits.open(p) as h:
        hd=h[1];names=hd.columns.names;time=np.asarray(hd.data['TIME'],float)
        fluxcol=next((x for x in ['PDCSAP_FLUX','KSPSAP_FLUX','SAP_FLUX'] if x in names),None)
        if fluxcol is None:raise ValueError('No supported flux column: '+str(names))
        flux=np.asarray(hd.data[fluxcol],float);errcol=next((x for x in [fluxcol+'_ERR','KSPSAP_FLUX_ERR','SAP_FLUX_ERR'] if x in names),None)
        err=np.asarray(hd.data[errcol],float) if errcol else np.full_like(flux,np.nan)
        qualcol=next((x for x in ['QUALITY','SAP_QUALITY'] if x in names),None);qual=np.asarray(hd.data[qualcol],np.int64) if qualcol else np.zeros(len(time),dtype=np.int64)
        finite=np.isfinite(time)&np.isfinite(flux);usable=finite&(qual==0)
        assert usable.any()
        bjdref=float(hd.header.get('BJDREFI',0))+float(hd.header.get('BJDREFF',0))
        meta={k:h[0].header.get(k) for k in ['TELESCOP','MISSION','OBJECT','KEPLERID','TICID','RA_OBJ','DEC_OBJ','QUARTER','SECTOR']}
        meta['CROWDSAP']=hd.header.get('CROWDSAP');meta['FLFRCSAP']=hd.header.get('FLFRCSAP')
        if r['obs_collection']=='Kepler':
            expected=int(r['target_name'].replace('kplr',''));assert int(meta['KEPLERID'])==expected
        elif meta['TICID'] is not None:assert str(meta['TICID']).lstrip('0')==str(r['target_name']).lstrip('0')
    out=ROOT/'lightcurves';out.mkdir(exist_ok=True);npz=out/f'{gid}.npz'
    np.savez_compressed(npz,time_days=time,time_bjd=time+bjdref,flux=flux,flux_error=err,quality=qual,finite=finite,usable=usable)
    return {'gaia_dr3_source_id':gid,'field':star['field'],'status':'downloaded','observation':r,'product':product,'data_uri':uri,
        'fits_file':str(p.resolve()),'arrays_file':str(npz.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
        'flux_column':fluxcol,'flux_error_column':errcol,'quality_column':qualcol,'bjd_reference':bjdref,'metadata':meta,
        'samples':len(time),'usable_samples':int(usable.sum()),'time_span_days':float(np.ptp(time[finite])),
        'caution':'One segment only; catalog position/target-ID checked, but aperture blending and transit suitability not assessed. MAST target position compared at Gaia reference epoch, not fully epoch-corrected.'}

def main():
    stars=json.loads((ROOT/'availability.json').read_text());results=[]
    with ThreadPoolExecutor(4) as pool:
        fs={pool.submit(download,s):s for s in stars}
        for f in as_completed(fs):
            try:r=f.result()
            except Exception as e:r={'gaia_dr3_source_id':fs[f]['gaia_dr3_source_id'],'status':'error','error':str(e)}
            results.append(r);print(r['gaia_dr3_source_id'],r['status'],r.get('samples',''),r.get('error',''),flush=True)
    write(ROOT/'lightcurve_manifest.json',results)
if __name__=='__main__':main()
