"""Download a measured LAMOST spectrum and a real HST FITS cutout per accepted star."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import json,requests,hashlib,io
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from validate_matches import ROOT,DATA,save
from retrieve_matched_stars import fetch_spectrum,retrieve

def cutout(m):
    """Get a 101-pixel HLA cutout using the catalog's astrometric correction."""
    dest=DATA/'hst_cutouts'/f"{m['gaia_dr3_source_id']}.fits";dest.parent.mkdir(exist_ok=True)
    d=m['best_hst_detection'];params={'red':d['ImageName'],'ra':d['SourceRA'],'dec':d['SourceDec'],
        'size':101,'format':'fits','applyomega':'true'}
    url=requests.Request('GET','https://hla.stsci.edu/cgi-bin/fitscut.cgi',params=params).prepare().url
    if not dest.exists():
        r=requests.get(url,timeout=(30,120));r.raise_for_status()
        with fits.open(io.BytesIO(r.content)) as h:
            assert any(x.data is not None and x.data.ndim==2 for x in h)
        tmp=dest.with_suffix('.partial');tmp.write_bytes(r.content);tmp.replace(dest)
    with fits.open(dest) as h:
        hd=next(x for x in h if x.data is not None and x.data.ndim==2)
        a=np.array(hd.data,float);w=WCS(hd.header).celestial
        # The HLA ApplyOmega adjustment may center using corrected coordinates while
        # retaining original WCS headers; a cutout is visual evidence, not a new fit.
        finite=int(np.isfinite(a).sum());assert finite>0
    return {'gaia_dr3_source_id':m['gaia_dr3_source_id'],'hst_image':d['ImageName'],'path':str(dest.resolve()),
            'url':url,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'shape':list(a.shape),'finite_pixels':finite,
            'note':'Published HSC coordinates provide the measured centroid. Cutout downloaded for visual review; no new independent centroid claimed.'}

def main():
    matches=json.loads((ROOT/'matches.json').read_text());fail=[];cuts=[]
    with ThreadPoolExecutor(4) as pool:
        tasks={}
        for m in matches:
            tasks[pool.submit(fetch_spectrum,m['lamost_obsid'],DATA/'spectra')]=('spectrum',m)
            tasks[pool.submit(cutout,m)]=('cutout',m)
        for f in as_completed(tasks):
            kind,m=tasks[f]
            try:
                result=f.result()
                if kind=='cutout':cuts.append(result)
            except Exception as e:fail.append({'kind':kind,'gaia_dr3_source_id':m['gaia_dr3_source_id'],'error':str(e)})
    save(ROOT/'cutout_manifest.json',cuts);save(ROOT/'download_failures.json',fail)
    result=retrieve(ROOT/'matches.json',ROOT/'selected_lamost.xml',DATA/'spectra',ROOT/'retrieved')
    print(json.dumps({'cutouts':len(cuts),'spectra_retrieval':result,'download_failures':fail}),flush=True)
if __name__=='__main__':main()
