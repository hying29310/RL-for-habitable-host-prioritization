"""Real Hubble FITS -> candidate point sources -> WCS sky coordinates -> JSON.
Engineering smoke test, not validated stellar classification or planet search.
"""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR', str(ROOT / '.mplconfig'))
import json
import hashlib
from datetime import datetime, timezone
import numpy as np
import requests
import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.stats import SigmaClip
from astropy.visualization import simple_norm
from scipy.ndimage import binary_dilation
import photutils
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAME = 'icdm0a070_drc.fits'
URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/' + NAME

def save(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def main():
    path = ROOT / 'data' / NAME
    path.parent.mkdir(exist_ok=True)
    if not path.exists():
        # Stream into a temporary file so interrupted downloads are never cached.
        temporary = path.with_suffix('.partial')
        with requests.get(URL, stream=True, timeout=(30,180)) as response:
            response.raise_for_status()
            with temporary.open('wb') as out:
                for chunk in response.iter_content(1024*1024):
                    out.write(chunk)
        with fits.open(temporary) as hd:
            assert hd[0].header['TELESCOP'] == 'HST'
            hd.verify('exception')
        temporary.replace(path)
    with fits.open(path) as hd:
        hd.verify('exception')
        primary = hd[0].header.copy()
        sci_header = hd['SCI'].header.copy()
        assert primary['TELESCOP'] == 'HST'
        science = hd['SCI'].data
        full_shape = list(science.shape)
        wcs = WCS(sci_header, fobj=hd).celestial
        # Fixed central crop selected independently of any detected sources.
        cutout = Cutout2D(science, ((science.shape[1]-1)/2,
                         (science.shape[0]-1)/2), (1024,1024), wcs=wcs)
        data = np.array(cutout.data, dtype=float)
        weight = np.array(hd['WHT'].data[cutout.slices_original], dtype=float)
        overview = np.array(science[::4,::4], dtype=float)
    mask = (~np.isfinite(data)) | (~np.isfinite(weight)) | (weight <= 0)
    background = Background2D(data, (64,64), filter_size=(3,3), mask=mask,
                              sigma_clip=SigmaClip(sigma=3),
                              bkg_estimator=MedianBackground())
    residual = data - background.background
    rms = background.background_rms
    bad = mask | (~np.isfinite(rms)) | (rms <= 0)
    # Exclude detections near missing coverage and the crop edge.
    detection_mask = binary_dilation(bad, iterations=8)
    detection_mask[:8,:] = detection_mask[-8:,:] = True
    detection_mask[:,:8] = detection_mask[:,-8:] = True
    detection_image = np.divide(residual, rms, out=np.zeros_like(data), where=~bad)
    finder = DAOStarFinder(threshold=5.0, fwhm=2.5, exclude_border=True)
    sources = finder(detection_image, mask=detection_mask)
    assert sources is not None and len(sources) > 0, 'No candidate point sources detected'
    positions = np.column_stack((sources['xcentroid'], sources['ycentroid']))
    aperture = CircularAperture(positions, r=3)
    sums = aperture_photometry(residual, aperture, mask=bad)['aperture_sum']
    # Flux is a small-aperture estimate, not aperture-corrected total flux.
    x,y = positions.T
    ra,dec = cutout.wcs.all_pix2world(x,y,0)
    backx,backy = cutout.wcs.all_world2pix(ra,dec,0)
    roundtrip = float(np.max(np.hypot(backx-x,backy-y)))
    x0,y0 = cutout.origin_original
    rows=[]
    for i in range(len(sources)):
        rows.append({'source_id':f'HST_DET_{i+1:04d}',
                     'x_cutout_zero_based':float(x[i]), 'y_cutout_zero_based':float(y[i]),
                     'x_full_image_zero_based':float(x[i]+x0),
                     'y_full_image_zero_based':float(y[i]+y0),
                     'ra_deg':float(ra[i]), 'dec_deg':float(dec[i]),
                     'aperture_radius_pixels':3,
                     'background_subtracted_aperture_flux':float(sums[i]),
                     'flux_unit':sci_header.get('BUNIT'),
                     'sharpness':float(sources['sharpness'][i]),
                     'roundness1':float(sources['roundness1'][i]),
                     'temperature_K':None,'metallicity_Fe_H_dex':None,
                     'catalog_match_status':'not_attempted',
                     'classification':'candidate_point_source_not_confirmed_star'})
    checks={'verified_HST_header':primary['TELESCOP']=='HST',
            'celestial_WCS':bool(cutout.wcs.has_celestial),
            'finite_valid_pixels':bool(np.isfinite(data[~mask]).all()),
            'sources_detected':len(rows)>0,
            'finite_source_coordinates':bool(np.isfinite(ra).all() and np.isfinite(dec).all()),
            'WCS_roundtrip_under_1e_minus_4_pixel':roundtrip<1e-4,
            'finite_aperture_fluxes':bool(np.isfinite(sums).all())}
    header = cutout.wcs.to_header()
    header['TELESCOP']='HST'
    header['INSTRUME']=primary['INSTRUME']
    header['FILTER']=primary['FILTER']
    header['BUNIT']=sci_header.get('BUNIT','UNKNOWN')
    header['ORIGFILE']=NAME
    fits.HDUList([fits.PrimaryHDU(data.astype('float32'), header=header),
                  fits.ImageHDU(weight.astype('float32'),name='WHT')]).writeto(ROOT/'science_cutout.fits',overwrite=True)
    # Verify that serialized WCS preserves measured sky coordinates.
    with fits.open(ROOT/'science_cutout.fits') as hd:
        rr,dd=WCS(hd[0].header).all_pix2world(x,y,0)
    checks['saved_cutout_WCS_preserved']=bool(np.allclose(rr,ra,atol=1e-9,rtol=0) and np.allclose(dd,dec,atol=1e-9,rtol=0))
    output={'status':'PASS' if all(checks.values()) else 'FAIL',
            'run_utc':datetime.now(timezone.utc).isoformat(),
            'observation':{k:primary.get(k) for k in ['TELESCOP','INSTRUME','DETECTOR','TARGNAME','FILTER','DATE-OBS','EXPTIME']},
            'source_url':URL,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'source_bytes':path.stat().st_size, 'original_shape_yx':full_shape,
            'cutout_shape_yx':list(data.shape),'cutout_origin_xy':[int(x0),int(y0)],
            'candidate_source_count':len(rows),'checks':checks,
            'max_WCS_roundtrip_error_pixels':roundtrip,
            'settings':{'background_box_pixels':64,'DAO_fwhm_pixels':2.5,
                'DAO_threshold_on_local_RMS_normalized_image':5.0,'aperture_radius_pixels':3},
            'limitations':['Candidate point sources may include clusters, compact galaxies, blends, and artifacts.',
                'Threshold is a detection setting, not calibrated statistical significance.',
                'Fluxes lack aperture corrections, uncertainty estimates, and crowding corrections.',
                'WCS roundtrip checks numerical consistency, not absolute astrometric accuracy.',
                'No Gaia/LAMOST match attempted; temperature and metallicity deliberately null.',
                'No transit search, habitability inference, or RL training performed.'],
            'versions':{'numpy':np.__version__,'astropy':astropy.__version__,'photutils':photutils.__version__}}
    save('results.json',output)
    save('sources.json',{'observation':output['observation'],'pixel_origin':0,'sources':rows})
    assert len(json.loads((ROOT/'sources.json').read_text())['sources'])==len(rows)
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    axes[0].imshow(overview,origin='lower',cmap='gray',norm=simple_norm(overview,'asinh',percent=99.5))
    axes[0].add_patch(plt.Rectangle((x0/4,y0/4),1024/4,1024/4,fill=False,edgecolor='#00d6df',lw=1.5))
    axes[0].set(title='Hubble WFC3/UVIS · F814W\nFull field; cyan box = analyzed crop',xlabel='Overview column (4× downsampled)',ylabel='Overview row (4× downsampled)')
    axes[1].imshow(data,origin='lower',cmap='gray',norm=simple_norm(data[~mask],'asinh',percent=99.5))
    axes[1].scatter(x,y,s=26,facecolors='none',edgecolors='#00e5ef',lw=.7)
    axes[1].set(title=f'{len(rows)} candidate point sources\n1024 × 1024 pixel crop',xlabel='Cutout column (zero-based)',ylabel='Cutout row (zero-based)')
    fig.suptitle(f"Real Hubble image: {primary['TARGNAME']} · {primary['DATE-OBS']}\nImage → source detection → sky coordinates → JSON",fontsize=13)
    fig.savefig(ROOT/'diagnostic.png',dpi=160)
    plt.close(fig)
    print(json.dumps(output,indent=2))
    assert all(checks.values()), 'Inspect failed checks in results.json'

if __name__=='__main__':
    main()
