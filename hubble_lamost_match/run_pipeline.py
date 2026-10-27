"""Real HST image detection joined to a quality-filtered LAMOST DR11 spectrum.
Run with Python 3.12 and requirements-lock.txt. Missing inputs are downloaded.
"""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplconfig'))
import json, hashlib
from datetime import datetime,timezone
import requests
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.stats import SigmaClip
from astropy.visualization import simple_norm
import astropy.units as u
from photutils.background import Background2D, MedianBackground
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture, aperture_photometry
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA=ROOT/'data'
HST_URL='https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/icmp09011_drz.fits'
CAT_URL='https://www.lamost.org/dr11/v2.0/voservice/conesearch?RA=132.825&DEC=11.8&SR=0.05'
GAIA_QUERY='select source_id,ra,dec,ref_epoch,pmra,pmdec,parallax,ruwe from gaiadr3.gaia_source where source_id in (604917526275831040,604917285757663872,604917491916095872)'

def fetch(name,url,params=None):
    p=DATA/name
    if not p.exists():
        temp=p.with_suffix('.partial')
        with requests.get(url,params=params,stream=True,timeout=(30,180)) as r:
            r.raise_for_status()
            with temp.open('wb') as f:
                for chunk in r.iter_content(1024*1024): f.write(chunk)
        temp.replace(p)
    return p

def save(name,obj):
    (ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def main():
    DATA.mkdir(exist_ok=True)
    hst_path=fetch('icmp09011_drz.fits',HST_URL)
    cat_path=fetch('lamost_cone.xml',CAT_URL)
    gaia_path=fetch('gaia.json','https://gea.esac.esa.int/tap-server/tap/sync',
        {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':GAIA_QUERY})
    cat=parse_single_table(cat_path).to_table(use_names_over_ids=True)
    # Conservative demonstration filter; not a universal survey selection function.
    clean=cat[(cat['fibermask']==0)&(cat['snrg']>20)&(cat['teff']>0)&(cat['feh']>-5)]
    assert len(clean)>0
    clean.sort('snrg',reverse=True)
    row=clean[0]
    obsid=int(row['obsid'])
    spectrum_url=f'https://www.lamost.org/dr11/v2.0/spectrum/fits/{obsid}'
    spec_path=fetch(f'lamost_{obsid}.fits.gz',spectrum_url)
    with fits.open(hst_path) as h:
        h.verify('exception')
        hdr=h[0].header.copy()
        sh=h['SCI'].header.copy()
        image=np.array(h['SCI'].data,dtype=float)
        weight=np.array(h['WHT'].data,dtype=float)
        wcs=WCS(sh,fobj=h).celestial
    mask=(~np.isfinite(image))|(~np.isfinite(weight))|(weight<=0)
    bkg=Background2D(image,(64,64),filter_size=(3,3),mask=mask,
        sigma_clip=SigmaClip(sigma=3),bkg_estimator=MedianBackground())
    resid=image-bkg.background
    rms=bkg.background_rms
    z=np.divide(resid,rms,out=np.zeros_like(resid),where=(rms>0)&~mask)
    detections=DAOStarFinder(threshold=5,fwhm=2.5,exclude_border=True)(z,mask=mask)
    assert detections is not None
    x=np.asarray(detections['xcentroid']);y=np.asarray(detections['ycentroid'])
    ra,dec=wcs.all_pix2world(x,y,0)
    sky=SkyCoord(ra*u.deg,dec*u.deg)
    gj=json.loads(gaia_path.read_text())
    gaia=[dict(zip([c['name'] for c in gj['metadata']],r)) for r in gj['data']]
    g=next(r for r in gaia if int(r['source_id'])==int(row['gaia_source_id']))
    # Gaia pmra is mu_alpha*cos(delta). Positive parallax supplies a distance;
    # radial velocity is set to zero for this short (~1 year) propagation.
    gc=SkyCoord(ra=g['ra']*u.deg,dec=g['dec']*u.deg,
        pm_ra_cosdec=g['pmra']*u.mas/u.yr,pm_dec=g['pmdec']*u.mas/u.yr,
        distance=Distance(parallax=g['parallax']*u.mas),radial_velocity=0*u.km/u.s,
        obstime=Time(g['ref_epoch'],format='jyear'))
    at_hst=gc.apply_space_motion(new_obstime=Time(hdr['EXPSTART'],format='mjd'))
    sep=sky.separation(at_hst).arcsec
    order=np.argsort(sep);k=int(order[0])
    rawsep=sky[k].separation(SkyCoord(float(row['ra'])*u.deg,float(row['dec'])*u.deg)).arcsec
    aperture=CircularAperture((x[k],y[k]),r=3)
    flux=float(aperture_photometry(resid,aperture,mask=mask)['aperture_sum'][0])
    neighbors=sky.separation(sky[k]).arcsec
    neighbor_count=int(np.sum((neighbors>0.001)&(neighbors<1.65)))
    with fits.open(spec_path) as h:
        # Released file contains stray NAXIS1/NAXIS2 in an empty primary HDU.
        # Fix that metadata in memory only; preserve original bytes on disk.
        h.verify('silentfix+exception')
        sph=h[0].header.copy();s=h['COADD'].data[0]
        wave=np.array(s['WAVELENGTH'],float);sf=np.array(s['FLUX'],float)
        iv=np.array(s['IVAR'],float);orm=np.array(s['ORMASK'],int)
        good=np.isfinite(wave)&np.isfinite(sf)&(iv>0)&(orm==0)
    # Plot relative flux; do not assert absolute flux calibration.
    sfrel=sf/np.median(sf[good])
    source='HST_DET_'+str(k+1).zfill(4)
    match={'hst_source_id':source,'lamost_designation':str(row['designation']),
        'lamost_obsid':obsid,'lamost_release':'DR11 v2.0 LRS',
        'gaia_dr3_source_id':str(g['source_id']),
        'hst_centroid_xy_zero_based':[float(x[k]),float(y[k])],
        'hst_ra_deg':float(ra[k]),'hst_dec_deg':float(dec[k]),
        'gaia_at_hst_ra_deg':float(at_hst.ra.deg),'gaia_at_hst_dec_deg':float(at_hst.dec.deg),
        'hst_gaia_epoch_corrected_separation_arcsec':float(sep[k]),
        'hst_lamost_catalog_separation_arcsec':float(rawsep),
        'second_nearest_detection_to_gaia_arcsec':float(sep[order[1]]),
        'other_detections_within_1p65_arcsec':neighbor_count,
        'gaia_ruwe':g['ruwe'],'temperature_K':float(row['teff']),
        'temperature_catalog_error_K':float(row['teff_err']),
        'Fe_H_dex':float(row['feh']),'Fe_H_catalog_error_dex':float(row['feh_err']),
        'logg_dex':float(row['logg']),'logg_catalog_error_dex':float(row['logg_err']),
        'spectral_subclass':str(row['subclass']),'snrg':float(row['snrg']),
        'fibermask':int(row['fibermask']),'lamost_observation_date':str(row['obsdate']),
        'aperture_flux_electrons_per_second':flux,'spectrum_url':spectrum_url,
        'stellar_parameter_origin':'LAMOST catalog; not estimated from Hubble pixels',
        'uncertainty_note':'Quoted catalog fit errors do not include all systematic uncertainty.'}
    bx,by=wcs.all_world2pix(ra,dec,0)
    spectrum_position=SkyCoord(sph['RA']*u.deg,sph['DEC']*u.deg)
    checks={'image_is_HST':hdr['TELESCOP']=='HST','in_LAMOST_declination_region':float(dec[k])>-10,
        'detected_source_match_under_0p3_arcsec':bool(sep[k]<0.3),
        'one_detection_within_0p3_arcsec':int(np.sum(sep<0.3))==1,
        'no_additional_detections_within_1p65_arcsec':neighbor_count==0,
        'zero_LAMOST_fibermask':int(row['fibermask'])==0 and sph['FIB_MASK']==0,
        'catalog_spectrum_obsid_agree':int(sph['OBSID'])==obsid,
        'catalog_spectrum_position_agree':spectrum_position.separation(SkyCoord(float(row['ra'])*u.deg,float(row['dec'])*u.deg)).arcsec<0.01,
        'spectrum_has_over_3000_good_bins':int(good.sum())>3000,
        'finite_stellar_parameters':bool(np.isfinite([row['teff'],row['feh']]).all()),
        'WCS_roundtrip_under_1e_minus_4_pixel':bool(np.max(np.hypot(bx-x,by-y))<1e-4)}
    checks={key:bool(value) for key,value in checks.items()}
    sources=[]
    for i in range(len(x)):
        sources.append({'source_id':'HST_DET_'+str(i+1).zfill(4),'x':float(x[i]),'y':float(y[i]),
            'ra_deg':float(ra[i]),'dec_deg':float(dec[i]),
            'temperature_K':match['temperature_K'] if i==k else None,
            'Fe_H_dex':match['Fe_H_dex'] if i==k else None,
            'lamost_obsid':obsid if i==k else None,
            'match_status':'validated_example' if i==k else 'not_validated'})
    save('sources.json',sources)
    save('matched_star.json',match)
    np.savez_compressed(ROOT/'spectrum_arrays.npz',wavelength_angstrom=wave,
        catalog_flux=sf,inverse_variance=iv,ormask=orm,good=good)
    result={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,
        'run_utc':datetime.now(timezone.utc).isoformat(),'hst_file':hst_path.name,
        'hst_target':hdr['TARGNAME'],'hst_filter':hdr['FILTER'],'hst_date':hdr['DATE-OBS'],
        'hst_exposure_seconds':hdr['EXPTIME'],'image_shape_yx':list(image.shape),
        'detection_count':len(x),'lamost_catalog_rows':len(cat),'quality_filtered_catalog_rows':len(clean),
        'catalog_query_url':CAT_URL,'hst_url':HST_URL,'gaia_query':GAIA_QUERY,
        'match':match,'spectrum_range_angstrom':[float(wave.min()),float(wave.max())],
        'usable_spectrum_bins':int(good.sum()),
        'settings':{'background_box_pixels':64,'DAO_fwhm_pixels':2.5,'threshold_local_rms':5},
        'limitations':['One validated association, not a complete field catalog.',
            'LAMOST empty-primary-HDU NAXIS1/NAXIS2 metadata fixed in memory; downloaded bytes preserved.',
            'Image detections can include artifacts; absence of detected neighbors does not prove no unresolved companion.',
            'No HST saturation/DQ audit or precision photometric calibration performed.',
            'Only zero-fibermask, SNRg>20 spectra with finite plausible parameters are used for this demonstration.',
            'No stellar parameter fitting, habitability inference, transit search, or RL training.'],
        'files':[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [hst_path,cat_path,gaia_path,spec_path]]}
    save('results.json',result)
    fig=plt.figure(figsize=(13,9),layout='constrained');gs=fig.add_gridspec(2,2)
    a=fig.add_subplot(gs[0,0]);a.imshow(image,origin='lower',cmap='gray',norm=simple_norm(image[~mask],'asinh',percent=99.5));a.scatter([x[k]],[y[k]],s=160,facecolors='none',edgecolors='cyan');a.set(title='Hubble M67 · F621M · 2014-11-09',xlabel='Pixel column',ylabel='Pixel row')
    a=fig.add_subplot(gs[0,1]);a.imshow(image,origin='lower',cmap='gray',norm=simple_norm(image[~mask],'asinh',percent=99.9));a.set_xlim(x[k]-45,x[k]+45);a.set_ylim(y[k]-45,y[k]+45)
    a.scatter([x[k]],[y[k]],s=150,facecolors='none',edgecolors='cyan',label='Detected Hubble source')
    gx,gy=wcs.all_world2pix(at_hst.ra.deg,at_hst.dec.deg,0);a.scatter([gx],[gy],marker='+',s=100,color='orange',label='Gaia at Hubble epoch');a.legend(fontsize=8)
    a.set(title=f"LAMOST {obsid} match: {sep[k]:.3f} arcsec",xlabel='Pixel column',ylabel='Pixel row')
    a=fig.add_subplot(gs[1,:]);a.plot(wave,np.where(good,sfrel,np.nan),lw=.65,color='#176b87');a.set(xlabel='Wavelength (Angstrom)',ylabel='Flux / median valid flux',title=f"Downloaded LAMOST spectrum · {row['designation']}\nTeff = {row['teff']:.0f} K; [Fe/H] = {row['feh']:.3f} dex; S/N(g) = {row['snrg']:.1f}")
    fig.suptitle('Verified Hubble image → detected star → LAMOST spectrum and stellar parameters',fontsize=14)
    fig.savefig(ROOT/'diagnostic.png',dpi=160);plt.close(fig)
    print(json.dumps({'status':result['status'],'checks':checks,'match':match,'detection_count':len(x)},indent=2))
    assert all(checks.values())

if __name__=='__main__':main()
