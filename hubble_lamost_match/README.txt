Hubble → LAMOST: successful real-data match

This test downloaded a real Hubble image of M67, detected sources, matched one source using Gaia astrometry, downloaded its actual LAMOST DR11 v2.0 spectrum, and attached LAMOST temperature and metallicity measurements. All 11 checks passed. No stellar parameters were inferred from the Hubble image.

Data and result

Item: Value
Hubble product: icmp09011_drz.fits
Instrument / filter: WFC3/UVIS / F621M
Hubble observation: 2014-11-09; 10 seconds total exposure
Image dimensions: 2,063 × 2,192 pixels
LAMOST star: J085119.17+114754.6
LAMOST observation ID: 420808077
Gaia DR3 source ID: 604917526275831040
LAMOST observation: 2016-02-06
Hubble–Gaia separation at Hubble epoch: 0.0399 arcseconds
Hubble–LAMOST catalog separation: 0.0683 arcseconds
Effective temperature: 6105.19 ± 16.18 K
Metallicity [Fe/H]: −0.145 ± 0.010 dex
Surface gravity log g: 4.241 ± 0.020 dex
LAMOST spectral subclass: F6
g-band spectral S/N: 217.48
LAMOST fiber mask: 0
Downloaded spectral range: Approximately 3700–9099 Å

Uncertainties above are the catalog's reported errors, not a complete uncertainty budget. Survey calibration, stellar modeling, and other systematic errors can be larger.

How the field and star were selected

The LAMOST public cone search at RA=132.825°, Dec=11.8°, radius=0.05° returned 22 observation rows representing three Gaia-linked stars. MAST returned overlapping M67 Hubble images. The selected short-exposure F621M image includes all three stars. Short exposure reduces, but does not establish absence of, saturation; a full HST detector-quality audit was not performed.

The example applies a conservative filter: fiber mask zero, S/N(g)>20, and available plausible temperature and metallicity. One observation passed. The other two stars' available candidate records have nonzero fiber flags; two inspected examples have fibermask=64, the documented WHOPPER flag. They are not promoted as clean matches here. The original cone response is included for audit.

Processing and validation

1. Read Hubble SCI, WHT, and celestial WCS; mask invalid or uncovered pixels.
2. Estimate local background and RMS with 64-pixel boxes.
3. Run DAOStarFinder (FWHM=2.5 pixels, threshold=5 on RMS-normalized residuals), yielding 526 candidate detections. Many can be artifacts; this is not a scientifically vetted star catalog.
4. Retrieve Gaia DR3 astrometry using the Gaia identifier already present in LAMOST. Propagate its position from epoch 2016.0 to the Hubble exposure epoch using proper motion. Radial velocity is set to zero for this short propagation; the positive Gaia parallax is used only for space-motion propagation, not a scientific distance estimate.
5. Find the nearest independently detected Hubble source. Exactly one detection is within 0.3 arcsec; the second nearest is 6.88 arcsec away. No additional pipeline detections lie within 1.65 arcsec of the matched source. Undetected or unresolved companions remain possible.
6. Download the LAMOST spectrum and verify the observation ID and coordinates against the catalog. The file has 3,909 usable spectral bins under the applied finite-flux, positive-inverse-variance, and zero-ORMASK criteria.
7. Save the association, parameters, detected sources, spectral arrays, checks, and plot. All other detections retain null stellar parameters.

The LAMOST FITS primary header contains stray NAXIS1/NAXIS2 cards for an empty primary HDU. Astropy fixes these in memory before verification. The original downloaded bytes are preserved and hashed. Spectral data are read from the COADD table. The plotted spectrum is divided by its median for display and is not presented as absolute calibrated flux.

Run

On this workspace:

cd /Users/mac/Documents/Codex/2026-09-28/new-chat/outputs/hubble_lamost_match
/Users/mac/Documents/Codex/2026-09-28/new-chat/work/venv/bin/python run_pipeline.py
open diagnostic.png

Elsewhere, with Python 3.12 installed:

python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python run_pipeline.py

The ZIP includes the four required original inputs, so cached reruns are offline. Missing inputs download automatically. This is a reproducible demonstration for the selected field, not a general all-sky search program: the M67 image and Gaia query are fixed to the field chosen during archive exploration.

Outputs

- matched_star.json: the verified example and catalog parameters.
- sources.json: 526 image detections, one joined to clean LAMOST data.
- results.json: checks, URLs, settings, timestamps, and SHA-256 provenance.
- spectrum_arrays.npz: wavelengths, original flux, inverse variance, pixel mask, and usable-bin mask.
- diagnostic.png: Hubble field, close-up of the positional match, and LAMOST spectrum.
- run_pipeline.py: standalone processing code.
- data/: Hubble image, LAMOST cone result, Gaia astrometry, and matched LAMOST spectrum.

This establishes a working image-to-spectroscopic-catalog connection. It does not establish a planet, habitability, membership of M67, a full precision-photometry solution, or performance of an RL policy.

Direct sources

- Hubble FITS: https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/icmp09011_drz.fits
- LAMOST spectrum record: https://www.lamost.org/dr11/v2.0/spectrum/view?obsid=420808077
- LAMOST spectrum FITS: https://www.lamost.org/dr11/v2.0/spectrum/fits/420808077
- LAMOST DR11 data definitions and masks: https://www.lamost.org/dr11/v2.0/doc/lr-data-production-description
- Gaia archive: https://gea.esac.esa.int/archive/


STEP-BY-STEP CODE WALKTHROUGH
================================

The excerpts below are taken directly from run_pipeline.py, except the
explicitly labeled archive-exploration example. They explain the tested
pipeline, not a new implementation. Variables are shared across steps;
individual excerpts are not standalone programs. Run run_pipeline.py to
execute the complete workflow in its implemented order.

The explanation separates discovery from the reproducible run: archive
exploration selected the M67 field and product once; the saved script
reuses those fixed choices. This is not an automatic all-sky matcher.

RUN THE COMPLETE TEST ON THIS MAC

cd /Users/mac/Documents/Codex/2026-09-28/new-chat/outputs/hubble_lamost_match
/Users/mac/Documents/Codex/2026-09-28/new-chat/work/venv/bin/python run_pipeline.py
open diagnostic.png

Expected output: status PASS and 11 true checks.

SETUP AND SHARED HELPERS

ROOT points to the script directory. fetch() downloads only missing files,
using a temporary file before renaming it. save() writes strict JSON.
The required libraries and fixed archive queries are included below.

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


STEP 1. QUERY LAMOST AND READ THE CATALOG

The cone is centered at RA=132.825 degrees, Dec=11.8 degrees, with a
0.05-degree radius. The VOTable is read using its human-readable column
names. The saved response contains 22 observation rows for three distinct
Gaia-linked stars; repeated observations are not separate stars.

    cat_path=fetch('lamost_cone.xml',CAT_URL)
    cat=parse_single_table(cat_path).to_table(use_names_over_ids=True)


STEP 2. FIND THE OVERLAPPING HUBBLE IMAGE

During archive exploration, I queried MAST with the following service and
parameters. This is the discovery code, not part of run_pipeline.py.
The product listing for observation 26095712 contained icmp09011_drz.fits.
The selected image is WFC3/UVIS F621M, 2014-11-09, 10 seconds exposure.
A short exposure does not itself prove the star is unsaturated.

    import requests, json

    response = requests.post(
        'https://mast.stsci.edu/api/v0/invoke',
        data={'request': json.dumps({
            'service': 'Mast.Caom.Cone',
            'params': {'ra': 132.825, 'dec': 11.8, 'radius': 0.06},
            'format': 'json', 'pagesize': 2000
        })}, timeout=45
    )
    response.raise_for_status()
    observations = response.json()

    response = requests.post(
        'https://mast.stsci.edu/api/v0/invoke',
        data={'request': json.dumps({
            'service': 'Mast.Caom.Products',
            'params': {'obsid': '26095712'},
            'format': 'json'
        })}, timeout=45
    )
    response.raise_for_status()
    products = response.json()
    for product in products.get('data', []):
        if product.get('productSubGroupDescription') in ['DRC', 'DRZ']:
            print(product['productFilename'], product['dataURI'])

The saved pipeline downloads the selected image using:

    hst_path=fetch('icmp09011_drz.fits',HST_URL)


STEP 3. FILTER THE LAMOST OBSERVATIONS

Apply the actual demonstration cuts: fiber mask zero, S/N(g)>20,
temperature greater than zero, and [Fe/H]>-5. Sort by S/N(g). One
observation passes in this cached sample: 420808077. These simple cuts
are not a full scientific quality assessment or a general validity range.

    clean=cat[(cat['fibermask']==0)&(cat['snrg']>20)&(cat['teff']>0)&(cat['feh']>-5)]
    assert len(clean)>0
    clean.sort('snrg',reverse=True)
    row=clean[0]
    obsid=int(row['obsid'])


STEP 4. READ THE HUBBLE IMAGE AND DETECT SOURCES

Read SCI and WHT, estimate local background/RMS, and run DAOStarFinder.
The normalization divides the residual image by estimated local RMS.
The threshold of 5 is a detection setting, not a calibrated false-alarm
probability. This run returns 526 candidate detections, including artifacts.
The final lines convert zero-based pixel centroids to RA/Dec using WCS.

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


STEP 5. RETRIEVE GAIA ASTROMETRY AND PROPAGATE TO THE HUBBLE EPOCH

The Gaia query uses the three source IDs discovered in the LAMOST cone
result. The selected row is joined by Gaia source ID. Proper motion
propagates the position from epoch 2016.0 to Hubble EXPSTART.
Radial velocity is assumed zero for this short propagation.

    gaia_path=fetch('gaia.json','https://gea.esac.esa.int/tap-server/tap/sync',
        {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'json','QUERY':GAIA_QUERY})

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


STEP 6. MATCH THE DETECTED SOURCE AND MEASURE APERTURE BRIGHTNESS

Select the closest detected centroid. Matching acceptance is checked in
Step 9, not by nearest-neighbor selection alone. The Gaia-to-Hubble
separation is 0.0399 arcseconds; the next closest detection is 6.88
arcseconds away. Brightness is a preliminary 3-pixel-radius aperture sum.
No aperture correction or precision photometry claim is made.

    sep=sky.separation(at_hst).arcsec
    order=np.argsort(sep);k=int(order[0])
    rawsep=sky[k].separation(SkyCoord(float(row['ra'])*u.deg,float(row['dec'])*u.deg)).arcsec
    aperture=CircularAperture((x[k],y[k]),r=3)
    flux=float(aperture_photometry(resid,aperture,mask=mask)['aperture_sum'][0])
    neighbors=sky.separation(sky[k]).arcsec
    neighbor_count=int(np.sum((neighbors>0.001)&(neighbors<1.65)))


STEP 7. DOWNLOAD AND READ THE LAMOST SPECTRUM

The spectrum is downloaded by observation ID; its primary header and
COADD table are read. Stray NAXIS cards are fixed in memory only;
the original downloaded file is unchanged. Valid spectral samples
require finite wavelength/flux, positive inverse variance and ORMASK=0.
The display flux is divided by its median, not continuum normalized.

    spectrum_url=f'https://www.lamost.org/dr11/v2.0/spectrum/fits/{obsid}'
    spec_path=fetch(f'lamost_{obsid}.fits.gz',spectrum_url)

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


STEP 8. ATTACH THE LAMOST STELLAR PARAMETERS

This builds the joined record. Temperature, [Fe/H], and log g come
directly from the LAMOST catalog row. We did not fit new parameters
from the spectrum or infer them from Hubble pixels. Catalog errors
do not include every systematic uncertainty.

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


STEP 9. RUN THE ELEVEN CONSISTENCY CHECKS

Check the HST identity, positional uniqueness, local detected neighbors,
LAMOST fiber mask, spectrum identity, usable spectral samples, finite
parameters and WCS numerical roundtrip. These engineering checks do
not demonstrate absence of unresolved companions or habitability.

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


STEP 10. EXPORT JSON, SPECTRAL ARRAYS AND PROVENANCE

Only the validated example receives LAMOST parameters; other detections
keep null values. Save raw spectral arrays, all checks, settings, source
URLs and SHA-256 hashes. The downloadable JSON is the pipeline output
that a later model can consume after feature/label separation.

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


STEP 11. PLOT THE HUBBLE FIELD, MATCH AND SPECTRUM

The plot shows the full field, a close-up comparing the detection and
epoch-corrected Gaia position, and the downloaded LAMOST spectrum.
The last assertion fails the run if any consistency check is false.

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


COMPLETE EXECUTABLE SOURCE

Use run_pipeline.py in this same folder. No RL model or transformer is
trained in this demonstration. The full source is kept separate so the
plain-text README remains an explanation rather than a second program.
