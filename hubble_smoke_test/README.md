# Hubble image pipeline test — PASS

Executed on a real public Hubble science image on 2026-10-04. Eight engineering checks passed. This validates data ingestion, candidate detection, pixel-to-sky conversion, aperture summation, and JSON export; it does not validate a planet-search or star-classification model.

## Observation

- Telescope/instrument: HST WFC3/UVIS.
- FITS product: `icdm0a070_drc.fits` (CTE-corrected drizzled image).
- FITS target name: `NGC1512-2`.
- Filter: F814W. Observation date: 2013-09-02. Exposure time: 977 seconds.
- Original dimensions: 4,181 columns × 4,462 rows; approximately 224 MB.
- Analyzed central crop: 1,024 × 1,024 pixels.
- [Original public FITS download](https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/icdm0a070_drc.fits)
- This filename is used in [STScI's Hubble astrometry example](https://hst-docs.stsci.edu/hstdhb/4-hst-data-analysis/4-6-analyzing-hst-images).

## Executed workflow

1. Download the scientific FITS image and verify its HST header and FITS structure.
2. Read the SCI pixels, WHT coverage, and celestial WCS.
3. Extract a fixed central crop, preserving its sky-coordinate mapping.
4. Mask missing coverage and estimate a spatially varying background using 64-pixel boxes.
5. Run Photutils DAOStarFinder on a local-background-RMS-normalized image, with FWHM 2.5 pixels and threshold 5.
6. Measure background-subtracted 3-pixel-radius aperture sums in electrons/second.
7. Convert detected centroids to RA/Dec using the image WCS.
8. Export source records to JSON and save an annotated image and FITS crop.

## Results

**183 candidate point sources** were returned. All eight checks passed: HST identity, celestial WCS, valid-pixel finiteness, nonempty detections, finite sky coordinates, coordinate round-trip accuracy, finite aperture sums, and preservation of WCS in the saved FITS crop. Maximum pixel → sky → pixel discrepancy was approximately **2.7 × 10⁻⁹ pixels**. This is numerical consistency, not a claim about absolute sky-position accuracy.

The annotated image shows detections around a bright source as well as isolated sources. Bright-star artifacts, blends, star clusters, and compact galaxies can contaminate this initial detection list. **183 detections does not mean 183 confirmed stars.** No completeness or false-positive rate was measured. The detection settings require tuning and validation before scientific use.

## What remains unknown

Every source has `temperature_K: null` and `metallicity_Fe_H_dex: null`, with `catalog_match_status: "not_attempted"`. These are intentionally missing: the pipeline has not retrieved spectra or matched Gaia/LAMOST. The field contains an external galaxy, so many detections are unsuitable for your intended nearby-star planet search. This image was chosen to validate Hubble handling, not to select optimal exoplanet targets.

The aperture sums are preliminary: no total-flux aperture corrections, uncertainty estimates, crowding corrections, or comprehensive artifact rejection were applied. A single-band image supplies neither a transit light curve nor reliable temperature/metallicity by itself. No transformer, RL training, habitability score, or planet detection was performed.

## Run on the existing workspace

```bash
cd /Users/mac/Documents/Codex/2026-09-28/new-chat/outputs/hubble_smoke_test
/Users/mac/Documents/Codex/2026-09-28/new-chat/work/venv/bin/python run_pipeline.py
open diagnostic.png
```

## Run elsewhere

With Python 3.12 installed, in the extracted project directory:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python run_pipeline.py
```

The portable ZIP excludes the 224 MB original image to keep the package small; the script downloads it when missing. The original is already cached in this workspace, so local reruns are offline. The downloaded file's URL and SHA-256 are recorded in `results.json`. Derived outputs are overwritten on reruns.

## Files

- `run_pipeline.py`: full runnable code.
- `sources.json`: all 183 candidate records, including sky positions and aperture fluxes.
- `results.json`: provenance, settings, versions, checks, and limitations.
- `science_cutout.fits`: science crop with WCS and a WHT extension.
- `diagnostic.png`: full field and annotated crop.
- `requirements-lock.txt`: exact installed dependencies.
- `data/icdm0a070_drc.fits`: original cached data, available locally but omitted from the portable ZIP.

For the thesis, the next scientific step is selecting a suitable Milky Way stellar field, filtering artifacts, cross-matching detections against stellar catalogs with epoch-aware astrometry, and joining time-series observations. This test establishes the Hubble image-to-JSON component only.
