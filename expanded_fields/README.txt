EXPANDED HUBBLE / LAMOST / GAIA FIELD TEST
Run date: 2026-10-25 UTC

RESULT
One new positionally validated star outside M67, in the NGC 6791 field.
Together with the three previous M67 matches, four distinct Gaia sources.
These are positional associations, not evidence of planets or habitability.
Field names identify image regions, not confirmed cluster membership.

Field        Catalog candidates*   New accepted image matches
NGC 6791     12                    1
NGC 6819     37                    0 (three or more ambiguous detections)
NGC 2420     24                    Not processed: no eligible star in the
                                  returned modern-camera footprint examined
NGC 188       0                    Not processed
*Unique Gaia IDs with LAMOST parameter entries in radius 0.08 degree cones.
These counts are not numbers of Hubble matches or all-sky overlap estimates.
After STAR, SNR(g)>20, finite temperature/[Fe/H] filtering, NGC 6791 had
10 eligible sources and NGC 6819 had 32. Some lie outside actual pixels.
NGC 188's zero applies only to this small query, not the whole cluster/sky.

NEW STAR
Gaia DR3: 2051105376456445824
LAMOST DR11 v2.0 LRS observation: 664016211
Hubble: ibc603zoq_drc.fits, WFC3/UVIS, F606W, 30 s, 2010-09-02
Gaia-motion-corrected separation: 0.00836 arcsec
Temperature: 7150.92 +/- 47.77 K (catalog statistical error)
[Fe/H]: +0.430 +/- 0.041 dex; log(g): 3.889 +/- 0.068 dex
LAMOST spectral subclass: F0; SNR(g): 54.93
Spectrum: 3901 wavelength samples, 3699.99--9082.39 Angstrom.
IMPORTANT: LAMOST fiber mask is 192, not zero. Parameters/spectrum require
quality review. Four other image detections within 1.65 arcsec also warrant
blending review; some detections can be artifacts. HST saturation is not
fully audited. No new zero-fiber-mask star was added in this test.

PROCEDURE AND CODE
1. discover_fields.py: queries LAMOST cone search and MAST for four centers.
   Cached raw results and ranked-image lists are in data/<field>/.
   Footprint ranking was exploratory; chosen image products are explicit
   in acquire_fields.py and process_fields.py, not automatically optimized.
2. acquire_fields.py: downloads selected calibrated HST images, LAMOST
   VOTables and Gaia DR3 proper-motion measurements, reusing cached files.
   Network is needed for missing inputs; two selected images total 434 MB.
3. process_fields.py, at_epoch(): propagates Gaia positions to Hubble dates.
4. process_fields.py, process(): checks image bounds and weight coverage;
   extracts 129 x 129 pixel cutouts near catalog positions; detects sources
   using DAOStarFinder at 5 times the local clipped noise and FWHM 2.5 pixels.
   Accepts only one detected centroid within 0.3 arcsec of propagated Gaia.
   These are targeted detections, not a blind-source completeness survey.
   LAMOST-to-Gaia associations use the Gaia identifier supplied by LAMOST.
5. retrieve_matched_stars.py: selects catalog entries by exact observation
   ID, downloads FITS spectra, validates the observation ID, and exports
   wavelength, flux, inverse variance, uncertainty and quality-mask arrays.
   Parameters come from LAMOST; they are not inferred from Hubble images.
6. process_fields.py, main(): collects match, property and audit JSON files;
   makes full-field and matched-source PNG previews. Failed matches are
   retained in audit.json, not silently accepted.

NGC 6819: tried iexu17quq_drc.fits (F350LP, 40 s) and ib2o06020_drc.fits
(F814W, 10 s). Neither passed the unique-centroid criterion. The first
attempt is preserved in NGC6819_initial_40s; the second is in NGC6819.
Multiple peaks can arise from crowding, saturation or detection artifacts;
these causes have not been distinguished. Zero accepted is not proof of
zero genuine overlap. More images or improved PSF/deblending can be tested.

RUN FROM THE PROJECT DIRECTORY
work/venv/bin/python outputs/expanded_fields/acquire_fields.py
work/venv/bin/python outputs/expanded_fields/process_fields.py

Dependencies: numpy, scipy, astropy, photutils, matplotlib, requests.
The optional discovery script repeats remote archive queries; rerunning
processing uses local inputs and retrieves only missing spectra.

OUTPUTS
all_new_matches.json: new accepted source associations.
all_new_stellar_properties.json: parameters, errors, quality and spectra links.
combined_matches_including_M67.json: four distinct sources, including old M67.
expansion_summary.json: counts, image metadata, SHA256 and retrieval statuses.
NGC6791/retrieved/: stellar_properties.json, NPZ spectrum, retrieval manifest.
data/NGC6791/: original LAMOST FITS spectrum and catalog; Hubble/Gaia inputs.
NGC6791/field.png and matches.png: new field and accepted-source cutout.
NGC6819/field.png: second new field, with no accepted match marked.

This pilot samples selected images. It does not exhaust all Hubble images,
LAMOST releases, search radii, or stars without catalog stellar parameters.
