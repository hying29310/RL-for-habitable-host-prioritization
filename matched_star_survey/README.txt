HUBBLE + LAMOST + GAIA: EXPANDED MATCHED-STAR SAMPLE

WHAT THIS DOES
Finds stars outside M67 that have both a measured Hubble source position and
LAMOST stellar spectroscopy, connected through a Gaia DR3 source identifier.
It downloads the spectrum and a real Hubble image cutout for each accepted star.

FINAL RESULT
31 distinct matched stars outside M67, in seven fields:
  NGC 6819: 16    NGC 6791: 5    M3: 4    M53: 2
  NGC 2420: 2     NGC 2158: 1    M15: 1
With the three previous M67 stars: 34 distinct stars in total.
One NGC 6791 star was already known from the earlier pilot, so this expansion
adds 30 previously unreported stars to our existing sample of four.

Each accepted star has:
- A Gaia DR3 identifier, position, proper motion, G magnitude and RUWE.
- An exact LAMOST DR11 v2.0 LRS observation identifier.
- Temperature, metallicity [Fe/H], surface gravity, radial velocity and
  catalog errors, plus spectral class/subclass and signal-to-noise ratios.
- The original LAMOST FITS spectrum and exported numerical spectrum arrays.
- Hubble Source Catalog v3 detection records, image names, dates and flags.
- A downloaded Hubble Legacy Archive FITS image cutout for visual review.

IMPORTANT DIFFERENCE FROM THE FIRST PILOT
The original M67 and two-image pilot ran DAOStarFinder on downloaded images.
This larger search uses published Hubble Source Catalog (HSC) detections.
We check those measured positions and download image cutouts; we do NOT claim
that all 31 stars passed our earlier DAOStarFinder image-extraction pipeline.
The HSC includes other images and visits. That is why NGC 6819 now has matches
even though the two particular images in the first pilot had ambiguous peaks.

HOW A STAR ENTERS THE FINAL LIST
1. Search LAMOST: retain stars with valid Gaia IDs, finite temperature,
   metallicity and gravity, temperature >0, [Fe/H]>-5 and g-band SNR >20.
   Choose one spectrum per star: prefer fiber mask zero, then higher SNR.
2. Find potential HSC objects within 3 arcseconds of the LAMOST coordinates.
   This is only a discovery radius, NOT the final acceptance tolerance.
3. Retrieve the Gaia source by its exact identifier; do not convert IDs to
   floating-point numbers. Propagate its proper motion to each Hubble date.
4. Examine HSC objects with at least two images, an absolute astrometric
   correction and a summary position within 0.6 arcseconds of predicted Gaia.
5. Require detections in at least two distinct filter images within 0.3
   arcseconds of predicted Gaia, and a median detection offset below 0.3.
   Exactly one tested HSC object must pass for each Gaia star. An HSC object
   cannot be assigned to multiple Gaia stars.
6. Query nearby Gaia sources. The assigned source must be closest to the
   selected Hubble detection, with no other Gaia source within 0.6 arcseconds.
7. Also compare LAMOST coordinates with Gaia at the LAMOST observation date:
   the assigned Gaia star must be nearest and within 1.5 arcseconds.
8. Download spectra by exact observation ID, verify FITS IDs and fiber masks,
   preserve wavelength/flux/inverse variance/masks, and export catalog errors.
9. Count unique Gaia IDs, never repeated images or repeated spectra.

34 stars passed the Hubble/Gaia stage. The extra LAMOST-position check rejected
three: two in M53 and one in NGC 2158 had another Gaia star closer to the
LAMOST coordinates. They are not included in the final 31. See neighbor_audit.json.

QUALITY AND SCIENTIFIC LIMITS
- 17 of the 31 selected LAMOST spectra have fiber mask zero; 14 have flags.
- Six stars have at least one HSC detection with flag zero. Only three have
  both a zero LAMOST fiber mask and a zero-flag HSC detection.
- Most Hubble detections have flags, often for saturation/extended shape.
  They can support a positional association but their brightness measurements
  must not be treated as reliable without further work.
- Several filter measurements may share the same white-light source detection.
  The field distinct_white_light_detections records this; two filter images
  do not necessarily mean two independent detections or two separate visits.
- Zero flags are not a guarantee of accuracy, purity, or absence of blending.
  HSC candidate preselection and Gaia incompleteness can miss close companions.
  Missing-proper-motion neighbors are checked at their catalog positions.
- Proper-motion propagation ignores radial velocity/perspective acceleration.
  No probabilistic association confidence or false-match rate was estimated.
- The cutouts are for visual review. HLA astrometric corrections are requested
  using ApplyOmega; cutout centers are not newly fitted source centroids.
- These are catalog matches, not proof of planets or habitability. Field names
  do not establish cluster membership. Catalog errors omit some systematics.
- This sample is selected by available spectra/images and quality cuts. It is
  not representative of all stars and is not an all-sky completeness estimate.

SEARCH AREA
The reproducible default search uses 0.25-degree-radius cones around NGC 6791,
NGC 6819 and NGC 2420, and 0.08-degree-radius cones around M13, M3, M92, M15,
NGC 2158, NGC 2264, Orion, NGC 5466, M53 and NGC 6229 (13 regions).
search_plan.json gives exact coordinates. Other exploratory LAMOST cones and
old pilot inputs are retained in data/; not every exploratory cone completed
an HSC search. --wide requests broader searches of all 24 exploratory regions.

FILES TO OPEN FIRST
RESULTS.txt                      Simple results and one-line list of stars.
star_catalog.json                Final 31 stars, measurements and file links.
summary.json                     Verified counts and quality breakdown.
hubble_contact_sheet.png         Actual Hubble cutouts for all 31 stars.
spectra_contact_sheet.png        Actual downloaded LAMOST spectra.
all_fields_including_M67.json     34 unique stars including the old M67 matches.
matched_stars_workflow.ipynb      Notebook with explanations and runnable steps.

CODE AND WHAT EACH FILE DOES
search_fields.py      Downloads LAMOST cone catalogs; searches HSC candidates.
validate_matches.py   Checks Gaia motion and HSC individual image detections.
check_neighbors.py    Checks competing Gaia stars and the LAMOST positions.
download_evidence.py  Downloads spectra and Hubble FITS cutouts for accepted stars.
retrieve_matched_stars.py  Reads FITS spectral arrays and catalog properties.
make_report.py        Verifies counts and makes the final catalog and previews.

RUN THE PIPELINE FROM THE PROJECT DIRECTORY
work/venv/bin/python outputs/matched_star_survey/search_fields.py
work/venv/bin/python outputs/matched_star_survey/validate_matches.py
work/venv/bin/python outputs/matched_star_survey/check_neighbors.py
work/venv/bin/python outputs/matched_star_survey/download_evidence.py
work/venv/bin/python outputs/matched_star_survey/make_report.py

Run these sequentially, not concurrently. Existing cached responses are reused.
Internet is required for missing inputs. A different archive response or newly
added search region can change the sample. make_report.py checks that at least
30 stars have all required outputs; it does not pad the count to meet a target.
Dependencies: Python 3, numpy, scipy, astropy, photutils (earlier pilot only),
matplotlib, requests. See requirements.txt for the versions used here.

AUDIT TRAIL
preliminary_matches.json: HSC/Gaia candidates before LAMOST-position rejection.
matching_audit.json: HSC tests, offsets and rejected candidate reasons.
neighbor_audit.json: nearby Gaia stars and final positional decisions.
selected_lamost.xml: catalog rows for preliminary matches (a superset of final).
retrieved/stellar_properties.json: spectra/parameters for final accepted stars.
retrieved/retrieval_manifest.json: retrieval outcome and input catalog checksum.
cutout_manifest.json: Hubble download URLs and SHA256 hashes.
data/hsc_details/: original per-image HSC measurements.
data/gaia_neighbors/: original Gaia neighborhood query responses.
data/spectra/: original spectra. Extra files for rejected preliminary candidates
may remain cached; only IDs in matches.json count toward the final sample.

OFFICIAL REFERENCES
Hubble catalog API:
https://catalogs.mast.stsci.edu/docs/hsc.html
Hubble cutout service (including ApplyOmega):
https://hla.stsci.edu/fitscutcgi_interface.html
HSC v3 documentation:
https://outerspace.stsci.edu/spaces/HSC/pages/421462969/HSCv3+FAQ+Hubble+Source+Catalog
LAMOST DR11 v2.0:
https://www.lamost.org/dr11/v2.0/
Gaia archive:
https://gea.esac.esa.int/archive/

Based on observations made with the NASA/ESA Hubble Space Telescope and obtained
from the Hubble Legacy Archive. Parameters and spectra are from LAMOST DR11 v2.0;
astrometry is from Gaia DR3.
