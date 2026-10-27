70 ADDITIONAL HUBBLE--LAMOST--GAIA STARS

WHAT THIS DELIVERS
70 new distinct stars, excluding all 34 previously known IDs (31 in the earlier
expanded sample plus 3 original M67 stars). These are catalog associations checked
against HST detection positions and Gaia proper motion, not planet detections.
No light curves are used. Field labels identify search cones; a star in the M31
field is not thereby a member of M31, and cluster membership is not established.

START HERE
RANKED_STARS.txt             Human-readable table with all 70 stars and scores.
ranked_stars.json            Full parameters, errors, quality, spectrum paths,
                            score components and sensitivity ranges.
RANKING_METHOD.txt          Formula, assumptions, interpretation and references.
verification.json          Actual counts and automated verification results.
selected_70_matches.json    HSC detections, image names, epochs and match evidence.
retrieved/                  70 extracted wavelength/flux/error/mask NPZ spectra,
                            plus retrieval manifest and catalog audit.
data/spectra/              Original LAMOST FITS (also contains earlier cached data).
selection_method.json       Exact rule used to select 70 complete stars.

SOURCES AND PROCEDURE
1. search_fields.py / expand_search.py / search_extra.py:
   Query LAMOST DR11 v2.0 LRS cones around HST-covered fields. Require STAR class,
   catalog temperature, metallicity and Gaia ID, and g-band SNR >20. Keep an
   eligible observation per Gaia ID; prefer zero fiber mask and higher SNR.
   Crossmatch positions with Hubble Source Catalog v3 within 3 arcsec.
   https://www.lamost.org/dr11/v2.0/
   https://catalogs.mast.stsci.edu/docs/hsc.html
2. validate_matches.py:
   Remove every previously known Gaia ID. Propagate Gaia DR3 proper motions to
   HST observation times. Require exactly one qualifying HSC association with
   >=2 distinct filter images within 0.3 arcsec and median separation <0.3 arcsec.
   Different filter images are not necessarily independent observing visits.
   Require corrected HSC astrometry. Reject HSC sources assigned to two stars.
3. check_neighbors.py:
   Check Gaia neighbors within 5 arcsec. The intended star must be closest to
   the selected HST detection, with other Gaia stars farther than 0.6 arcsec.
   Also require it to be closest to the LAMOST position at its observation date,
   within 1.5 arcsec. This does not exclude all unresolved blends or binaries.
4. get_stellar_parameters.py:
   Retrieve Gaia DR3 astrophysical_parameters by exact integer source ID.
   Use FLAME radius/luminosity and lower/upper bounds. These are model-derived
   estimates, not direct measurements; preserve flags and systematic caveats.
   https://gea.esac.esa.int/archive/documentation/GDR3/
5. finalize_sample.py / retrieve_matched_stars.py:
   Require zero LAMOST fiber mask, then select 70 complete stars by HST flag
   quality, then LAMOST g-band SNR, then ID.
   Do not select using the resulting host-priority score. Download each spectrum
   by exact LAMOST obsid; verify FITS observation ID and fiber mask. Preserve
   wavelength, catalog flux, inverse variance and masks, deriving flux errors
   as 1/sqrt(inverse variance). Keep original FITS and SHA-256 checksums.
   Temperature, [Fe/H] and spectral subclass are published LAMOST estimates;
   this pipeline does not refit spectra. LAMOST obsid identifies an observation,
   while Gaia DR3 source_id identifies a source in that release.
6. rank_hosts.py:
   Apply the explicit heuristic in RANKING_METHOD.txt. The spectrum supplies
   temperature, metallicity and classification; it is not counted a second time
   as independent evidence. Only the five requested physical input categories
   contribute. Quality metadata are reported separately.
7. export_results.py:
   Verify 70 unique Gaia/HSC IDs, no previous IDs, complete finite inputs,
   valid extracted spectral arrays, scores in [0,100], and basic formula behavior.
   Write plain-text results and verification.json.

RUN AGAIN
Use Python 3.12 with numpy, scipy, astropy, requests, matplotlib and photutils.
From this directory, with that environment activated:
  python expand_search.py --phase all
  python search_extra.py
  python validate_matches.py
  python check_neighbors.py
  python get_stellar_parameters.py
  python finalize_sample.py
  python rank_hosts.py
  python export_results.py
Queries use cached responses when available. Search jobs must finish before
validation for a full rerun; later-arriving cones can increase the candidate pool
and change which 70 are selected. The delivered selected_70_matches.json records
the actual cohort. Rerun just rank_hosts.py and export_results.py to reproduce its
scores offline. JSON paths to spectra are absolute for this workspace; update
those paths if relocating files. Identifiers should be kept as strings, never
rounded by spreadsheet software.

SCIENTIFIC LIMITS
A high score means priority under chosen dwarf-star/temperature preferences;
it is not a measured probability of a planet, habitability or life. A low score
does not rule out planets. The score is a baseline requiring calibration and
validation before any probability claim. Radius and luminosity are correlated;
their quoted bounds do not represent independent errors. Catalog parameter
errors omit some systematics. Flagged HST detections may be saturated or extended.
Incomplete samples and HST/LAMOST coverage create strong selection biases.
Stars absent from this list may simply lack suitable catalog data.
Use sensitivity ranges to see how chosen assumptions change scores. They are
not statistical confidence intervals. No planet signals were searched here.

DELIVERED SAMPLE QUALITY
All 70 selected LAMOST spectra have zero fiber masks. Ten stars have no
zero-flag HST detection; their HST photometry should be reviewed. Two stars
(Gaia DR3 3950102317148387456 and 4275336138666632064) have T/R/L
consistency residuals above 0.3 dex; their scores remain provisional.
The physical estimates are not independent or guaranteed accurate.

Sample overview: sample_overview.png (generated by plot_results.py).
