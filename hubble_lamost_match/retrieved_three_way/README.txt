THREE-WAY HUBBLE / LAMOST / GAIA MATCHES IN THE TESTED M67 FIELD

3 distinct stars were positionally matched; spectroscopy and stellar
properties were extracted successfully for all 3. These are the stars in
the tested LAMOST cone and Hubble image, not an archive-wide count.

We propagated each Gaia DR3 position to the Hubble exposure epoch and
required exactly one image detection within 0.3 arcseconds. All 3 passed,
with offsets of 0.024 to 0.049 arcseconds. No two stars share a detection.

One spectrum per star was selected: prefer a zero fiber mask, then the
highest g-band signal-to-noise among available rows with usable parameters.

LAMOST OBSID   Teff (K)          [Fe/H] (dex)      Fiber mask
420808077      6105.19 +/-16.18  -0.145 +/-0.010   0
804808069      6629.59 +/-32.41   0.016 +/-0.027   0
1059308069     5650.34 +/-25.81   0.138 +/-0.021   64 (WHOPPER)

Consult stellar_properties.json for the exact catalog uncertainties and
all other properties. These errors do not include all systematic effects.
The mask=64 spectrum was retained for inspection, not declared clean.
Positional matching and spectral quality are separate judgments.

Correction to earlier reporting: the initial pipeline selected ONE
highest-S/N observation after filtering. That did not mean only one
star had any unflagged spectrum. The second star has 20 LAMOST observation
rows, including unflagged alternatives. We now selected observation
804808069 for that star. There are 22 observation rows for 3 unique stars;
this run extracts one representative observation per star, not all epochs.

Outputs:
- stellar_properties.json: all 3 stars, catalog properties and uncertainties,
  original match evidence, quality flags, original FITS paths and hashes.
- spectrum_<obsid>.npz: wavelength, flux, error, inverse variance, masks,
  normalization and usable-bin indicator for each star.
- retrieval_manifest.json: 3 retrieved, 0 failures.

Reproduce from the parent project folder:
python match_all_three_catalogs.py

Inputs and scripts are in the parent hubble_lamost_match folder.
all_three_way_matches.json stores the 3 accepted associations.
three_way_match_audit.json records nearest and second-nearest separations.
sources_three_way.json enriches the original detection list with these
3 matches while preserving the original one-example sources.json.

No new metallicity fits, habitable-planet claims, or RL training are made.
