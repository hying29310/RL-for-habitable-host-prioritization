# Kepler image-to-transit pipeline: successful smoke test

This is an executed engineering test using real public data, not an RL training result.

## What ran

1. Downloaded the Kepler Quarter 1 target-pixel sequence for **Kepler-8 / KIC 6922244**: 1,639 images, each 6 × 6 pixels.
2. Selected the six pixels in the archived optimal aperture and summed their background-subtracted flux. Added pixel errors in quadrature.
3. Retained 1,452 finite, positive-flux observations with QUALITY = 0, spanning 33.47 days.
4. Compared the extracted flux with the official Kepler SAP light curve: maximum relative difference **1.03 × 10⁻⁷**.
5. Removed slow trends with a segment-wise running median, and searched periods from 0.5–10 days using Astropy box least squares (BLS). The catalog period and epoch were not passed to the search.
6. Queried the NASA Exoplanet Archive for stellar measurements. Verified the named target's coordinates against the FITS metadata (separation **0.0725 arcsec**). This verifies a known identity; it is not a general catalog cross-matcher.
7. Ran a scripted `coarse → refine → submit` episode through a minimal RL-style environment. Checked budget enforcement and rejection of post-terminal actions.
8. Exported machine-readable results, the extracted light curve, a median image, and a diagnostic plot.

## Results

| Quantity | Result |
|---|---:|
| Recovered period | 3.521820 days |
| Published comparison period | 3.5224991 days |
| Relative period error | 0.0193% |
| Mid-transit epoch difference, modulo period | 2.50 minutes |
| Approximate box-fit depth | 0.792% |
| Catalog stellar effective temperature | 6,213 ± 150 K |
| Catalog [Fe/H] | −0.055 ± 0.030 dex |
| Automated checks | 7 / 7 passed |

The median image is a derived stack of the downloaded image sequence. The JSON uses catalog temperature and metallicity; those quantities were not inferred from the image. The reference for the stellar and planet parameters is preserved in `data/catalog.json` and `results.json`. Uncertainties are those reported by the catalog, not a new uncertainty analysis.

## Environment test

| Action | Abstract cost | Remaining budget | Reward |
|---|---:|---:|---:|
| coarse | 2 | 4 | −0.10 |
| refine | 3 | 1 | −0.15 |
| submit | 1 | 0 | +0.95 |

Initial budget = 6 units. Reward = one period-and-epoch match bonus minus 0.05 per cost unit; total = 0.70. These are illustrative accounting units, **not calibrated runtime costs**. Actual action wall times are recorded separately. Submission requires relative period error < 0.1% and epoch error < 0.1 day. This is a recovery sanity check, not full planet vetting. A failed match is not automatically a false positive.

There is **no trained policy, GRPO optimization, multi-star selection, habitable-zone scoring, or RL performance comparison** in this test. Actions are scripted. Kepler-8 b is a known giant planet (catalog radius about 15.87 Earth radii), selected to make this first test easy to verify. It is not a small habitable-zone test case. Refining a box fit improves that fit's objective but need not reduce its difference from a published physical-model period.

## Reproduce

Use Python 3.12. From this directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python run_pipeline.py
```

Raw inputs are included. Cached reruns require no network; missing inputs are downloaded from the recorded public endpoints. The script overwrites derived outputs. `requirements-lock.txt` records the exact package versions used in this execution; it may need adjustment on other platforms.

Main files:

- `run_pipeline.py`: download, extraction, transit search, and environment checks.
- `results.json`: measurements, provenance, checks, and action/reward trace.
- `diagnostic.png`: image, light curve, periodogram, folded transit.
- `extracted_lightcurve.npz`: numerical arrays for further analysis.
- `median_image.fits`: median image in local cutout coordinates; no sky WCS is attached.
- `data/`: original public FITS inputs and the archived catalog response.

## What this establishes, and what comes next

**Established:** real astronomical pixels can pass through calibrated flux extraction, a numerical transit search, JSON export, and a budgeted action interface.

**Still needed for the thesis:** a multi-star Kepler/LAMOST/Gaia sample; valid uncertainty and missing-data treatment; train/validation/test splits by star; transit-injection controls; false-positive vetting; a justified habitable-zone target definition; realistic action costs; baselines; and actual policy training. A short 33-day, 0.5–10-day-period example does not establish sensitivity to long-period habitable-zone signals. Hubble and TESS ingestion were not tested here.

For learned policies, remove target identifiers and bibliographic reference strings from their feature observations and keep labels only in a separate evaluator. The demo preserves provenance for auditing, and no policy is trained on those fields.

## Sources

- [Kepler download documentation](https://archive.stsci.edu/kepler/download_options.html)
- [Target-pixel archive directory](https://archive.stsci.edu/missions/kepler/target_pixel_files/0069/006922244/)
- [Astropy BLS documentation](https://docs.astropy.org/en/stable/timeseries/bls.html)
- [NASA Exoplanet Archive TAP documentation](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html)

Exact download URLs, query text, file sizes, SHA-256 hashes, and software versions are included in the saved results.
