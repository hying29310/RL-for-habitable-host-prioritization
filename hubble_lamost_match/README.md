# Hubble → LAMOST: successful real-data match

This test downloaded a real Hubble image of M67, detected sources, matched one source using Gaia astrometry, downloaded its actual LAMOST DR11 v2.0 spectrum, and attached LAMOST temperature and metallicity measurements. **All 11 checks passed.** No stellar parameters were inferred from the Hubble image.

## Data and result

| Item | Value |
|---|---|
| Hubble product | `icmp09011_drz.fits` |
| Instrument / filter | WFC3/UVIS / F621M |
| Hubble observation | 2014-11-09; 10 seconds total exposure |
| Image dimensions | 2,063 × 2,192 pixels |
| LAMOST star | J085119.17+114754.6 |
| LAMOST observation ID | 420808077 |
| Gaia DR3 source ID | 604917526275831040 |
| LAMOST observation | 2016-02-06 |
| Hubble–Gaia separation at Hubble epoch | 0.0399 arcseconds |
| Hubble–LAMOST catalog separation | 0.0683 arcseconds |
| Effective temperature | 6105.19 ± 16.18 K |
| Metallicity [Fe/H] | −0.145 ± 0.010 dex |
| Surface gravity log g | 4.241 ± 0.020 dex |
| LAMOST spectral subclass | F6 |
| g-band spectral S/N | 217.48 |
| LAMOST fiber mask | 0 |
| Downloaded spectral range | Approximately 3700–9099 Å |

**Uncertainties above are the catalog's reported errors, not a complete uncertainty budget.** Survey calibration, stellar modeling, and other systematic errors can be larger.

## How the field and star were selected

The LAMOST public cone search at RA=132.825°, Dec=11.8°, radius=0.05° returned 22 observation rows representing three Gaia-linked stars. MAST returned overlapping M67 Hubble images. The selected short-exposure F621M image includes all three stars. Short exposure reduces, but does not establish absence of, saturation; a full HST detector-quality audit was not performed.

The example applies a conservative filter: fiber mask zero, S/N(g)>20, and available plausible temperature and metallicity. One observation passed. The other two stars' available candidate records have nonzero fiber flags; two inspected examples have `fibermask=64`, the documented WHOPPER flag. They are not promoted as clean matches here. The original cone response is included for audit.

## Processing and validation

1. Read Hubble SCI, WHT, and celestial WCS; mask invalid or uncovered pixels.
2. Estimate local background and RMS with 64-pixel boxes.
3. Run DAOStarFinder (FWHM=2.5 pixels, threshold=5 on RMS-normalized residuals), yielding 526 candidate detections. Many can be artifacts; this is not a scientifically vetted star catalog.
4. Retrieve Gaia DR3 astrometry using the Gaia identifier already present in LAMOST. Propagate its position from epoch 2016.0 to the Hubble exposure epoch using proper motion. Radial velocity is set to zero for this short propagation; the positive Gaia parallax is used only for space-motion propagation, not a scientific distance estimate.
5. Find the nearest independently detected Hubble source. Exactly one detection is within 0.3 arcsec; the second nearest is 6.88 arcsec away. No additional pipeline detections lie within 1.65 arcsec of the matched source. Undetected or unresolved companions remain possible.
6. Download the LAMOST spectrum and verify the observation ID and coordinates against the catalog. The file has 3,909 usable spectral bins under the applied finite-flux, positive-inverse-variance, and zero-ORMASK criteria.
7. Save the association, parameters, detected sources, spectral arrays, checks, and plot. All other detections retain null stellar parameters.

The LAMOST FITS primary header contains stray NAXIS1/NAXIS2 cards for an empty primary HDU. Astropy fixes these in memory before verification. The original downloaded bytes are preserved and hashed. Spectral data are read from the COADD table. The plotted spectrum is divided by its median for display and is not presented as absolute calibrated flux.

## Run

On this workspace:

```bash
cd /Users/mac/Documents/Codex/2026-09-28/new-chat/outputs/hubble_lamost_match
/Users/mac/Documents/Codex/2026-09-28/new-chat/work/venv/bin/python run_pipeline.py
open diagnostic.png
```

Elsewhere, with Python 3.12 installed:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python run_pipeline.py
```

The ZIP includes the four required original inputs, so cached reruns are offline. Missing inputs download automatically. This is a reproducible demonstration for the selected field, not a general all-sky search program: the M67 image and Gaia query are fixed to the field chosen during archive exploration.

## Outputs

- `matched_star.json`: the verified example and catalog parameters.
- `sources.json`: 526 image detections, one joined to clean LAMOST data.
- `results.json`: checks, URLs, settings, timestamps, and SHA-256 provenance.
- `spectrum_arrays.npz`: wavelengths, original flux, inverse variance, pixel mask, and usable-bin mask.
- `diagnostic.png`: Hubble field, close-up of the positional match, and LAMOST spectrum.
- `run_pipeline.py`: standalone processing code.
- `data/`: Hubble image, LAMOST cone result, Gaia astrometry, and matched LAMOST spectrum.

This establishes a working image-to-spectroscopic-catalog connection. It does not establish a planet, habitability, membership of M67, a full precision-photometry solution, or performance of an RL policy.

## Direct sources

- [Hubble FITS](https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/icmp09011_drz.fits)
- [LAMOST spectrum record](https://www.lamost.org/dr11/v2.0/spectrum/view?obsid=420808077)
- [LAMOST spectrum FITS](https://www.lamost.org/dr11/v2.0/spectrum/fits/420808077)
- [LAMOST DR11 data definitions and masks](https://www.lamost.org/dr11/v2.0/doc/lr-data-production-description)
- [Gaia archive](https://gea.esac.esa.int/archive/)
