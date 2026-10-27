RETRIEVE SPECTRA AND STELLAR PROPERTIES FOR EXISTING MATCHES

Result: 1 saved association processed successfully; 525 detections without a
LAMOST observation ID skipped. No new matches were made. Cached original
LAMOST files were reused, so this run required no network request.

Outputs:
  stellar_properties.json — one record per source/observation association,
    with temperature, [Fe/H], log g, radial velocity, reported uncertainties,
    spectral class/subclass, S/N, quality flags, original match and provenance.
  spectrum_420808077.npz — wavelength, original flux, flux error derived from
    inverse variance, inverse variance, ANDMASK, ORMASK, catalog normalization,
    and a usable-bin mask. Flux remains on the original LAMOST scale.
  retrieval_manifest.json — requested/successful/failed counts and input paths.

Run on this Mac:
cd /Users/mac/Documents/Codex/2026-09-28/new-chat/outputs/hubble_lamost_match
/Users/mac/Documents/Codex/2026-09-28/new-chat/work/venv/bin/python retrieve_matched_stars.py

Process another already-matched sample:
python retrieve_matched_stars.py --matches /path/to/matches.json --catalog /path/to/lamost_catalog.xml --output /path/to/results

Required match record:
  {"source_id": "your_star_id", "lamost_obsid": 420808077}

A JSON list, a dictionary containing a sources list, or one matched record
is accepted. hst_source_id can substitute for source_id. The supplied
catalog must be a LAMOST DR11 v2.0 LRS VOTable containing the requested
observation IDs and named columns. Other releases/formats need adaptation.
This step assumes you already checked the associations. It does not create
new matches or silently assign data to unmatched image detections.

Missing spectra download automatically; existing spectra are reused and
checked for observation-ID consistency. Missing/duplicate catalog rows or
bad spectra are logged as failures while other records continue. A run with
failures exits with status 1. Repeated observations stay separate; no spectra
or stellar measurements are averaged. Outputs with matching names overwrite
previous outputs. Inspect the manifest rather than counting old files.

The code includes explanatory comments beside each function. The annotated
M67 notebook also has a retrieval section at the end.

Scientific limits:
These are existing LAMOST parameter estimates, not parameters newly fitted
from the spectra. Catalog errors do not include all systematic effects.
Quality flags are retained; retrieval alone does not establish scientific
suitability or habitability. This is retrieval, not RL training.
