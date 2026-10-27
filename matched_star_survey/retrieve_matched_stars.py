"""Retrieve LAMOST DR11 v2.0 LRS data for existing matches; no new matching.
Example: python retrieve_matched_stars.py --matches sources.json
Supply a same-release VOTable with --catalog when processing another sample.
"""
from pathlib import Path
import argparse
import hashlib
import json
from datetime import datetime, timezone
import numpy as np
import requests
from astropy.io import fits
from astropy.io.votable import parse_single_table

ROOT = Path(__file__).resolve().parent


def scalar(value):  # Convert a catalog scalar into JSON-safe data, preserving missing values.
    if np.ma.is_masked(value):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path, obj):  # Save readable JSON without nonstandard NaN/Infinity tokens.
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def load_matches(path):  # Accept a source list, a {'sources': [...]} wrapper, or one matched record.
    obj = json.loads(path.read_text())
    records = obj if isinstance(obj, list) else obj.get('sources', [obj])
    selected, seen = [], set()
    for record in records:
        obsid = record.get('lamost_obsid')
        if obsid is None:
            continue
        obsid = int(obsid)
        if obsid <= 0:
            raise ValueError('LAMOST observation IDs must be positive')
        release = record.get('lamost_release', 'DR11 v2.0 LRS')
        if release != 'DR11 v2.0 LRS':
            raise ValueError('This downloader supports DR11 v2.0 LRS only')
        source_id = record.get('source_id', record.get('hst_source_id'))
        if not source_id:
            raise ValueError('Each matched record needs source_id or hst_source_id')
        key = (source_id, obsid)
        if key not in seen:
            seen.add(key)
            selected.append((record, source_id, obsid))
    return records, selected


def fetch_spectrum(obsid, cache):  # Reuse an existing FITS spectrum or download it atomically.
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f'lamost_{obsid}.fits.gz'
    url = f'https://www.lamost.org/dr11/v2.0/spectrum/fits/{obsid}'
    if not path.exists():
        tmp = path.with_suffix('.partial')
        try:
            with requests.get(url, stream=True, timeout=(30, 120)) as r:
                r.raise_for_status()
                with tmp.open('wb') as out:
                    for chunk in r.iter_content(1024 * 1024):
                        out.write(chunk)
            with fits.open(tmp) as hd:
                assert int(hd[0].header['OBSID']) == obsid
                assert 'COADD' in hd
            tmp.replace(path)
        finally:
            if tmp.exists():
                tmp.unlink()
    return path, url


def extract_spectrum(path, obsid, output):  # Preserve measured arrays and produce a usable-bin mask.
    with fits.open(path) as hd:
        # Fix known empty-primary NAXIS metadata in memory; do not modify cached original.
        hd.verify('silentfix+exception')
        header = hd[0].header
        if int(header['OBSID']) != obsid:
            raise ValueError('Spectrum observation ID does not match the requested ID')
        data = hd['COADD'].data[0]
        wavelength = np.array(data['WAVELENGTH'], dtype=float)
        flux = np.array(data['FLUX'], dtype=float)
        ivar = np.array(data['IVAR'], dtype=float)
        andmask = np.array(data['ANDMASK'], dtype=np.int64)
        ormask = np.array(data['ORMASK'], dtype=np.int64)
        normalization = np.array(data['NORMALIZATION'], dtype=float)
        fiber_mask = int(header['FIB_MASK'])
    assert len(wavelength) == len(flux) == len(ivar)
    assert np.isfinite(wavelength).all() and np.all(np.diff(wavelength) > 0)
    good = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0) & (ormask == 0)
    if not good.any():
        raise ValueError('No usable spectral bins under the applied mask')
    sigma = np.full_like(ivar, np.nan)
    valid_ivar = np.isfinite(ivar) & (ivar > 0)
    sigma[valid_ivar] = 1 / np.sqrt(ivar[valid_ivar])
    target = output / f'spectrum_{obsid}.npz'
    np.savez_compressed(target, wavelength_angstrom=wavelength, catalog_flux=flux,
        catalog_flux_error=sigma, inverse_variance=ivar, andmask=andmask,
        ormask=ormask, catalog_normalization=normalization, usable=good)
    return {'arrays_file': str(target.resolve()), 'original_fits': str(path.resolve()),
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'sample_count':len(flux),
        'usable_sample_count':int(good.sum()), 'wavelength_min_angstrom':float(wavelength[0]),
        'wavelength_max_angstrom':float(wavelength[-1]), 'spectrum_fibermask':fiber_mask,
        'flux_description':'Original LAMOST catalog flux scale; no absolute calibration assumed.'}


def retrieve(matches_path, catalog_path, cache, output):  # Process every supplied association by exact obsid.
    output.mkdir(parents=True, exist_ok=True)
    inputs, selected = load_matches(matches_path)
    if not selected:
        raise ValueError('No existing LAMOST matches found; matching must happen before retrieval')
    catalog = parse_single_table(catalog_path).to_table(use_names_over_ids=True)
    properties, failures = [], []
    fields = {
        'designation':'designation', 'observation_date':'obsdate',
        'gaia_dr3_source_id':'gaia_source_id', 'temperature_K':'teff',
        'temperature_catalog_error_K':'teff_err', 'Fe_H_dex':'feh',
        'Fe_H_catalog_error_dex':'feh_err', 'logg_dex':'logg',
        'logg_catalog_error_dex':'logg_err', 'radial_velocity_km_s':'rv',
        'radial_velocity_catalog_error_km_s':'rv_err', 'spectral_class':'class',
        'spectral_subclass':'subclass', 'snrg':'snrg','snrr':'snrr',
        'fibermask':'fibermask', 'catalog_ra_deg':'ra', 'catalog_dec_deg':'dec'}
    for original, source_id, obsid in selected:
        try:
            indices = np.flatnonzero(np.asarray(catalog['obsid']) == obsid)
            if len(indices) != 1:
                raise ValueError(f'Expected exactly one catalog row for obsid {obsid}, got {len(indices)}; supply the correct --catalog')
            row = catalog[int(indices[0])]
            star = {'source_id':source_id, 'lamost_obsid':obsid,
                    'release':'DR11 v2.0 LRS', 'original_match_record':original}
            for outkey, col in fields.items():
                star[outkey] = scalar(row[col]) if col in catalog.colnames else None
            if star['gaia_dr3_source_id'] is not None:
                star['gaia_dr3_source_id'] = str(star['gaia_dr3_source_id'])
            path, url = fetch_spectrum(obsid, cache)
            star['spectrum'] = extract_spectrum(path, obsid, output)
            star['spectrum']['download_url'] = url
            if star['fibermask'] != star['spectrum']['spectrum_fibermask']:
                raise ValueError('Catalog and spectrum fiber masks disagree')
            star['quality'] = {'fiber_mask_zero':star['fibermask']==0,
                'snrg_above_20':star['snrg'] is not None and star['snrg']>20,
                'stellar_parameters_available':all(star[k] is not None for k in ['temperature_K','Fe_H_dex','logg_dex']),
                'note':'Preserved supplied association; no new position match or scientific quality certification.'}
            properties.append(star)
        except Exception as exc:
            failures.append({'source_id':source_id,'lamost_obsid':obsid,'error':str(exc)})
    write_json(output/'stellar_properties.json', properties)
    manifest = {'run_utc':datetime.now(timezone.utc).isoformat(),
        'input_file':str(matches_path.resolve()),'catalog_file':str(catalog_path.resolve()),
        'catalog_sha256':hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
        'input_records':len(inputs),'unique_associations_requested':len(selected),
        'records_without_lamost_obsid':sum(r.get('lamost_obsid') is None for r in inputs),
        'retrieved':len(properties),'failed':failures,
        'status':'PASS' if not failures else 'PARTIAL_FAILURE',
        'notes':['Parameters are catalog estimates, not newly fitted from spectra.',
                 'Catalog parameter errors omit some systematic uncertainty.',
                 'Spectrum flux-error array is derived as 1/sqrt(inverse variance).',
                 'Repeated observations remain distinct; no averaging or coaddition performed.']}
    write_json(output/'retrieval_manifest.json',manifest)
    return manifest


def main():  # Parse paths and run retrieval; default inputs are the existing M67 matches.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matches',type=Path,default=ROOT/'sources.json')
    parser.add_argument('--catalog',type=Path,default=ROOT/'data/lamost_cone.xml')
    parser.add_argument('--cache',type=Path,default=ROOT/'data')
    parser.add_argument('--output',type=Path,default=ROOT/'retrieved')
    args = parser.parse_args()
    result = retrieve(args.matches,args.catalog,args.cache,args.output)
    print(json.dumps(result,indent=2))
    if result['failed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
