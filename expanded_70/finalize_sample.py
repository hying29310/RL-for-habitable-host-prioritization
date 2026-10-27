"""Select 70 complete new matches by data quality, retrieve spectra, assemble parameters."""
from pathlib import Path
import json
from concurrent.futures import ThreadPoolExecutor,as_completed
from astropy.io.votable import parse_single_table
from retrieve_matched_stars import retrieve,fetch_spectrum
ROOT=Path(__file__).resolve().parent

def main():
    ms=json.loads((ROOT/'matches.json').read_text())
    ap=json.loads((ROOT/'gaia_stellar_parameters.json').read_text())
    cat=parse_single_table(ROOT/'selected_lamost.xml').to_table(use_names_over_ids=True)
    snr={int(r['obsid']):float(r['snrg']) for r in cat}
    eligible=[m for m in ms if m['fibermask']==0 and all(ap.get(m['gaia_dr3_source_id'],{}).get(k) is not None and ap[m['gaia_dr3_source_id']][k]>0 for k in ['radius_flame','lum_flame'])]
    eligible.sort(key=lambda m:(not m['hst_has_zero_flag_detection'],-snr[m['lamost_obsid']],m['gaia_dr3_source_id']))
    print('Eligible DR3 complete:',len(eligible),flush=True)
    if len(eligible)<70:raise RuntimeError('Need 70 complete DR3 stars; continue search')
    selected=eligible[:70]
    (ROOT/'selected_70_matches.json').write_text(json.dumps(selected,indent=2)+'\n')
    with ThreadPoolExecutor(4) as pool:
        fs=[pool.submit(fetch_spectrum,m['lamost_obsid'],ROOT/'data/spectra') for m in selected]
        for i,f in enumerate(as_completed(fs),1):
            f.result()
            if i%10==0:print('Spectra cached',i,flush=True)
    result=retrieve(ROOT/'selected_70_matches.json',ROOT/'selected_lamost.xml',ROOT/'data/spectra',ROOT/'retrieved')
    if result['failed']:raise RuntimeError(result)
    props=json.loads((ROOT/'retrieved/stellar_properties.json').read_text());out=[]
    for s in props:
        m=s['original_match_record'];a=ap[m['gaia_dr3_source_id']]
        fields=['gaia_dr3_source_id','lamost_obsid','designation','observation_date','temperature_K','temperature_catalog_error_K','Fe_H_dex','Fe_H_catalog_error_dex','spectral_class','spectral_subclass','spectrum']
        row={k:s[k] for k in fields}
        row.update(field=m['field'],hsc_match_id=m['hsc_match_id'],ra_deg=m['ra_deg'],dec_deg=m['dec_deg'],radius_solar=a['radius_flame'],radius_lower_solar=a['radius_flame_lower'],radius_upper_solar=a['radius_flame_upper'],luminosity_solar=a['lum_flame'],luminosity_lower_solar=a['lum_flame_lower'],luminosity_upper_solar=a['lum_flame_upper'],temperature_metallicity_source='LAMOST DR11 v2.0 LRS',radius_luminosity_source='Gaia DR3 FLAME',quality={'lamost_snrg':s['snrg'],'lamost_fibermask':s['fibermask'],'gaia_flags_flame':a['flags_flame'],'hst_has_zero_flag_detection':m['hst_has_zero_flag_detection'],'hst_detection_flags':m['hst_detection_flags'],'median_hst_match_separation_arcsec':m['median_epoch_corrected_separation_arcsec']})
        out.append(row)
    (ROOT/'complete_star_parameters.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    (ROOT/'selection_method.json').write_text(json.dumps({'accepted_new_matches':len(ms),'complete_dr3_matches':len(eligible),'selected':len(out),'rule':'Complete Gaia DR3 FLAME R/L and zero LAMOST fiber mask; prioritize at least one zero-flag HST detection, then descending LAMOST g SNR, then Gaia ID. No host score used for selection.'},indent=2)+'\n')
    keep={Path(s['spectrum']['arrays_file']).name for s in out}
    for path in (ROOT/'retrieved').glob('spectrum_*.npz'):
        if path.name not in keep:path.unlink()
    print('Complete stars',len(out),flush=True)
if __name__=='__main__':main()
