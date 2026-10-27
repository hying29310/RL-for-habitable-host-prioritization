"""Create a readable plain-text results table and verify the delivered sample."""
from pathlib import Path
from collections import Counter
import json,numpy as np
from rank_hosts import terms,dwarf_reference
ROOT=Path(__file__).resolve().parent

def main():
    stars=json.loads((ROOT/'ranked_stars.json').read_text())
    known={s['gaia_dr3_source_id'] for s in json.loads((ROOT/'previously_known_matches.json').read_text())}
    assert len(stars)==70
    assert len({s['gaia_dr3_source_id'] for s in stars})==70
    assert len({s['hsc_match_id'] for s in stars})==70
    assert not known.intersection(s['gaia_dr3_source_id'] for s in stars)
    for s in stars:
        assert all(np.isfinite(s[k]) for k in ['temperature_K','Fe_H_dex','radius_solar','luminosity_solar'])
        with np.load(s['spectrum']['arrays_file']) as a:
            assert len(a['wavelength_angstrom'])==s['spectrum']['sample_count']
            assert np.sum(a['usable'])==s['spectrum']['usable_sample_count']>0
        assert Path(s['spectrum']['original_fits']).exists()
        assert 0<=s['ranking']['score']<=100
    ref=dwarf_reference()
    assert terms(5772,0,1,1,ref)['score']>terms(5772,0,10,100,ref)['score']
    assert terms(5200,0,0,1,ref) is None
    assert terms(5200,-3,.8,.4,ref)['score']/terms(5200,.5,.8,.4,ref)['score']>.97
    summary={'new_unique_stars':70,'overlap_with_previous_34':0,'spectra_downloaded_and_checked':70,'radius_luminosity_source':'Gaia DR3 FLAME for all 70','field_counts':dict(Counter(s['field'] for s in stars)),'hst_without_zero_flag_detection':sum(not s['quality']['hst_has_zero_flag_detection'] for s in stars),'nonzero_fiber_mask':sum(s['quality']['lamost_fibermask']!=0 for s in stars),'large_T_R_L_disagreement_over_0p3_dex':sum(abs(s['ranking']['stefan_boltzmann_log10_residual'])>.3 for s in stars),'checks':'PASS'}
    (ROOT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['70 NEW HUBBLE--LAMOST--GAIA MATCHES','Score is a heuristic priority, NOT a percentage or a detected planet.','T: kelvin; Fe/H: dex; R and L: solar units. Field names indicate search pointing, NOT cluster/galaxy membership.','Detailed errors, quality flags, HST detections and spectra are in the JSON files.','',f"{'Rank':>4} {'Gaia DR3 source ID':>19} {'LAMOST obsid':>12} {'Field':<15} {'Type':<8} {'T(K)':>8} {'Fe/H':>7} {'R/Rsun':>9} {'L/Lsun':>10} {'Score':>7}"]
    for s in stars:
        lines.append(f"{s['rank']:4d} {s['gaia_dr3_source_id']:>19} {s['lamost_obsid']:12d} {s['field']:<15} {s['spectral_subclass']:<8} {s['temperature_K']:8.1f} {s['Fe_H_dex']:7.3f} {s['radius_solar']:9.3f} {s['luminosity_solar']:10.3f} {s['ranking']['score']:7.2f}")
    (ROOT/'RANKED_STARS.txt').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
