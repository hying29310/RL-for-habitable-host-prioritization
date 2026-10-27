"""Select 70 NEW matched stars with complete catalog T, [Fe/H], radius and luminosity.
Selection uses spectrum quality, not the subsequent host-priority score.
"""
import json
import numpy as np
from astropy.io.votable import parse_single_table
from validate_matches import ROOT,save

def main():
    matches=json.loads((ROOT/'matches.json').read_text())
    g=json.loads((ROOT/'gaia_sources.json').read_text());old={x['gaia_dr3_source_id'] for x in json.loads((ROOT/'excluded_existing_stars.json').read_text())}
    t=parse_single_table(ROOT/'selected_lamost.xml').to_table(use_names_over_ids=True);byobs={int(r['obsid']):r for r in t}
    complete=[];incomplete=[]
    for m in matches:
        gid=m['gaia_dr3_source_id'];assert gid not in old
        row=byobs[m['lamost_obsid']];ap=g[gid]
        keys=['radius_flame','lum_flame','radius_flame_lower','radius_flame_upper','lum_flame_lower','lum_flame_upper']
        if not all(ap.get(k) is not None and np.isfinite(ap[k]) and ap[k]>0 for k in keys):
            incomplete.append({'gaia_dr3_source_id':gid,'field':m['field'],'reason':'missing_or_invalid_FLAME_radius_luminosity_or_bounds'});continue
        m['gaia_astrophysical_parameters']={k:ap[k] for k in keys+['flags_flame','teff_gspphot']}
        m['snrg']=float(row['snrg']);complete.append(m)
    complete.sort(key=lambda m:(m['fibermask']!=0,-m['snrg'],m['gaia_dr3_source_id']))
    save(ROOT/'complete_matched_pool.json',complete);save(ROOT/'incomplete_parameters.json',incomplete)
    print('New matched:',len(matches),'complete:',len(complete),'zero fiber mask:',sum(m['fibermask']==0 for m in complete),flush=True)
    if len(complete)<70:raise SystemExit('Continue searching: fewer than 70 complete stars.')
    chosen=complete[:70];save(ROOT/'selected_70_matches.json',chosen)
    assert len({m['gaia_dr3_source_id'] for m in chosen})==70
    print('Selected 70 by zero fiber mask first, then SNR(g), then Gaia ID.',flush=True)
if __name__=='__main__':main()
