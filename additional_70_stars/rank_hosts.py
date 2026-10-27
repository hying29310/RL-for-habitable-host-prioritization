"""Transparent, UNCALIBRATED search-priority score; never a planet probability.
Only T, [Fe/H], the measured spectrum/its quality, L and R enter this score.
"""
from pathlib import Path
import json,math,itertools,hashlib,re
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parent
CONFIG={'temperature_plateau_K':[4000.,6000.],'temperature_tail_width_K':1000.,
        'dwarf_deviation_width_dex':0.30,'metallicity_max_penalty':0.10,
        'metallicity_midpoint_dex':-0.5,'metallicity_width_dex':0.3,
        'spectrum_snr_saturation':50.,'nonzero_fiber_mask_factor':0.5,
        'reference':'Pecaut & Mamajek 2013 / Mamajek updated dwarf sequence, version 2024.05.15',
        'not_probability':True}

def sequence():
    """Read author-supplied dwarf T, radius and luminosity sequence; no extrapolation."""
    rows=[]
    for line in (ROOT/'data/dwarf_sequence.txt').read_text().splitlines():
        a=line.split()
        if len(a)<7 or not re.match(r'^[OBAFGKMLTY][0-9].*V$',a[0]):continue
        try:temp,logl,r=float(a[1]),float(a[4]),float(a[6])
        except ValueError:continue
        if temp>0 and r>0:rows.append((temp,np.log10(r),logl))
    rows=sorted(set(rows));assert len(rows)>30
    return np.array(rows)
SEQ=sequence()

def components(T,Z,R,L,snr,usable_fraction,fibermask,metal_penalty=None,dwarf_width=None):
    """Compute all score factors, preserving physical preference versus data quality."""
    if not np.isfinite([T,Z,R,L,snr,usable_fraction]).all() or R<=0 or L<=0:return None
    if not SEQ[0,0]<=T<=SEQ[-1,0]:return None
    lo,hi=CONFIG['temperature_plateau_K'];dt=max(lo-T,0,T-hi)
    ft=math.exp(-0.5*(dt/CONFIG['temperature_tail_width_K'])**2)
    log_rms=float(np.interp(T,SEQ[:,0],SEQ[:,1]));log_lms=float(np.interp(T,SEQ[:,0],SEQ[:,2]))
    dr=abs(math.log10(R)-log_rms);dl=0.5*abs(math.log10(L)-log_lms)
    # Max rather than product: luminosity and radius are strongly correlated.
    deviation=max(dr,dl);width=CONFIG['dwarf_deviation_width_dex'] if dwarf_width is None else dwarf_width
    fd=math.exp(-0.5*(deviation/width)**2)
    penalty=CONFIG['metallicity_max_penalty'] if metal_penalty is None else metal_penalty
    # Deliberately weak, adjustable hypothesis: not an empirical occurrence law.
    logistic=1/(1+math.exp(-(Z-CONFIG['metallicity_midpoint_dex'])/CONFIG['metallicity_width_dex']))
    fz=1-penalty+penalty*logistic
    qs=float(np.clip(usable_fraction,0,1))*float(np.clip(snr/CONFIG['spectrum_snr_saturation'],0,1))
    if fibermask!=0:qs*=CONFIG['nonzero_fiber_mask_factor']
    physical=100*ft*fd*fz
    return {'priority_score':physical*qs,'physical_preference_score':physical,'temperature_factor':ft,
        'dwarf_factor':fd,'metallicity_factor':fz,'spectrum_quality_factor':qs,
        'reference_radius_solar':10**log_rms,'reference_luminosity_solar':10**log_lms,
        'dwarf_deviation_dex':deviation,
        'stefan_boltzmann_log10_residual':math.log10(L/(R*R*(T/5772.)**4))}

def main():
    ms=json.loads((ROOT/'selected_70_matches.json').read_text());props=json.loads((ROOT/'retrieved_70/stellar_properties.json').read_text());ps={p['gaia_dr3_source_id']:p for p in props}
    cuts={r['gaia_dr3_source_id']:r for r in json.loads((ROOT/'cutout_manifest_70.json').read_text())};out=[]
    assert len(ms)==len(ps)==70
    for m in ms:
        gid=m['gaia_dr3_source_id'];p=ps[gid];a=m['gaia_astrophysical_parameters'];spec=p['spectrum']
        T=p['temperature_K'];Z=p['Fe_H_dex'];R=a['radius_flame'];L=a['lum_flame'];frac=spec['usable_sample_count']/spec['sample_count'];snr=p['snrg'];mask=p['fibermask']
        c=components(T,Z,R,L,snr,frac,mask);assert c is not None
        terr=p['temperature_catalog_error_K'];zerr=p['Fe_H_catalog_error_dex']
        # Sensitivity envelope of marginal intervals, NOT a joint confidence interval.
        grid=[]
        for tt,zz,rr,ll in itertools.product([T-terr,T+terr],[Z-zerr,Z+zerr],[a['radius_flame_lower'],a['radius_flame_upper']],[a['lum_flame_lower'],a['lum_flame_upper']]):
            cc=components(tt,zz,rr,ll,snr,frac,mask)
            if cc is not None:grid.append(cc['priority_score'])
        no_z=components(T,Z,R,L,snr,frac,mask,metal_penalty=0)['priority_score']
        wide=components(T,Z,R,L,snr,frac,mask,dwarf_width=.45)['priority_score']
        warnings=[]
        if mask:warnings.append('LAMOST fiber quality flags present')
        if abs(c['stefan_boltzmann_log10_residual'])>.2:warnings.append('T, L and R disagree with Stefan-Boltzmann relation by >0.2 dex; review parameter estimates')
        if not 2600<=T<=7200:warnings.append('Outside usual Kopparapu 2014 HZ-temperature range; score is not an HZ calculation')
        if m['gaia_ruwe'] is not None and m['gaia_ruwe']>1.4:warnings.append('Elevated Gaia RUWE; not used in score')
        out.append({'gaia_dr3_source_id':gid,'field':m['field'],'lamost_obsid':m['lamost_obsid'],
            'temperature_K':T,'temperature_error_K':terr,'Fe_H_dex':Z,'Fe_H_error_dex':zerr,
            'radius_solar':R,'radius_16th_percentile_solar':a['radius_flame_lower'],'radius_84th_percentile_solar':a['radius_flame_upper'],
            'luminosity_solar':L,'luminosity_16th_percentile_solar':a['lum_flame_lower'],'luminosity_84th_percentile_solar':a['lum_flame_upper'],
            'radius_luminosity_source':'Gaia DR3 FLAME','gaia_flags_flame':a['flags_flame'],
            'spectral_subclass':p['spectral_subclass'],'spectrum':spec,'spectrum_snrg':snr,'fiber_mask':mask,
            'usable_spectrum_fraction':frac,'score_components':c,
            'sensitivity_envelope_low':min(grid+[c['priority_score']]),'sensitivity_envelope_high':max(grid+[c['priority_score']]),
            'score_without_metallicity_preference':no_z,'score_with_wider_dwarf_tolerance':wide,
            'warnings':warnings,'match_evidence':m,'hubble_cutout':cuts.get(gid),
            'interpretation':'Uncalibrated search priority for this chosen stellar target class, not probability of planets, habitability or life.'})
    out.sort(key=lambda r:(-r['score_components']['priority_score'],r['gaia_dr3_source_id']))
    for i,r in enumerate(out,1):r['rank']=i
    primary=[r['score_components']['priority_score'] for r in out];alt=[r['score_without_metallicity_preference'] for r in out]
    top10={r['gaia_dr3_source_id'] for r in out[:10]};alt10={r['gaia_dr3_source_id'] for r in sorted(out,key=lambda r:-r['score_without_metallicity_preference'])[:10]}
    sensitivity={'metallicity_off_spearman':float(spearmanr(primary,alt).statistic),'metallicity_off_top10_overlap':len(top10&alt10),
                 'envelope_warning':'Parameter-box sensitivity only. Marginal Gaia R/L bounds are correlated; this is not a confidence interval or probability.'}
    for name,obj in [('ranking_config.json',CONFIG),('ranked_70_stars.json',out),('ranking_sensitivity.json',sensitivity)]:
        (ROOT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    print('Ranked',len(out),'new stars. Top score',out[0]['score_components']['priority_score']);print(sensitivity)
if __name__=='__main__':main()
