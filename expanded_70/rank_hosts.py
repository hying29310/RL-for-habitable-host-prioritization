"""Transparent, uncalibrated stellar search-priority baseline.
Only T_eff, [Fe/H], spectrum-derived information, luminosity and radius are used.
The score is NOT a host probability or evidence for a planet. Spectrum quality is
reported separately, not multiplied into an astrophysical host likelihood.
"""
from pathlib import Path
import json,math,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parent

def dwarf_reference():
    """Read the versioned Pecaut/Mamajek dwarf sequence; interpolate in log space."""
    rows=[]
    for line in (ROOT/'reference/mamajek_dwarf_sequence.txt').read_text().splitlines():
        if not line or line.startswith('#'):continue
        x=line.split()
        if len(x)<7 or not x[0].endswith('V'):continue
        try:t=float(x[1]);logl=float(x[4]);r=float(x[6])
        except ValueError:continue
        if t>0 and r>0 and 2000<=t<=15000:rows.append((t,r,10**logl))
    rows.sort();return np.asarray(rows,float)

def terms(t,z,r,l,ref,center=5200,width=1500,rscale=.30,lscale=.60):
    """Evaluate explicit design preferences, without pretending they are fitted laws."""
    if not all(np.isfinite([t,z,r,l])) or min(t,r,l)<=0:return None
    if not ref[0,0]<=t<=ref[-1,0]:return None
    rmain=10**np.interp(np.log10(t),np.log10(ref[:,0]),np.log10(ref[:,1]))
    lmain=10**np.interp(np.log10(t),np.log10(ref[:,0]),np.log10(ref[:,2]))
    dt=(t-center)/width
    dr=math.log10(r/rmain)/rscale;dl=math.log10(l/lmain)/lscale
    d=max(abs(dr),abs(dl))
    log_st=-.5*dt*dt;log_sms=-.5*d*d
    sz=.9+.1*math.tanh((z+.5)/.5)
    logscore=.4*log_st+.5*log_sms+.1*math.log(sz)
    return {'score':100*math.exp(max(-745,logscore)),
            'temperature_preference':math.exp(max(-745,log_st)),
            'main_sequence_similarity':math.exp(max(-745,log_sms)),
            'weak_metallicity_factor':sz,'reference_radius_solar':rmain,
            'reference_luminosity_solar':lmain,'main_sequence_distance':d,
            'log_temperature_term':log_st,'log_main_sequence_term':log_sms,
            'stefan_boltzmann_log10_residual':math.log10(l/(r*r*(t/5772)**4))}

def evaluate(star,ref):
    """Rank complete stars and retain uncertainties/flags without fabricating missing values."""
    t=star.get('temperature_K');z=star.get('Fe_H_dex');r=star.get('radius_solar');l=star.get('luminosity_solar')
    if any(x is None for x in [t,z,r,l]):return {'rank_status':'missing_required_parameter','score':None}
    v=terms(t,z,r,l,ref)
    if v is None:return {'rank_status':'outside_reference_domain','score':None}
    # Sensitivity bounds, NOT posterior confidence intervals. Include zero Fe/H weight.
    variants=[]
    for center in [4800,5200,5600]:
        for scale in [.2,.3,.4]:
            a=terms(t,z,r,l,ref,center=center,rscale=scale,lscale=2*scale)
            variants.append(a['score'])
            variants.append(100*math.exp(max(-745,(.4*a['log_temperature_term']+.5*a['log_main_sequence_term'])/.9)))
    v.update(rank_status='heuristic_not_probability',sensitivity_min=min(variants),sensitivity_max=max(variants))
    return v

def main():
    records=json.loads((ROOT/'complete_star_parameters.json').read_text());ref=dwarf_reference()
    out=[]
    for s in records:out.append({**s,'ranking':evaluate(s,ref)})
    out.sort(key=lambda s:-(s['ranking']['score'] if s['ranking']['score'] is not None else -1))
    for i,s in enumerate(out,1):s['rank']=i if s['ranking']['score'] is not None else None
    (ROOT/'ranked_stars.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('Ranked',sum(s['rank'] is not None for s in out),'of',len(out),'stars')
if __name__=='__main__':main()
