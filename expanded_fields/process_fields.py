"""Targeted source verification in new HST fields using existing LAMOST stars.
Unlike an all-sky survey, selected fields and archive products are explicit.
Input acquisition: acquire_fields.py. Cached reruns need no network except missing spectra.
"""
from pathlib import Path
import os,json,gc,warnings,hashlib
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mplconfig'))
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord,Distance
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
import astropy.units as u
from photutils.detection import DAOStarFinder
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from retrieve_matched_stars import retrieve

IMAGES={'NGC6791':'ibc603zoq_drc.fits','NGC6819':'ib2o06020_drc.fits'}

def write(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def at_epoch(g,epoch):
    kw=dict(ra=g['ra']*u.deg,dec=g['dec']*u.deg,
        pm_ra_cosdec=g['pmra']*u.mas/u.yr,pm_dec=g['pmdec']*u.mas/u.yr,
        obstime=Time(g['ref_epoch'],format='jyear'))
    if g['parallax'] is not None and g['parallax']>0:
        kw.update(distance=Distance(parallax=g['parallax']*u.mas),radial_velocity=0*u.km/u.s)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='.*distance overridden.*')
        return SkyCoord(**kw).apply_space_motion(new_obstime=epoch)

def process(name,filename):
    p=ROOT/'data'/name;out=ROOT/name;out.mkdir(exist_ok=True)
    t=parse_single_table(p/'lamost.xml').to_table(use_names_over_ids=True)
    gj=json.loads((p/'gaia.json').read_text());cols=[c['name'] for c in gj['metadata']]
    gaia={str(row[0]):dict(zip(cols,row)) for row in gj['data']}
    eligible={}
    for row in t:
        keys=['gaia_source_id','teff','feh','snrg']
        if any(np.ma.is_masked(row[k]) for k in keys):continue
        if not np.isfinite([row['teff'],row['feh'],row['snrg']]).all():continue
        if row['teff']<=0 or row['feh']<=-5 or row['snrg']<=20 or str(row['class'])!='STAR':continue
        eligible.setdefault(str(row['gaia_source_id']),[]).append(row)
    for rows in eligible.values():rows.sort(key=lambda r:(int(r['fibermask'])!=0,-float(r['snrg'])))
    matches=[];audit=[];thumbs=[]
    with fits.open(p/filename) as h:
        ph=h[0].header;hdr={k:ph.get(k) for k in ['TELESCOP','INSTRUME','DETECTOR','TARGNAME','FILTER','DATE-OBS','EXPTIME','EXPSTART']}
        assert hdr['TELESCOP']=='HST'
        data=h['SCI'].data;weight=h['WHT'].data;w=WCS(h['SCI'].header,fobj=h).celestial
        epoch=Time(ph['EXPSTART'],format='mjd');overview=np.array(data[::4,::4],float)
        for gid,rows in eligible.items():
            a={'gaia_dr3_source_id':gid,'available_parameter_observations':len(rows)}
            if gid not in gaia or any(gaia[gid].get(k) is None for k in ['pmra','pmdec']):
                a['status']='missing_Gaia_astrometry';audit.append(a);continue
            g=gaia[gid];pred=at_epoch(g,epoch);px,py=w.all_world2pix(pred.ra.deg,pred.dec.deg,0);px=float(px);py=float(py)
            a.update(predicted_x=px,predicted_y=py)
            if not(10<px<data.shape[1]-10 and 10<py<data.shape[0]-10):
                a['status']='outside_image';audit.append(a);continue
            x0=max(0,int(px)-64);x1=min(data.shape[1],int(px)+65)
            y0=max(0,int(py)-64);y1=min(data.shape[0],int(py)+65)
            cut=np.array(data[y0:y1,x0:x1],float);wt=np.array(weight[y0:y1,x0:x1],float)
            mask=~np.isfinite(cut)|~np.isfinite(wt)|(wt<=0)
            if mask[int(py)-y0,int(px)-x0]:
                a['status']='no_image_coverage';audit.append(a);continue
            _,med,std=sigma_clipped_stats(cut,mask=mask,sigma=3)
            if not np.isfinite(std) or std<=0:
                a['status']='invalid_noise';audit.append(a);continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                found=DAOStarFinder(threshold=5*std,fwhm=2.5,exclude_border=True)(cut-med,mask=mask)
            if found is None:
                a['status']='no_detected_source';audit.append(a);continue
            xx=np.asarray(found['xcentroid'])+x0;yy=np.asarray(found['ycentroid'])+y0
            ra,dec=w.all_pix2world(xx,yy,0);sky=SkyCoord(ra*u.deg,dec*u.deg);sep=sky.separation(pred).arcsec;order=np.argsort(sep);i=int(order[0])
            near=int(np.sum(sep<0.3));a.update(nearest_separation_arcsec=float(sep[i]),detections_within_0p3_arcsec=near,
                second_nearest_arcsec=float(sep[order[1]]) if len(order)>1 else None)
            if sep[i]>=0.3 or near!=1:
                a['status']='rejected_offset_or_ambiguous';audit.append(a);continue
            row=rows[0];neighbors=sky.separation(sky[i]).arcsec
            neighbor_count=int(np.sum((neighbors>0.001)&(neighbors<1.65)))
            obsid=int(row['obsid']);source=f'{name}_Gaia_{gid}'
            m={'source_id':source,'field':name,'gaia_dr3_source_id':gid,'lamost_obsid':obsid,
               'lamost_release':'DR11 v2.0 LRS','match_status':'positionally_validated',
               'hst_product':filename,'hst_date':hdr['DATE-OBS'],'ra_deg':float(ra[i]),'dec_deg':float(dec[i]),
               'hst_x':float(xx[i]),'hst_y':float(yy[i]),'epoch_corrected_match_arcsec':float(sep[i]),
               'second_nearest_arcsec':a['second_nearest_arcsec'],'other_detections_within_1p65_arcsec':neighbor_count,
               'fibermask':int(row['fibermask']),'gaia_ruwe':g['ruwe'],
               'caution':'Crowding/artifacts and HST saturation not fully audited; cluster membership not established.'}
            matches.append(m);a['status']='accepted';a['selected_lamost_obsid']=obsid;audit.append(a)
            thumbs.append((cut,xx[i]-x0,yy[i]-y0,px-x0,py-y0,source))
    write(out/'matches.json',matches);write(out/'audit.json',audit)
    retrieval=retrieve(out/'matches.json',p/'lamost.xml',p,out/'retrieved') if matches else None
    fig,ax=plt.subplots(figsize=(8,8),layout='constrained');ax.imshow(overview,origin='lower',cmap='gray',norm=simple_norm(overview,'asinh',percent=99.5))
    for i,m in enumerate(matches):
        ax.scatter(m['hst_x']/4,m['hst_y']/4,s=100,facecolors='none',edgecolors='cyan');ax.annotate(str(i+1),(m['hst_x']/4+8,m['hst_y']/4+8),color='cyan')
    ax.set(title=f"{name} field · {hdr['FILTER']} · {hdr['DATE-OBS']}\n{len(matches)} Hubble–LAMOST–Gaia positional matches",xlabel='Column (4× downsampled)',ylabel='Row (4× downsampled)')
    fig.savefig(out/'field.png',dpi=140);plt.close(fig)
    if thumbs:
        fig,axs=plt.subplots(1,len(thumbs),figsize=(4*len(thumbs),4),squeeze=False,layout='constrained')
        for i,(cut,x,y,px,py,source) in enumerate(thumbs):
            ax=axs[0,i];ax.imshow(cut,origin='lower',cmap='gray',norm=simple_norm(cut,'asinh',percent=99.5));ax.scatter(x,y,s=100,facecolors='none',edgecolors='cyan');ax.scatter(px,py,marker='+',s=80,color='orange');ax.set_xlim(x-25,x+25);ax.set_ylim(y-25,y+25);ax.set_title(f"Match {i+1}: {matches[i]['epoch_corrected_match_arcsec']:.3f} arcsec")
        fig.suptitle('Cyan: detection; orange: Gaia at Hubble epoch');fig.savefig(out/'matches.png',dpi=140);plt.close(fig)
    summary={'field':name,'observation':hdr,'image':filename,'eligible_unique_stars_in_cone':len(eligible),
        'positionally_matched_stars':len(matches),'retrieval':retrieval,
        'audit_status_counts':{status:sum(a['status']==status for a in audit) for status in sorted(set(a['status'] for a in audit))},
        'sha256_image':hashlib.sha256((p/filename).read_bytes()).hexdigest()}
    write(out/'summary.json',summary);print(json.dumps(summary),flush=True)
    gc.collect();return matches,summary

def main():
    all_matches=[];summaries=[]
    for name,filename in IMAGES.items():
        m,s=process(name,filename);all_matches+=m;summaries.append(s)
    write(ROOT/'all_new_matches.json',all_matches)
    props=[]
    for name in IMAGES:
        f=ROOT/name/'retrieved/stellar_properties.json'
        if f.exists():props+=json.loads(f.read_text())
    write(ROOT/'all_new_stellar_properties.json',props)
    write(ROOT/'expansion_summary.json',{'run_utc':datetime.now(timezone.utc).isoformat(),
        'new_unique_stars':len(set(m['gaia_dr3_source_id'] for m in all_matches)),
        'fields':summaries,'caveat':'Targeted verification of preselected LAMOST stars; not blind image completeness or cluster-membership analysis.'})

if __name__=='__main__':main()
