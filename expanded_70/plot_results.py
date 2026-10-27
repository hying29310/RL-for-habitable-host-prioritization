"""Plot the actual 70-star catalog and the spectrum of its highest-ranked target."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rank_hosts import dwarf_reference
ROOT=Path(__file__).resolve().parent
stars=json.loads((ROOT/'ranked_stars.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
a=axes[0];ref=dwarf_reference();a.plot(ref[:,0],ref[:,2],c='gray',lw=1,label='Dwarf reference')
c=a.scatter([s['temperature_K'] for s in stars],[s['luminosity_solar'] for s in stars],c=[s['ranking']['score'] for s in stars],cmap='viridis',vmin=0,vmax=100,s=35,edgecolors='black',linewidths=.3)
a.set(yscale='log',xlim=(max(s['temperature_K'] for s in stars)+300,3500),ylim=(.04,1000),xlabel='LAMOST temperature (K)',ylabel='Gaia luminosity (solar units)',title='70 new matched stars');a.legend(fontsize=8)
fig.colorbar(c,ax=a,label='Heuristic priority (not probability)')
s=stars[0];a=axes[1]
with np.load(s['spectrum']['arrays_file']) as x:
 good=x['usable'];a.plot(x['wavelength_angstrom'][good],x['catalog_flux'][good],lw=.45)
a.set(xlabel='Wavelength (angstrom)',ylabel='Catalog flux (not absolutely calibrated)',title=f"Highest priority: LAMOST obsid {s['lamost_obsid']}")
fig.savefig(ROOT/'sample_overview.png',dpi=160)
