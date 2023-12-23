import numpy as np
import FSIN as sim
import os
import brian2 as b2
import sys
import matplotlib.pyplot as plt
TigerFish = True
if TigerFish: 	# For Tigerfish use only
	cache_dir = os.path.expanduser('~/.cython/brian-pid-{}'.format(os.getpid()))
	b2.prefs.codegen.runtime.cython.cache_dir = cache_dir
	b2.prefs.codegen.runtime.cython.multiprocess_safe = False

method = "rk4"; # methods "rk4" and "gsl_rk8pd" give similar results
dt = 1.e-3;

Ess = [-55.];#-75.

gjs=True
read_ggaps=False
mod_gL=False
ggaps= 2.

for factor in [1,2,3,5]:
    gm_scale = factor;

    Nnets = 30;
    Ncycs = 30;

    # Total simulated time. Theta frequency is 8 Hz, i.e. period of 125 ms.
    # sim_time given as multiple of period, i.e. number of simulated cycles.
    # Only the last two are plotted. To remove effects of transients.
    sim_time = (4.+Ncycs)*250.; 

    # We apply a downsampling to the computed smoothed population rates.
    # We save values at every dt_dt ms, i.e. the sampling rate is 1.0/dt_ds kHz.
    dt_ds = .1;
    kds = int(np.round(dt_ds/dt));
    fsin=8.
    gsin=7.

    for Es in Ess:
        conn_mat = int(sys.argv[1]) 
       #for conn_mat in range(Nnets):
        rates = sim.gewnet(read_gms=False,scale_gm=True,gm_scale = gm_scale, ggaps=ggaps,gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=mod_gL,read_ggaps=read_ggaps,gjs=gjs,method=method,dt=dt,connectivity_matrix=conn_mat,Es=Es)[0][1];
        rates_ds = rates[::kds];
        tmp = np.savetxt("RatesFig7/NonPhysrates_Esm%d_cm%d_%dHzgChR%d_Fact%d.dat" % (int(np.round(-Es)),conn_mat,int(fsin),int(gsin), int(factor)), rates_ds );
