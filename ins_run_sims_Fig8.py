import numpy as np
import FSIN as sim
import sys

method = "rk4"; # methods "rk4" and "gsl_rk8pd" give similar results
dt = 1.e-3;

#Ess = [-75.,-55.];
Ess = [-75.];

Nnets = 30;
Ncycs = 30;

# Total simulated time. Theta frequency is 8 Hz, i.e. period of 125 ms.
# sim_time given as multiple of period, i.e. number of simulated cycles.
# Only the last two are plotted. To remove effects of transients.
sim_time = (4.+Ncycs)*125.; 

# We apply a downsampling to the computed smoothed population rates.
# We save values at every dt_dt ms, i.e. the sampling rate is 1.0/dt_ds kHz.
dt_ds = .1;
kds = int(np.round(dt_ds/dt));
fsin=8.
gsin=7.

for Es in Ess:
   #for conn_mat in range(Nnets):
    conn_mat = int(sys.argv[1])
    ratesYSTD = sim.gewnet(gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=True,gjs=True,method=method,dt=dt,connectivity_matrix=conn_mat,Es=Es)[0][1];
    ratesNSTD = sim.gewnet(gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=True,gjs=True,method=method,dt=dt,connectivity_matrix=conn_mat,Es=Es,std=False)[0][1];
    ratesYSTD_ds = ratesYSTD[::kds];
    ratesNSTD_ds = ratesNSTD[::kds];
    tmp = np.savetxt("RatesFig8/rates_Esm%d_cm%d_%dHzgChR%d_YSTD.dat" % (int(np.round(-Es)),conn_mat,int(fsin),int(gsin)), ratesYSTD_ds );
    tmp = np.savetxt("RatesFig8/rates_Esm%d_cm%d_%dHzgChR%d_NSTD.dat" % (int(np.round(-Es)),conn_mat,int(fsin),int(gsin)), ratesNSTD_ds );
