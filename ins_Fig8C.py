import FSIN as sim
import numpy as np
from scipy.signal import cwt, morlet2
import matplotlib.pyplot as plt
from funcs import *

method = "rk4"; # methods "rk4" and "gsl_rk8pd" give similar results
dt = 1.e-3;

conn_mat = 0;
sim_time = 5.*125.;
fsin = 8
gsin=7.

rates = sim.gewnet(gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=True,gjs=True,method=method,dt=dt,connectivity_matrix=conn_mat,Es=-75.)[0][1];

np.savetxt("RatesFig8C1.dat",rates)
rates = np.loadtxt("RatesFig8C1.dat")

dt_rec = 1e-1; # Down sample the rates

fs = np.arange(50.,450.1,3.);
Ntcyc = int(np.round(125./dt));
cwtPow = compPWT(rates[-Ntcyc::int(np.round(dt_rec/dt))],dt=dt_rec*1e-3,fs=fs);

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

tmp = plot1(cwtPow,Ncycs=1,dt=dt_rec,fs=fs,nlevels=30);

tmp = plt.savefig("Figures/Fig8C1Pow_%dHzgChR%d.eps" % (int(fsin),int(gsin)), dpi=300); # To save plot in eps file
tmp = plt.savefig("Figures/Fig8C1Pow_%dHzgChR%d.png" %  (int(fsin),int(gsin)), dpi=300); tmp = plt.clf(); #


cwtPow = compPhaseWT(rates[-Ntcyc::int(np.round(dt_rec/dt))],dt=dt_rec*1e-3,fs=fs);

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

tmp = plotPhase(cwtPow,Ncycs=1,dt=dt_rec,fs=fs,nlevels=30);

tmp = plt.savefig("Figures/Fig8C1Phase_%dHzgChR%d.eps" % (int(fsin),int(gsin)), dpi=300); #tmp = plt.clf(); # To save plot in eps file
tmp = plt.savefig("Figures/Fig8C1Phase_%dHzgChR%d.png" %  (int(fsin),int(gsin)), dpi=300); tmp = plt.clf(); #
