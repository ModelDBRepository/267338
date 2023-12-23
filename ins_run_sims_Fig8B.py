# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import cwt, morlet2
import matplotlib.pyplot as plt
from funcs import *

dt = .1;
thresh = .65;

Ncycs = 30;
Nnets = 15;

fsin = 8
gsin=7.

Ess = [-75.];

Ntcyc = int(np.round(125./dt)); # Number of points per theta cycle
phases = np.linspace(-np.pi,+np.pi,Ntcyc); # Values of phases within ones cycle. -PI and +PI for minima of externl input, 0 for the maxima.
fs = np.arange(50.,450.1,3.);


# Computing the FO onset and offset phases: With STD

rateses = [ [ np.loadtxt("RatesFig8/rates_Esm%d_cm%d_%dHzgChR%d_YSTD.dat" % (int(np.round(-Es)),kcm,int(fsin),int(gsin) ) ) for kcm in range(Nnets) ] for Es in Ess ];
cwtPowses = [ [ compPWT(rateses[kEs][kcm][-Ncycs*Ntcyc:],dt=dt*1.e-3,fs=fs) for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];
max_powers = [ [ np.max(cwtPowses[kEs][kcm]) for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];
kcycs_act = [ [ [ kcyc for kcyc in range(Ncycs) if np.max(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc])>.3*max_powers[kEs][kcm] ] for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];

phmins = [ [ np.min( phases[np.where(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc] > thresh*np.max(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc]))[1]]) for kcm in range(Nnets) for kcyc in kcycs_act[kEs][kcm] ] for kEs in range(len(Ess)) ];
phmaxs = [ [ np.max( phases[np.where(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc] > thresh*np.max(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc]))[1]]) for kcm in range(Nnets) for kcyc in kcycs_act[kEs][kcm] ] for kEs in range(len(Ess)) ];

tmp = [ np.savetxt("phmins_Esm%d_YSTD.dat" % int(np.round(-Ess[kEs])), phmins[kEs]) for kEs in range(len(Ess)) ];
tmp = [ np.savetxt("phmaxs_Esm%d_YSTD.dat" % int(np.round(-Ess[kEs])), phmaxs[kEs]) for kEs in range(len(Ess)) ];



# Computing the FO onset and offset phases: Without STD

rateses = [ [ np.loadtxt("RatesFig8/rates_Esm%d_cm%d_%dHzgChR%d_NSTD.dat" % (int(np.round(-Es)),kcm,int(fsin),int(gsin)) ) for kcm in range(Nnets) ] for Es in Ess ];
cwtPowses = [ [ compPWT(rateses[kEs][kcm][-Ncycs*Ntcyc:],dt=dt*1.e-3,fs=fs) for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];
max_powers = [ [ np.max(cwtPowses[kEs][kcm]) for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];
kcycs_act = [ [ [ kcyc for kcyc in range(Ncycs) if np.max(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc])>.3*max_powers[kEs][kcm] ] for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];

phmins = [ [ np.min( phases[np.where(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc] > thresh*np.max(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc]))[1]]) for kcm in range(Nnets) for kcyc in kcycs_act[kEs][kcm] ] for kEs in range(len(Ess)) ];
phmaxs = [ [ np.max( phases[np.where(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc] > thresh*np.max(cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc]))[1]]) for kcm in range(Nnets) for kcyc in kcycs_act[kEs][kcm] ] for kEs in range(len(Ess)) ];

tmp = [ np.savetxt("phmins_Esm%d_NSTD.dat" % int(np.round(-Ess[kEs])), phmins[kEs]) for kEs in range(len(Ess)) ];
tmp = [ np.savetxt("phmaxs_Esm%d_NSTD.dat" % int(np.round(-Ess[kEs])), phmaxs[kEs]) for kEs in range(len(Ess)) ];
