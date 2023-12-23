# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import cwt, morlet2
import matplotlib.pyplot as plt
from funcs import *


Nnets = 30;
Ncycs = 30;
fsin = 8
gsin=7.

dt = .1;

Ess = [-75.]#[-75.,-55.];

Ntcyc = int(np.round(125./dt)); # Number of points per theta cycle
phases = np.linspace(-np.pi,+np.pi,Ntcyc); # Values of phases within ones cycle. -PI and +PI for minima of externl input, 0 for the maxima.
fs = np.arange(50.,450.1,3.);

rateses = [ [ np.loadtxt("RatesFig7/rates_Esm%d_cm%d_%dHzgChR%d.dat" % (int(np.round(-Es)),kcm,int(fsin),int(gsin))) for kcm in range(Nnets) ] for Es in Ess ];
cwtPowses = [ [ compPWT(rateses[kEs][kcm][-Ncycs*Ntcyc:],dt=dt*1.e-3,fs=fs) for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];
max_powers = [ [ np.max(cwtPowses[kEs][kcm]) for kcm in range(Nnets) ] for kEs in range(len(Ess)) ];
maxPowses = np.array([ np.array([ np.array([ np.max( cwtPowses[kEs][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc] ) for kcyc in range(Ncycs) ]) for kcm in range(Nnets) ]) for kEs in range(len(Ess)) ]);
kcycs_act_hyp = [ [ kcyc for kcyc in range(Ncycs) if np.max(cwtPowses[0][kcm][:,kcyc*Ntcyc:kcyc*Ntcyc+Ntcyc])>.3*max_powers[0][kcm] ] for kcm in range(Nnets) ];


fpmaxes_hyp = [ [ np.min( fs[np.where( cwtPowses[0][kcm][:,kcyc*Ntcyc:(kcyc+1)*Ntcyc] > (1.-1.e-6)*np.max(cwtPowses[0][kcm][:,kcyc*Ntcyc:(kcyc+1)*Ntcyc])  )[0]] ) for kcyc in kcycs_act_hyp[kcm] ] for kcm in range(Nnets) ];
#fpmaxes_shu = [ [ np.min( fs[np.where( cwtPowses[1][kcm][:,kcyc*Ntcyc:(kcyc+1)*Ntcyc] > (1.-1.e-6)*np.max(cwtPowses[1][kcm][:,kcyc*Ntcyc:(kcyc+1)*Ntcyc])  )[0]] ) for kcyc in range(Ncycs) ] for kcm in range(Nnets) ];

meansfs_hyp = [ np.mean([ fpmaxes_hyp[kcm][kcyc] for kcyc in range(len(kcycs_act_hyp[kcm])) ]) for kcm in range(Nnets) ];
stdsfs_hyp = [ np.std([ fpmaxes_hyp[kcm][kcyc] for kcyc in range(len(kcycs_act_hyp[kcm])) ]) for kcm in range(Nnets) ];

#meansfs_shu = [ np.mean([ fpmaxes_shu[kcm][kcyc] for kcyc in range(Ncycs) ]) for kcm in range(Nnets) ];
#stdsfs_shu = [ np.std([ fpmaxes_shu[kcm][kcyc] for kcyc in range(Ncycs) ]) for kcm in range(Nnets) ];

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

# Plot the maximum wavelet powers
tmp1 = plt.subplot(2,2,1);
tmp = plt.errorbar(range(Nnets),np.mean(maxPowses[0],axis=1), yerr=np.std(maxPowses[0],axis=1), marker="o", color="b", ls="");
tmp = plt.ylim(0.,3.e+6);
tmp1.spines['top'].set_visible(False)
tmp1.spines['right'].set_visible(False)
tmp1.spines['bottom'].set_visible(True)
tmp1.spines['left'].set_visible(True)
tmp = plt.ylabel("$P_{max}$ (spike$^2$/s$^2$)")
tmp1.set_xticks([10,20,30])
tmp1.set_xticklabels([])
tmp1.set_yticks([0.,1.e+6,2.e+6,3.e+6])
tmp1.set_yticklabels(["0","","","3"])

# Plot the frequencies with maximum wavelet power

# For hyperpolarizing synapses
tmp1 = plt.subplot(2,2,2);
tmp = plt.errorbar(range(len(meansfs_hyp)), meansfs_hyp, yerr=stdsfs_hyp, marker="o", color="b", ls="");
tmp = plt.ylim(0.,300.);
tmp1.spines['top'].set_visible(False)
tmp1.spines['right'].set_visible(False)
tmp1.spines['bottom'].set_visible(True)
tmp1.spines['left'].set_visible(True)
tmp = plt.ylabel("$F_{max}$ (Hz)")
tmp = plt.xlabel("network index")
tmp1.set_xticks([10,20,30])
tmp1.set_xticklabels([])
tmp1.set_yticks([0,200,300])
# For shunting synapses
#tmp = plt.subplot(2,2,4);
#tmp = plt.errorbar(range(len(meansfs_shu)), meansfs_shu, yerr=stdsfs_shu, marker="o", color="b", ls="");
#tmp = plt.ylim(0.,450.);

plt.savefig("Figures/Fig7A_%dHzgChR%d.eps" % (int(fsin),int(gsin)) );plt.savefig("Figures/Fig7A_%dHzgChR%d.png" % (int(fsin),int(gsin)) );

