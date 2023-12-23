import numpy as np
from scipy.signal import cwt, morlet2
import matplotlib.pyplot as plt
from funcs import *

"""#Uncomment for redo simulations
Nnets = 15;
Ncycs = 30;
for fsin in [4,8]:
	v=[]
	for gsin in [0.5,1,2,3,4,5,6,7,8,9,10]:

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

		v.append([gsin,np.mean(maxPowses[0],axis=1)[0], np.std(maxPowses[0],axis=1)[0], meansfs_hyp[0], stdsfs_hyp[0] ])

	np.savetxt("PowFreq_%dHz.dat" % int(fsin), v)
"""

es = 0
for fsin in [8,4]:
	Vector = np.loadtxt("PowFreq_%dHz.dat" % int(fsin))       
	tmp = plt.subplot(2,2,1+2*es);
	tmp = plt.title("Freq %d Hz" % int(fsin))
	tmp = plt.errorbar(Vector[:,0],Vector[:,3],yerr=Vector[:,4],fmt="x-");
	tmp = plt.xlabel("gChr (nS)"); tmp = plt.ylabel("Max Amplitude Network Frequency (Hz)");
	tmp = plt.ylim([0,300])
	tmp = plt.subplot(2,2,2+2*es);
	tmp = plt.title("Freq %d Hz" % int(fsin))
	tmp = plt.errorbar(Vector[:,0],Vector[:,1],yerr=Vector[:,2],fmt="x-");
	tmp = plt.xlabel("gChr (nS)"); tmp = plt.ylabel("Max Amplitude Power (u.a.)");
	tmp = plt.ylim([0,2.5e+6])
	es+=1
tmp = plt.tight_layout()

tmp = plt.savefig("Figures/Fig7BC.png")
tmp = plt.savefig("Figures/Fig7BC.eps")
