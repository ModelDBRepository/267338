import numpy as np
import matplotlib.pyplot as plt
import FSIN as sim

dt_rec = .1;
Is = np.arange(300.,425.1,25.);

dataTraces = [ sim.gewnet(sim_time=800.,step_current=True,sin_cond=False,gjs=False,act_syns=False,Ips=I*np.ones(100),N=100,Nrec=100,mod_gL=False,dt_rec=dt_rec,return_state=True) for I in Is ];

tmin = 370.;
tmax = 470.;

kneu = 34
kI0 = 0; NIs = 6;
width = 5.2/2 # in inches
height = 7.*2/3.# in inches
tmp = plt.figure(figsize=[width,height]);

for kI in range(NIs):
 tmp1 = plt.subplot(NIs,1,NIs-kI);
 tmp = plt.plot( np.arange(tmin,tmax,dt_rec), dataTraces[kI+kI0][-1].v[kneu][int(np.round(tmin/dt_rec)):int(np.round(tmax/dt_rec))]/sim.mV , "k" );
 tmp = plt.xlim(tmin,tmax);
 tmp = plt.ylim(-90.,10.);
 if kI!=0: tmp = plt.xticks([]); tmp = plt.yticks([]);
 tmp1.spines['top'].set_visible(False)
 tmp1.spines['right'].set_visible(False)
 tmp1.spines['bottom'].set_visible(False)
 tmp1.spines['left'].set_visible(False)

tmp = plt.savefig("Figures/Fig1A2.eps", dpi=300); #tmp = plt.clf(); # To save plot in eps file
tmp = plt.savefig("Figures/Fig1A2.png", dpi=300); tmp = plt.clf(); # To save plot in png file
