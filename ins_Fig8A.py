import numpy as np
import matplotlib.pyplot as plt
import FSIN as sim

method = "rk4";
dt = 1.e-3;
fsin = 8
gsin=7.

connectivity_matrix=0; # Index for the connectivity matrix to be used among the ones predefined (loaded from file in the params folder)
Ess = [-75.]#[-75.,-55.];

# Total simulated time. Theta frequency is 8 Hz, i.e. period of 125 ms.
# sim_time given as multiple of period, i.e. number of simulated cycles.
# Only the last two are plotted. To remove effects of transients.
sim_time = 6.*125.; 

datas_hyp = [ sim.gewnet(gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=True,gjs=True,method=method,dt=dt,connectivity_matrix=connectivity_matrix,Es=-75.,std=std) for std in [True,False] ];

spikeses_hyp = [ db[1].t/sim.ms for db in datas_hyp ];

Title = ["A","A"]
SubTitle = ["Phys","Phys"]
STD =["","No"]

Ess1 =[-75,-75]

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

for kcond in range(2):
   tmp1 = plt.subplot(2,1,kcond+1);
   tmp = plt.hist(spikeses_hyp[kcond], np.linspace(sim_time-2.*125.,sim_time,2*125+1), color='k');
   tmp = plt.xlim(sim_time-2.*125.,sim_time); #tmp = plt.ylim(0.,85.);
   tmp1.spines['top'].set_visible(False)
   tmp1.spines['right'].set_visible(False)
   tmp1.spines['bottom'].set_visible(False)
   tmp1.spines['left'].set_visible(True)
   tmp1.set_xticks([])
   tmp1.set_xticklabels([])
   tmp1.set_yticks([0,40,80])
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s1. %s STD ($E_{syn}$=%d), %s. $g_{GAP}$" % (Title[kcond],STD[kcond],Ess1[kcond],SubTitle[kcond]),fontsize=10)




#tmp1.set_ylabel("")
#tmp3.set_ylabel("")
#tmp1.set_yticks([0,40,80])
#tmp1.set_yticklabels([])
#tmp3.set_yticks([0,40,80])
#tmp3.set_yticklabels([])
tmp = plt.tight_layout(pad=0.0, w_pad=1.5, h_pad=1.0)

tmp = plt.savefig("Figures/Fig8A_%dHzgChR%d.eps" % (int(fsin),int(gsin)), dpi=300); #tmp = plt.clf(); # To save plot in eps file
tmp = plt.savefig("Figures/Fig8A_%dHzgChR%d.png" %  (int(fsin),int(gsin)), dpi=300); tmp = plt.clf(); # To save plot in png file
