import FSIN as sim
import numpy as np
import matplotlib.pyplot as plt
import sys

method = "euler";
dt = 5.e-4;
sim_time=100.;

gms=1.65

kneu = 43

data_rics_hyp = sim.gewnet(kneu=kneu,sim_time=sim_time,return_state=True,Nrec=100,dt=dt,method=method,act_syns=True,mod_gL=False, gjs=False,gext=7.,sin_cond=False,const_cond=True,tau_rise=.3,homo_neurons=True,gms=gms,dist_syns=False,read_delays=False,                  
read_gms=False,dt_rec=.1,dmin=.8,dmax=.8,connectivity="FID",return_syns=False,read_syns=False,USEm=1.,std=False, Es=-75.,randics=True,one_perturbed=False);

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

tmp1 = plt.subplot(3,1,1); tmp = [ plt.plot(data_rics_hyp[1].spike_trains()[kneu]/sim.ms, kneu*np.ones(len(data_rics_hyp[1].spike_trains()[kneu])), "ok", ms=3.) for kneu in range(20) ]; tmp = plt.xlim(0.,100.);

ISIses_heteact = sim.b2.diff(data_rics_hyp[1].spike_trains()[0]/sim.ms);

print(ISIses_heteact[-1]-ISIses_heteact[-2])

plt.title("B. $E_{syn}$=-75, $\delta$=0.8 ms, random")
tmp1.spines['top'].set_visible(False)
tmp1.spines['right'].set_visible(False)
tmp1.spines['bottom'].set_visible(False)
tmp1.spines['left'].set_visible(False)
tmp1.set_yticks([])
tmp1.set_yticklabels([])
plt.savefig("Figures/Fig4B.png");
plt.savefig("Figures/Fig4B.eps"); plt.clf();
