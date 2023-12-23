import FSIN as sim
import numpy as np
import matplotlib.pyplot as plt


method = "euler";
dt = 5.e-4;
sim_time=100.;
kneu=43

data_1pert_shu = sim.gewnet(kneu=kneu,sim_time=sim_time,return_state=True,Nrec=100,dt=dt,method=method,act_syns=True,mod_gL=False, gjs=False,gext=7.,sin_cond=False,const_cond=True,tau_rise=.3,homo_neurons=True,gms=1.65,dist_syns=False,read_delays=False,                   
read_gms=False,dt_rec=.1,dmin=.8,dmax=.8,connectivity="FID",return_syns=False,read_syns=False,USEm=1.,std=False, Es=-55.,randics=False,one_perturbed=True);
data_rics_shu_d16ms = sim.gewnet(kneu=kneu,sim_time=sim_time,return_state=True,Nrec=100,dt=dt,method=method,act_syns=True,mod_gL=False, gjs=False,gext=7.,sin_cond=False,const_cond=True,tau_rise=.3,homo_neurons=True,gms=1.65,dist_syns=False,read_delays=False,                   
read_gms=False,dt_rec=.1,dmin=1.6,dmax=1.6,connectivity="FID",return_syns=False,read_syns=False,USEm=1.,std=False, Es=-55.,randics=True,one_perturbed=False);

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

tmp1 = plt.subplot(3,1,2); tmp = [ plt.plot(data_1pert_shu[1].spike_trains()[kneu]/sim.ms, kneu*np.ones(len(data_1pert_shu[1].spike_trains()[kneu])), "ok", ms=3.) for kneu in range(20) ]; tmp = plt.xlim(0.,100.);
plt.title("B. $E_{syn}$=-55, $\delta$=0.8 ms, perturbed")
tmp1.spines['top'].set_visible(False)
tmp1.spines['right'].set_visible(False)
tmp1.spines['bottom'].set_visible(False)
tmp1.spines['left'].set_visible(False)
tmp1.set_yticks([])
tmp1.set_yticklabels([])
tmp2 = plt.subplot(3,1,3); tmp = [ plt.plot(data_rics_shu_d16ms[1].spike_trains()[kneu]/sim.ms, kneu*np.ones(len(data_rics_shu_d16ms[1].spike_trains()[kneu])), "ok", ms=3.) for kneu in range(20) ]; tmp = plt.xlim(0.,100.);
plt.title("C. $E_{syn}$=-55, $\delta$=1.6 ms, random")
tmp2.spines['top'].set_visible(False)
tmp2.spines['right'].set_visible(False)
tmp2.spines['bottom'].set_visible(False)
tmp2.spines['left'].set_visible(False)
tmp2.set_yticks([])
tmp2.set_yticklabels([])

plt.tight_layout()
plt.savefig("Figures/Fig5BC.png");
plt.savefig("Figures/Fig5BC.eps"); plt.clf();
