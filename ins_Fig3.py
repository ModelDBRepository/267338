import numpy as np
import matplotlib.pyplot as plt
import FSIN as sim
import brian2 as b2
from scipy.signal import find_peaks

method = "euler"; # methods "rk4" and "gsl_rk8pd" give similar results
dt = 5.e-4;


Ess = [-75.,-55.];
gmsSim = 1.65 
# Total simulated time. Theta frequency is 8 Hz, i.e. period of 125 ms.
# sim_time given as multiple of period, i.e. number of simulated cycles.
# Only the last two are plotted. To remove effects of transients.
sim_time = 6.*125.; 
kneu=33
width=.55
ConnMatrx = 0
gsin=7.
fsine=8
v_thresh=-30.
dt_rec=dt

factorTau=1

datas_homogeneous = [ sim.gewnet(factorTau=factorTau,dt_rec=dt_rec,v_thresh=v_thresh,return_state=True,sim_time=sim_time, mod_gL=False,gjs=False,homo_neurons=True,kneu=kneu,dist_syns=False,read_syns=False,read_delays=False,read_gms=False,dmin=.8,dmax=.8,gms=gmsSim,USEm=1.,connectivity="FID",method=method,dt=dt,Es=Es,window_width=width*b2.ms) for Es in Ess ];
rate_homogeneous = [ db[0] for db in datas_homogeneous ];
spikeses_homogeneous = [ db[1].t/sim.ms for db in datas_homogeneous ];
spikeses_homogeneousIdx = [ db[1].i for db in datas_homogeneous ];

datas_heterogeneous_all = [ sim.gewnet(factorTau=factorTau,dt_rec=dt_rec,v_thresh=v_thresh,return_state=True, sim_time=sim_time,gms=gmsSim,mod_gL=False,gjs=False,connectivity_matrix=ConnMatrx, method=method,dt=dt,Es=Es,window_width=width*b2.ms) for Es in Ess ];
spikeses_heterogeneous_all = [ db[1].t/sim.ms for db in datas_heterogeneous_all ];

datas_heterogeneous_syns = [ sim.gewnet(factorTau=factorTau,dt_rec=dt_rec,v_thresh=v_thresh,gsin=gsin,return_state=True, sim_time=sim_time,mod_gL=False,gms=gmsSim,gjs=False,homo_neurons=True,kneu=kneu,connectivity_matrix=ConnMatrx, method=method,dt=dt,Es=Es,window_width=width*b2.ms) for Es in Ess ];
spikeses_heterogeneous_syns = [ db[1].t/sim.ms for db in datas_heterogeneous_syns ];

Npoints = int(2*125/2.)

bin_edges = np.linspace(sim_time-2.*125.,sim_time,Npoints,endpoint=True)

Title = ["A","B"]

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

for kEs in range(len(Ess)):
   tmp1 = plt.subplot(4,2,kEs+1);
   tmp1.spines['top'].set_visible(False)
   tmp1.spines['right'].set_visible(False)
   tmp1.spines['bottom'].set_visible(False)
   tmp1.spines['left'].set_visible(True)
   tmp1.set_xticks([])
   tmp1.set_xticklabels([])
   tmp1.set_yticks([0,50,100])
   tmp = plt.hist(spikeses_homogeneous[kEs], bin_edges, color = 'k');
   tmp = plt.xlim(sim_time-2.*125.,1.*sim_time); tmp = plt.ylim(0.,105.);
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s1. Homogeneous Network ($E_{syn}$=%d)" % (Title[kEs],Ess[kEs]),fontsize=10)

   tmp2 = plt.subplot(4,2,kEs+3);
   tmp2.spines['top'].set_visible(False)
   tmp2.spines['right'].set_visible(False)
   tmp2.spines['bottom'].set_visible(False)
   tmp2.spines['left'].set_visible(True)
   tmp2.set_xticks([])
   tmp2.set_xticklabels([])
   tmp2.set_yticks([0,50])
   tmp = plt.hist(spikeses_heterogeneous_all[kEs], bin_edges, color = 'k');
   tmp = plt.xlim(sim_time-2.*125.,1.*sim_time); tmp = plt.ylim(0.,70.);
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s2. Intrinsic Heterogeneity ($E_{syn}$=%d)" % (Title[kEs],Ess[kEs]),fontsize=10)

   tmp3 = plt.subplot(4,2,kEs+5);
   tmp3.spines['top'].set_visible(False)
   tmp3.spines['right'].set_visible(False)
   tmp3.spines['bottom'].set_visible(False)
   tmp3.spines['left'].set_visible(True)
   tmp3.set_yticks([0,50,100])
   tmp = plt.hist(spikeses_heterogeneous_syns[kEs], bin_edges, color = 'k');
   tmp = plt.xlim(sim_time-2.*125.,1.*sim_time); tmp = plt.ylim(0.,105.);
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s3. Synaptic Heterogeneity ($E_{syn}$=%d)" % (Title[kEs],Ess[kEs]),fontsize=10)

Npoints = int(sim_time/(dt_rec)) 
time = np.linspace(0.,sim_time,Npoints,endpoint=True)
tmp4 = plt.subplot(4,2,7);
plt.plot(time,gsin*(1.+np.sin(2.*np.pi*fsine*time/1000-np.pi/2.)),'k');tmp = plt.xlim(sim_time-2.*125.,1.*sim_time);
tmp4 = plt.subplot(4,2,8);
plt.plot(time,gsin*(1.+np.sin(2.*np.pi*fsine*time/1000-np.pi/2.)),'k');tmp = plt.xlim(sim_time-2.*125.,1.*sim_time);


tmp1.set_ylabel("")
tmp2.set_ylabel("")
tmp3.set_ylabel("")
tmp1.set_yticks([0,50,100])
tmp1.set_yticklabels([])
tmp2.set_yticks([0,50])
tmp2.set_yticklabels([])
tmp3.set_yticks([0,50,100])
tmp3.set_yticklabels([])
tmp = plt.tight_layout(pad=0.0, w_pad=2.0, h_pad=1.0)
tmp=plt.savefig("Figures/Fig3gChR%d.png" % gsin);tmp=plt.savefig("Figures/Fig3gChR%d.eps" % gsin); tmp=plt.clf();
