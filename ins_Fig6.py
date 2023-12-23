import numpy as np
import matplotlib.pyplot as plt
import FSIN as sim
import os
import brian2 as b2

TigerFish = True
if TigerFish: 	# For Tigerfish use only
	cache_dir = os.path.expanduser('~/.cython/brian-pid-{}'.format(os.getpid()))
	b2.prefs.codegen.runtime.cython.cache_dir = cache_dir
	b2.prefs.codegen.runtime.cython.multiprocess_safe = False


method = "rk4"; # methods "rk4" and "gsl_rk8pd" give similar results
dt = 4.e-3;
factor=1;
gm_scale = factor;
tau_fall=2.0
tau_rise=.3
FactorScaleKV3=1.



connectivity_matrix=0; # Index for the connectivity matrix to be used among the ones predefined (loaded from file in the params folder)
Ess = [-75.,-55.];

fsin=8.;
gsin=7.;

std=True

# Total simulated time. Theta frequency is 8 Hz, i.e. period of 125 ms.
# sim_time given as multiple of period, i.e. number of simulated cycles.
# Only the last two are plotted. To remove effects of transients.
sim_time = 6.*125.; 
dt_rec = 0.01

datas_nogjs = [ sim.gewnet(std=std,FactorScaleKV3=FactorScaleKV3,tau_fall=tau_fall,read_gms=False,scale_gm=False,gm_scale = gm_scale,read_delays=False,gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=False,gjs=False,method=method,dt=dt,connectivity_matrix=connectivity_matrix,Es=Es) for Es in Ess ];
datas_physiogjs = [ sim.gewnet(std=std, FactorScaleKV3=FactorScaleKV3,dt_rec=dt_rec,return_state=True,tau_fall=tau_fall,tau_rise=tau_rise,read_gms=False,scale_gm=False,gm_scale = gm_scale,read_delays=False,gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=True,gjs=True,method=method,dt=dt,connectivity_matrix=connectivity_matrix,Es=Es) for Es in Ess ];

#Npoints = int(sim_time/(dt_rec)) # Timestep 2. ms
#time = np.linspace(0.,sim_time,Npoints,endpoint=True)
#V = [ db[2].v[0]/sim.ms for db in datas_physiogjs ];
#plt.plot(time,V[0])
#tmp = plt.xlim(sim_time-2.*125.,sim_time);
#tmp =plt.yticks([-80,-75,-70,-60,-50,-40,-30,-20,-10,0])
#plt.savefig("20ChRModifiedTauVoltageTraceConn%s.png" % Ess[0]); plt.savefig("20ChRModifiedTauVoltageTraceConn%s.eps" % Ess[0])
#plt.clf()
#quit()
datas_unphysiogjs = [ sim.gewnet(std=std,FactorScaleKV3=FactorScaleKV3,tau_fall=tau_fall,tau_rise=tau_rise,read_gms=False,scale_gm=False,gm_scale = gm_scale,read_delays=False,gsin=gsin,fsin=fsin,sim_time=sim_time,mod_gL=False,gjs=True,read_ggaps=False,ggaps=2.,method=method,dt=dt,connectivity_matrix=connectivity_matrix,Es=Es) for Es in Ess ];

spikeses_nogjs = [ db[1].t/sim.ms for db in datas_nogjs ];
spikeses_physiogjs = [ db[1].t/sim.ms for db in datas_physiogjs ];
spikeses_unphysiogjs = [ db[1].t/sim.ms for db in datas_unphysiogjs ];


Title = ["A","B"]

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

for kEs in range(len(Ess)):
   tmp1 = plt.subplot(3,2,kEs+1);
   tmp = plt.hist(spikeses_nogjs[kEs], np.linspace(sim_time-2.*125.,sim_time,2*125+1),color = 'k');
   tmp = plt.xlim(sim_time-2.*125.,sim_time); tmp = plt.ylim(0.,55.);
   tmp1.spines['top'].set_visible(False)
   tmp1.spines['right'].set_visible(False)
   tmp1.spines['bottom'].set_visible(False)
   tmp1.spines['left'].set_visible(True)
   tmp1.set_xticks([])
   tmp1.set_xticklabels([])
   tmp1.set_yticks([0,50])
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s1. No Gap Junctions ($E_{syn}$=%d)" % (Title[kEs],Ess[kEs]),fontsize=10)

   tmp2 = plt.subplot(3,2,kEs+3);
   tmp = plt.hist(spikeses_physiogjs[kEs], np.linspace(sim_time-2.*125.,sim_time,2*125+1),color = 'k');
   tmp = plt.xlim(sim_time-2.*125.,sim_time); tmp = plt.ylim(0.,85.);
   tmp2.spines['top'].set_visible(False)
   tmp2.spines['right'].set_visible(False)
   tmp2.spines['bottom'].set_visible(False)
   tmp2.spines['left'].set_visible(True)
   tmp2.set_xticks([])
   tmp2.set_xticklabels([])
   tmp2.set_yticks([0,50])
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s1. Physiological $g_{GAP}$ ($E_{syn}$=%d)" % (Title[kEs],Ess[kEs]),fontsize=10)


   tmp3 = plt.subplot(3,2,kEs+5);
   tmp3.spines['top'].set_visible(False)
   tmp3.spines['right'].set_visible(False)
   tmp3.spines['bottom'].set_visible(False)
   tmp3.spines['left'].set_visible(True)
   tmp3.set_yticks([0,50,100])
   tmp = plt.hist(spikeses_unphysiogjs[kEs], np.linspace(sim_time-2.*125.,sim_time,2*125+1),color = 'k');
   tmp = plt.xlim(sim_time-2.*125.,sim_time); tmp = plt.ylim(0.,105.);
   tmp = plt.ylabel("Counts",fontsize=11)
   tmp = plt.title("%s3. Non-physiological $g_{GAP}$ ($E_{syn}$=%d)" % (Title[kEs],Ess[kEs]),fontsize=10)


tmp1.set_ylabel("")
tmp2.set_ylabel("")
tmp3.set_ylabel("")
tmp1.set_yticks([0,50])
tmp1.set_yticklabels([])
tmp2.set_yticks([0,40,80])
tmp2.set_yticklabels([])
tmp3.set_yticks([0,50,100])
tmp3.set_yticklabels([])
tmp = plt.tight_layout(pad=0.0, w_pad=2.0, h_pad=1.0)

tmp = plt.savefig("Figures/Fig6_f{f}gChR{g}.eps".format(f=fsin,g=gsin), dpi=300); #tmp = plt.clf(); # To save plot in eps file
tmp = plt.savefig("Figures/Fig6_f{f}gChR{g}.png".format(f=fsin,g=gsin), dpi=300); tmp = plt.clf(); # To save plot in png file



