import numpy as np
import matplotlib.pyplot as plt
import FSIN as sim

Ises = [np.loadtxt("fIs_exp/IsAct_kneu%d.dat" % kneu) for kneu in range(90)]
fses = np.array([np.loadtxt("fIs_exp/fsAct_kneu%d.dat" % kneu) for kneu in range(90)])

#RheoBaseExp = [Is[0] for Is in Ises]
#CutoffExp = np.array([fs[0] for fs in fses])
#DerivadaExp = [(fs[1]-fs[0])/(Is[1]-Is[0]) for (Is,fs) in zip(Ises,fses)]

Is = np.arange(50.,500.1,25.);

datas_heteact = [ sim.gewnet(dt=5.e-4,step_current=True,sin_cond=False,Ips=I*np.ones(100),act_syns=False,gjs=False,homo_act=False,kneu=29,sim_time=800.,mod_gL=False) for I in Is ];


ISIses_heteact = [ [ sim.b2.diff(datas_heteact[kI][1].spike_trains()[kneu]/sim.ms) for kI in range(len(Is)) ] for kneu in range(100) ];
fses_heteact = [ [ 1.e+3 / ISIses_heteact[kneu][kI][-1] for kI in range(len(Is)) if len(ISIses_heteact[kneu][kI])>15 ] for kneu in range(100) ];
Ises_heteact = [ [ Is[kI] for kI in range(len(Is)) if len(ISIses_heteact[kneu][kI])>15 ] for kneu in range(100) ];

#RheoBaseModel = [Is[0] for Is in Ises_heteact]
#CutoffModel = [fs[0] for fs in fses_heteact]

#IndexCutoffFreq = np.argsort(np.array(CutoffModel))

#ISIses_heteact = [ISIses_heteact[Index][:] for Index in IndexCutoffFreq]
#fses_heteact = [fses_heteact[Index][:] for Index in IndexCutoffFreq]

width = 5.2 # in inches
height = 7./3# in inches
tmp = plt.figure(figsize=[width,height]);

tmp = plt.subplot(1,2,1); tmp = [ plt.plot(Is,fs,"-o",Is[0],fs[0],"rd") for (Is,fs) in zip(Ises,fses) ]; tmp = plt.xlim(0.,500.); tmp = plt.ylim(0.,400.);
#print([ plt.plot(Is,fs,"-o",Is[0],fs[0],"rd") for (Is,fs) in zip(Ises_heteact, fses_heteact) ])
#quit()
tmp = plt.subplot(1,2,2); tmp = [ plt.plot(Is,fs,"-o",Is[0],fs[0],"rd") for (Is,fs) in zip(Ises_heteact, fses_heteact) ]; tmp = plt.xlim(0.,500.); tmp = plt.ylim(0.,400.);

tmp = plt.savefig("Figures/Fig1B.eps", dpi=300);
tmp = plt.savefig("Figures/Fig1B.png", dpi=300); tmp = plt.clf(); # To save plot in png file


#np.savetxt("IModel.txt",ISIses_heteact,fmt='%s')
#np.savetxt("fModel.txt",fses_heteact,fmt='%s')

