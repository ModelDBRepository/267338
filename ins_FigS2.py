import numpy as np
import matplotlib.pyplot as plt
import FSIN as sim

#Ises = [np.loadtxt("fIs_exp/IsAct_kneu%d.dat" % kneu) for kneu in range(90)]
#fses = np.array([np.loadtxt("fIs_exp/fsAct_kneu%d.dat" % kneu) for kneu in range(90)])


#DerivadaExp = [(fs[1]-fs[0])/(Is[1]-Is[0]) for (Is,fs) in zip(Ises,fses)]

CMat = 0
Is = np.arange(50.,500.1,50.);#[50.,200.,500.]#
simtime = 800.
Nneur=100

Nrec=100

datas = [ sim.gewnet(Nrec=Nrec,return_state=True,method="rk4",dt=1.e-3,step_current=True,sin_cond=False,Ips=I*np.ones(Nneur),act_syns=False,gjs=True, homo_act=False,kneu=29,sim_time=simtime,mod_gL=True,read_ggaps=True,connectivity_matrix=CMat) for I in Is ];

N=100

ISIses = [ [ sim.b2.diff(datas[kI][1].spike_trains()[kneu]/sim.ms) for kI in range(len(Is)) ] for kneu in range(N) ];
fses = [ [ 1.e+3 / min(ISIses[kneu][kI][-10:]) for kI in range(len(Is)) if len(ISIses[kneu][kI])>15 ] for kneu in range(N) ];
Ises = [ [ Is[kI] for kI in range(len(Is)) if len(ISIses[kneu][kI])>15 ] for kneu in range(N) ];

RheoBaseModel = [Is[0] for Is in Ises]
CutoffModel = np.array([fs[0] for fs in fses])

tmp=plt.figure(figsize=[6.4, 2.5])
tmp =plt.subplot(1,2,1)
tmp = [ plt.plot(Is,fs,"-o",Is[0],fs[0],'rd') for (Is,fs) in zip(Ises,fses) ]; tmp = plt.xlim(0.,500.); tmp = plt.ylim(0.,400.);
#tmp = [ plt.plot(Is,fs,"-") for (Is,fs) in zip(Ises_heteact[::5], fses_heteact[::5]) ]; tmp = plt.xlim(0.,500.); tmp = plt.ylim(0.,500.);
#plt.show();
tmp = plt.savefig("Figures/FigS2A.eps", dpi=300); #tmp = plt.clf(); # To save plot in eps file
tmp = plt.savefig("Figures/FigS2A.png", dpi=300); tmp = plt.clf(); # To save plot in png file



np.savetxt("ModgLIModel.txt",Ises,fmt='%s')
np.savetxt("ModgLfModel.txt",fses,fmt='%s')

