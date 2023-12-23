# Code adapted by Guillem Via from the original code implementing the Brunel (2000) network model from the exercise repository for
# the book by Wulfram Gerstner et al. (2014)

# Single model neuron equations are extracted in parts by Erisir et al, (1999), Golomb et al. (2007) and Economo and White (2012).
# Parameters are calibrated to reproduce passive properties and f/I curve heterogeneity from a dataset of 11 parvalbumin positive
# fast-spiking interneurons recorded in mouse dorsal Medial Entorhinal Cortex slices.



# Original description of the source code is as follows:

# This file is part of the exercise code repository accompanying
# the book: Neuronal Dynamics (see http://neuronaldynamics.epfl.ch)
# located at http://github.com/EPFL-LCN/neuronaldynamics-exercises.

# This free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License 2.0 as published by the
# Free Software Foundation. You should have received a copy of the
# GNU General Public License along with the repository. If not,
# see http://www.gnu.org/licenses/.

# Should you reuse and publish the code for your own purposes,
# please cite the book or point to the webpage http://neuronaldynamics.epfl.ch.

# Wulfram Gerstner, Werner M. Kistler, Richard Naud, and Liam Paninski.
# Neuronal Dynamics: From Single Neurons to Networks and Models of Cognition.
# Cambridge University Press, 2014.

import brian2 as b2
from brian2 import NeuronGroup, Synapses, PoissonInput, PoissonGroup, SpikeGeneratorGroup, Equations, TimedArray
from brian2.monitors import StateMonitor, SpikeMonitor, PopulationRateMonitor
from brian2 import Hz, nS, nF, Mohm, pA, mV, ms
import matplotlib.pyplot as plt
from time import time
import numpy as np
from numpy import random, exp, sqrt, log, pi
from random import shuffle, sample
import os
TigerFish = True
if TigerFish: 	# For Tigerfish use only
	cache_dir = os.path.expanduser('~/.cython/brian-pid-{}'.format(os.getpid()))
	b2.prefs.codegen.runtime.cython.cache_dir = cache_dir
	b2.prefs.codegen.runtime.cython.multiprocess_safe = False

def gewnet(
    eps=1.e-9, # Small parameter to avoid divide by zero in gating variable activation function
    connectivity_matrix=0, # To load connectivity matrices for chemical and electrical synapses from different files
    sim_time=1500., # Total simulated time (in ms)
    N=100, # Number of neurons
    Nrec=10, # Number of neurons to record state variables from (only if return_state==True)
    v_thresh=-30., # Threshold membrane potential whose crossing counts a spike and triggers a postsynaptic IPSC
    epsilon=.36, # Connection probability for chemical synapses. Used only if connectivity matrix is not read from file
    I=0., # External input current in the constant current protocol, i.e. if const_current==True (in pA)
    gext=7., # In case external drive is applied as a constant conductance, multiplied by a driving force with reversal Ese (in nS)
    Ips=0., # Value or array of values for the Ip applied to each neuron. If array, then need to be numpy.array() type.
    Isin=7.,
    gsin=7., # In case external drive is applied as a conductance (times a driving force). Conductance oscillates between 0 and 2*gsin.
    fsin=8.,
    pgj=.27,
    ggaps=1.,
    gLmin=1.5,
    dmin=.6,
    dmax=1.,
    gms=1.65, # Value for peak conductances if not heterogeneous, 1.65 by defaut, i.e. mean of the log-normal used for most simulations (in nS)
    gmmean=0., # if we assume Esyn!=-75, gms need to be rescaled by a factor a, i.e. lognormal dist mu->mu+log(a)
    gmstd=1., # If lognormal gm, its std
    gm_scale=1., # If scale-gm==True, scale peak conductance values gms by this factor gm_scale
    tau_fall=2.0, # Synaptic decay time constant. For the bi-exponential IPSC conductance
    tau_rise=.3, # Synaptic rise time constant. For the bi-exponential IPSC conductance
    Es=-75., # Reversal potential for the chemical synapses
    Ese=0.,  # Synaptic reversal for excitation. In case external through conductance.
    taudm=100.,
    USEm=.3,
    gm_distribution="lognormal",
    homo_neurons=False,
    homo_act=False,
    kneu=33,
    const_current=False,
    const_cond=False,
    step_current=False,
    sin_cond=True,
    sin_current=False,
    randics=True,
    one_perturbed=False,
    act_syns=True,
    std=True,
    scale_gm=False,
    gjs=True,
    mod_gL=True,
    read_syns=True,
    dist_syns=True,
    read_gms=True,
    read_delays=True,
    read_gjs=True,
    read_ggaps=True,
    return_state=False,
    return_syns=False,
    connectivity="ER",
    window_width=.5*b2.ms,
    windowType='flat',
    dt=.01,
    dt_rec=.1,
    factorTau=1,
    FactorScaleKV3=1,
    method="euler"):


    st=time(); # Start timer to compute simulation running time
    #b2.start_scope()

    b2.defaultclock.dt = dt * ms # Set brian simulator default clock


    ### LOADING PARAMETERS AND SETTING UNITS

    # Set fixed and homogeneous neuron parameters

    ENa=+50.*mV; EK=-90.*mV;

    km2=.1; kh1=.012; kh2=.2; kn2=FactorScaleKV3*.001; ka2=.02;  kn1 = FactorScaleKV3*1.;

    sigm1=4.*mV; sigh2=3.5*mV; sign1=12.*mV; sigm2=-13.*mV; sigh1=-20.*mV; sign2=-8.5*mV; siga1=12.*mV; siga2=-80.*mV;





    # Tuned parameters with heterogeneity and properties as observed in electrophysiological recordings

    # Load neuron parameters that can be heterogeneous (need to set folder and file names)
    params_act = np.transpose(np.loadtxt( "params/intrinsic/active/params_actRoman.dat" )); # Heterogeneous active parameters: Voltage-gated current peak conductances and thetas
    #params_act = np.loadtxt( "params/intrinsic/active/params_act.dat" ); # Heterogeneous active parameters: Voltage-gated current peak conductances and thetas
    if mod_gL==True:
        folder = "mod_gL/sigma04_higg12/gLmin%d/" % int(np.round(10.*gLmin));
        gLs = np.loadtxt( "params/intrinsic/passive/gLs/" + folder + "gLs%d.dat" % connectivity_matrix )[:N]; # Leakage conductances
        ELs = np.loadtxt( "params/intrinsic/passive/ELs/" + folder + "ELs%d.dat" % connectivity_matrix )[:N]; # Leakage reversal potentials
    else:
        folder = "original/";
        gLs = np.loadtxt( "params/intrinsic/passive/gLs/" + folder + "gLsRoman.dat" )[:N]; # Leakage conductances
        ELs = np.loadtxt( "params/intrinsic/passive/ELs/" + folder + "ELsRoman.dat" )[:N]; # Leakage reversal potentials
    taums = factorTau*np.loadtxt( "params/intrinsic/passive/taums/taumsRomanUniform.dat" )[:N]; # Membrane time constants

    Rins = 1.0e+3/gLs;
    Caps = taums/Rins;

    gNas, gKv3s, gKv1s, thm1s, thh2s, thn1s, tha1s = params_act[:,:N];

    #gLs = np.asarray([14.705184211981187])
    #ELs = np.asarray([-71.97508650437244])
    #taums = np.asarray([5.222681687884308310])#np.asarray([5.22])
    #Rins = 1.0e+3/gLs;
    #Caps = taums/Rins;
    #gNas, gKv3s, gKv1s, thm1s, thh2s, thn1s, tha1s = np.asarray([1.68049678e+04]), np.asarray([6.31700393e+02]), np.asarray([5.90431313e+01]),  np.asarray([-5.29951244e+01]), np.asarray([-5.57111801e+01]), np.asarray([5.90077045e+00]), np.asarray([5.13553887e+01])


    # SET HOMOGENEOUS INTRINSIC NEURONAL PARAMETERS IF DESIRED (ONLY FOR ACTIVE OR BOTH FOR ACTIVE AND PASSIVE)

    # If considering a homogeneous network, take passive parameters from a single model neuron
    if homo_neurons==True: gLs, ELs, taums, Rins, Caps = gLs[kneu]*np.ones(N), ELs[kneu]*np.ones(N), taums[kneu]*np.ones(N), Rins[kneu]*np.ones(N), Caps[kneu]*np.ones(N);

    # If considering a network with homogeneous intrinsic neuronal active parameters (from kneu set in the function call)
    if homo_neurons==True or homo_act==True: gNas, gKv3s, gKv1s, thm1s, thh2s, thn1s, tha1s = gNas[kneu]*np.ones(N), gKv3s[kneu]*np.ones(N), gKv1s[kneu]*np.ones(N), thm1s[kneu]*np.ones(N), thh2s[kneu]*np.ones(N), thn1s[kneu]*np.ones(N), tha1s[kneu]*np.ones(N);
    #print(kneu)
    #print(gLs[kneu], ELs[kneu], taums[kneu], Rins[kneu], Caps[kneu],gNas[kneu], gKv3s[kneu], gKv1s[kneu], thm1s[kneu], thh2s[kneu], thn1s[kneu], tha1s[kneu]); 
    #quit()


    # Import synaptic parameters

    # Chemical synapses
    if act_syns==True:

        if read_syns==True: syns = np.loadtxt( "params/chem_syns/conn_mat/p036/syns%d.dat" % connectivity_matrix, dtype=type(1))[:N];

        if read_delays==True: delays = np.loadtxt( "params/chem_syns/delays/p036/delays%d.dat" % connectivity_matrix)[:len(syns[0])];

        if read_gms==True: gms = np.loadtxt( "params/chem_syns/gms/p036/" + gm_distribution + "/gms%d.dat" % connectivity_matrix)[:len(syns[0])]; # Read peak conductances

    # Gap junctions
    if gjs==True:
        if read_gjs==True: synsgj = np.loadtxt( "params/gjs/conn_mat/" + folder + "synsgj%d.dat" % connectivity_matrix).astype(int);
        if read_ggaps==True: ggaps = np.loadtxt( "params/gjs/ggaps/" + folder + "ggaps%d.dat" % connectivity_matrix)[:len(synsgj[0])];





    # Set units

    I = I*pA; Ips = Ips*pA; Isin = Isin*pA; # External input (depending on protocol)

    gext = gext*nS; gsin = gsin*nS; Ese = Ese*mV; # External conductance (depending on protocol - emulates Chr optogenetic input - reversal of 0 mV)

    gNas = gNas*nS; gKv3s = gKv3s*nS; gKv1s = gKv1s*nS; # Peak conductances for voltage-gated currents

    thm1s = thm1s*mV; thh2s = thh2s*mV; thn1s = thn1s*mV; tha1s = tha1s*mV; eps = eps*mV;

    taums = taums*ms; gLs = gLs*nS; Caps = Caps*nF; Rins = Rins*Mohm; ELs = ELs*mV; # Passive properties
    v_thresh = v_thresh*mV;
    #print("gLs,ELs,taums,gNas, gKv3s, gKv1s, thm1s, thh2s, thn1s, tha1s Caps")
    #print(gLs[kneu],ELs[kneu],taums[kneu],gNas[kneu], gKv3s[kneu], gKv1s[kneu], thm1s[kneu], thh2s[kneu], thn1s[kneu], tha1s[kneu],Caps[kneu])
    #quit()  
    # Synaptic parameters
    if act_syns==True:
        tau_fall = tau_fall*ms; tau_rise = tau_rise*ms; Es = Es*mV; gms = gms*nS;
        if read_delays==True: delays = delays*ms;

    ggaps = ggaps*nS; # Gap junction peak conductances

    # Theta drive frequency
    fsin = fsin*Hz;

    # Short-term-depression recorvery time constant
    taudm = taudm*ms;




    # Normalizing factor for the IPSC conductance height
    c_fall = 1./tau_fall; c_rise = 1./tau_rise;
    norm_syn = 1./( exp(-c_fall*log(c_rise/c_fall)/(c_rise-c_fall)) - exp(-c_rise*log(c_rise/c_fall)/(c_rise-c_fall)) );
    









    ### SETTING EXTERNAL STIMULI

    # Applied step current (activated at simulation middle time)
    if step_current==True: intervals = [0.,1.]; stimulus = TimedArray(intervals,dt=sim_time/2.*b2.ms);







    ### SETTING NEURON EQUATIONS



    # Membrane-potential dynamics equations

    gat_vars =    """
                  alpham = -((v-thm1-eps)/sigm1)/(exp(-(v-thm1-eps)/sigm1)-1.)/ms : Hz
                  betam = km2*exp(v/sigm2)/ms : Hz
                  dm/dt = alpham*(1.-m) - betam*m : 1
                  alphah = kh1*exp(v/sigh1)/ms : Hz
                  betah = -kh2*(v-thh2)/(exp(-(v-thh2)/sigh2)-1.)/(ms*mV) : Hz
                  dh/dt = alphah*(1.-h) - betah*h : 1
                  alphan = -kn1*(v-thn1)/(exp(-(v-thn1)/sign1)-1.)/(ms*mV) : Hz
                  betan = kn2*exp(v/sign2)/ms : Hz
                  dn/dt = alphan*(1.-n) - betan*n : 1
                  alphaa = -(v-tha1)/(exp(-(v-tha1)/siga1)-1.)/(ms*mV) : Hz
                  betaa = ka2*exp(v/siga2)/ms : Hz
                  da/dt = alphaa*(1.-a) - betaa*a : 1
                  INa = - gNa*m*m*m*h*(v-ENa) : ampere
                  IKv3 = - gKv3*n*n*n*n*(v-EK) : ampere
                  IKv1 = - gKv1*a*a*a*a*(v-EK) : ampere
                  IL = - gL*(v-EL) : ampere
                  """

    # Heterogeneous parameters

    param_eqs =   """
                  Ip : ampere
                  EL : volt
                  taum : second
                  gL : siemens
                  Cap : farad
                  gNa : siemens
                  gKv3 : siemens
                  gKv1 : siemens
                  thm1 : volt
                  thh2 : volt
                  thn1 : volt
                  tha1 : volt
                  """

    # Gap-junction currents

    if gjs == True:
        param_eqs += """
                  Igap : ampere
                  """


    # Chemical synapse currents

    if act_syns == True:
        syn_eqs = """
                     Isyn = - (v-Es)*norm_syn*(g_fall-g_rise) : ampere
                     dg_fall/dt = -g_fall/tau_fall : siemens
                     dg_rise/dt = -g_rise/tau_rise : siemens
                     """
    else: syn_eqs = "";



    # External currents (conditional on the protocol: constant, step or theta-modulated; conductance or current-based)

    ext_eqs = "Iext = "

    if const_cond==True: ext_eqs += " - gext*(v-Ese) "
    if const_current==True: ext_eqs += " + I "
    if step_current==True: ext_eqs += " + Ip*stimulus(t) "
    if sin_cond==True: ext_eqs += " - gsin*(1.+sin(2.*pi*fsin*t-pi/2.))*(v-Ese) "
    if sin_current==True: ext_eqs += " + Isin*(1.+sin(2.*pi*fsin*t-pi/2.)) "	

    ext_eqs +=   " : ampere \n "








    # Master equation adding all above currents

    currs = " INa + IKv1 + IKv3 + IL + Iext "
    if act_syns == True: currs += " + Isyn "
    if gjs == True: currs += " + Igap "

    volt_eqs = "\n dv/dt = ( " + currs + " )/Cap : volt \n"

    equations = gat_vars + ext_eqs + syn_eqs + param_eqs + volt_eqs;

    print("Neuron equations are:\n %s\n\n" % ( equations ) );








    # Define neuron population. Use a threshold v_thresho to count a spike and trigger postsynaptic responses.
    # Since the spike has non-zero width a refractory period prevents the code from counting more than one spike while v>v_thresb
    refractory="v>=v_thresh"
    net = NeuronGroup(N, equations, threshold="v>v_thresh", refractory=refractory, method=method)

    net.Ip = Ips;

    net.EL, net.gL, net.taum, net.Cap = ELs, gLs, taums, Caps;
    net.gNa, net.gKv3, net.gKv1, net.thm1, net.thh2, net.thn1, net.tha1 = gNas, gKv3s, gKv1s, thm1s, thh2s, thn1s, tha1s;

    # Initialize state variable at random values
    if randics==True:
        net.v = ELs+3.*mV*np.random.randn(len(ELs)); # Normally distributed membrane potential
        net.m = .01+.1*np.random.rand(len(ELs)); # Uniformly distributed gating variables
        net.h = .99-.1*np.random.rand(len(ELs));
        net.n = .01+.1*np.random.rand(len(ELs));
        net.a = .01+.1*np.random.rand(len(ELs));
    else:
        net.v = -30.*mV;
        net.m = .465391;
        net.h = .184172;
        net.n = .582671;
        net.a = .399927;
    if one_perturbed==True:
        kpert = 0;
        print("The perturbed neuron is kneu=%d.\n\n" % kpert);
        net.v[kpert] = -58.592274*mV;
        net.m[kpert] = .470296;
        net.h[kpert] = .035708;
        net.n[kpert] = .910829;
        net.a[kpert] = .452683;



    # SYNAPSES

    # If chemical synapses are activated
    if act_syns==True:

        # Definition of synapses
        print("Setting synapses...")

        # Include STD at the chem syns
        if std == True:

            synapses_eqs = '''
            # Fraction of synaptic neurotransmitter resources available:
            dxs/dt = (1 - xs)/tau_d : 1 (event-driven)
            gmax : siemens
            U_SE : 1
            tau_d : second
            '''
            synapses_action = '''
            rs = U_SE * xs
            xs -= rs
            g_fall_post+=gmax*rs
            g_rise_post+=gmax*rs
            '''

        # If STD is not included
        else: synapses_eqs = "gmax : siemens"; synapses_action = "g_fall_post+=gmax;g_rise_post+=gmax";

        # Print the used synaptic equations
        print("Synapse equations are:\n %s\n\n" % synapses_eqs )
        print("Synapse actions are:\n %s\n\n" % synapses_action )

        # Define the equation population
        S = Synapses(net,net,model=synapses_eqs,on_pre=synapses_action,method=method)

        # Connectivity matrix can be read from file. Otherwise Erdos-Renyi (connectivity="ER") or fixed number of presynaptic partners, i.e. in-degree (connectivity="FID")
        if read_syns==True: S.connect(i=syns[0], j=syns[1])
        elif connectivity=="ER": S.connect(p=epsilon, condition="i!=j")
        elif connectivity=="FID":
            C = int(np.round(epsilon*N));
            for k in range(N): pres = [ kp for kp in range(N) if kp!=k ]; shuffle(pres); S.connect(i=pres[:C],j=k);
        else: print("ERROR: Unknown connectivity type"); return "ERROR: Unknown connectivity type";


        if read_syns==True and read_delays==True: S.delay = delays;
        else: S.delay = (dmin+(dmax-dmin)*np.random.rand(len(S)))*ms;
        if std:
            S.xs = 1.; # Initialize auxiliary variable for STD
            S.U_SE = USEm; S.tau_d = taudm;
        if dist_syns==True:
            if read_syns==False or read_gms==False: gms = np.random.lognormal(gmmean,gmstd,size=len(S))*nS; # Set them to lognormal if not read
            if std==True:
                if scale_gm==False: S.gmax = gms/S.U_SE; # In order to see the effects of upscaling gms preserving their distribution 
                else: S.gmax = gm_scale*gms/S.U_SE; # If STD is added, gm is divided by U_SE. Then results with short tau_d and without STD agree.
            else:
                if scale_gm==False: S.gmax = gms;
                else: S.gmax = gm_scale*gms;
        else:
            if std: S.gmax = gms/S.U_SE;
            else: S.gmax = gms;

    # GAP JUNCTIONS

    if gjs==True:

        gap_eqs = """
                 g_gap : siemens # gap junction conductance
                 Igap_post = g_gap * (v_pre - v_post) : ampere (summed)
                 """;

        print("Gap-junction equations are:\n %s\n\n" % gap_eqs );

        Sgj = Synapses(net, net, gap_eqs, method=method);

        if read_gjs==True: Sgj.connect(i=synsgj[0], j=synsgj[1])
        else: Sgj.connect(p=pgj,condition="i!=j"); # WRONG! THESE ARE NOT BIDIRECTIONAL, SYMMETRIC AND NON-RECTIFYING

        Sgj.g_gap = ggaps; # They can be read from file or set on function call. Can be an array of length synsgj or a single float value









    # MONITORS
    print("Setting monitors...")
    rMon = PopulationRateMonitor(net); # Population rate monitor to return population-averaged (and possibly time-smoothed) firing rates
    spMon = SpikeMonitor(net, record=True) # Monitor for spike times

    # Monitor for state variables. If return_state==True. If gap junctions are not activated, return membrane potential, synaptic current and external current.
    # If they are activated return also net gap junction current onto the neuron.
    # Record and return state variables for the first Nrec neurons
    if return_state==True:
        state_vars = ["v","Iext"];#["v","m","n","h","a","Iext"]#
        if gjs==True: state_vars.append("Igap");
        if act_syns==True: state_vars.append("Isyn");
        stMon = StateMonitor(net, state_vars, record=range(Nrec), dt=dt_rec*b2.ms)









    # RUN
    print("Running the simulation...")
    b2.run(sim_time*b2.ms)


    # PROCESS
    rates = [rMon.t/b2.ms,rMon.smooth_rate(window=windowType,width=window_width)/b2.Hz] # Save and return time-smoothed population-averaged firing rates
    # Return also time- and population-averaged spike-count rate. Averaged across the whole simulation.
    nu = np.sum(rates[1][-40000:])/len(rates[1][-40000:]); print( "\nnu = %.2lf Hz" % nu );

    # END
    print("Done (took %.2lf mins)" % ((time()-st)/60.) )
    if return_state==True: return rates, spMon, stMon
    elif return_syns==True: return rates, spMon, S
    else: return rates, spMon





