# -*- coding: utf-8 -*-
import numpy as np
from scipy.signal import cwt, morlet2, butter, sosfiltfilt
import matplotlib.pyplot as plt




def compPWT(signal,fs=np.arange(50.,450.1,3.),dt=1.e-4,word=5,tmin=.0):
  """
  Computes the Wavelet Power Spectra (modulus of complex wavelet transform) for a set of signals
  Input:
   signal (to be analyzed)
   fs (frequencies to compute cwt)
   dt
   word (wavelet order, default=5)
   nlevels (number of levels at contour plot)
   tmin (plot wavelet spectra from tmin to end of signal. in seconds)
  """

  ktmin = int(tmin/dt);

  # Perform Wavelet Analysis
  samp_rate = 1./dt; # Sampling rate is inverse of time step dt
  ss = word*samp_rate/(2.*np.pi*fs); # Compute corresponding wavelet scalings or widths (s).
  pwt = np.power(np.abs(cwt(signal[ktmin:]-np.mean(signal[ktmin:]),morlet2,widths=ss)),2.) # Compute wavelet power, i.e. squared modulus, using a Morlet wavelet

  return pwt

def compPhaseWT(signal,fs=np.arange(50.,450.1,3.),dt=1.e-4,word=5,tmin=.0):
  """
  Computes the Wavelet Power Phase (modulus of complex wavelet transform) for a set of signals
  Input:
   signal (to be analyzed)
   fs (frequencies to compute cwt)
   dt
   word (wavelet order, default=5)
   nlevels (number of levels at contour plot)
   tmin (plot wavelet spectra from tmin to end of signal. in seconds)
  """

  ktmin = int(tmin/dt);

  # Perform Wavelet Analysis
  samp_rate = 1./dt; # Sampling rate is inverse of time step dt
  ss = word*samp_rate/(2.*np.pi*fs); # Compute corresponding wavelet scalings or widths (s).
  #pwt = np.angle(cwt(signal[ktmin:]-np.mean(signal[ktmin:]),morlet2,widths=ss)) # Compute wavelet phase, using a Morlet wavelet
  pwt = np.angle(cwt(signal[ktmin:]-np.mean(signal[ktmin:]),morlet2,widths=ss))
  return pwt

def plot1(cwtPow,max_power=None,fs=np.arange(50.,450.1,3.),tmin=.0,dt=.1,nlevels=30,Ncycs=1,flims=[],xlims=[],colorbar=True,cmap="hot"):
    """
    plot the contour plot for the wavelet scalogram for one single simulation
    """
    ktmin = int(np.round(tmin/dt));
    Nt_pr_cyc = int(np.round(125./dt));
    if max_power==None: max_power = np.max(np.max(cwtPow));
    levels = np.linspace(0.,max_power+1.,nlevels);
    tmp = plt.contourf(np.linspace(-np.pi,-np.pi+2.*np.pi*Ncycs,int(np.round(Ncycs*Nt_pr_cyc))),fs,np.abs(cwtPow[:,ktmin:]),levels=levels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = plt.colorbar();

    if len(flims)==0: tmp = plt.ylim(fs[0],fs[-1]);
    else: tmp = plt.ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = plt.xlim(-np.pi,-np.pi+2.*np.pi*Ncycs);
    else: tmp = plt.xlim(xlims[0],xlims[1]);

def plotPhase(cwtPow,max_power=None,fs=np.arange(50.,450.1,3.),tmin=.0,dt=.1,nlevels=30,Ncycs=1,flims=[],xlims=[],colorbar=True,cmap="hot"):
    """
    plot the contour plot for the wavelet scalogram for one single simulation
    """
    ktmin = int(np.round(tmin/dt));
    Nt_pr_cyc = int(np.round(125./dt));
    #if max_power==None: max_power = np.max(np.max(cwtPow));
    levels = 30#np.linspace(0.,max_power+1.,nlevels);
    tmp = plt.contourf(np.linspace(-np.pi,-np.pi+2.*np.pi*Ncycs,int(np.round(Ncycs*Nt_pr_cyc))),fs,cwtPow[:,ktmin:],levels=levels,cmap=cmap,rasterized=True);

    if colorbar==True: tmp = plt.colorbar();

    if len(flims)==0: tmp = plt.ylim(fs[0],fs[-1]);
    else: tmp = plt.ylim(flims[0],flims[1]);

    if len(xlims)==0: tmp = plt.xlim(-np.pi,-np.pi+2.*np.pi*Ncycs);
    else: tmp = plt.xlim(xlims[0],xlims[1]);


