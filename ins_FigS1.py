import numpy as np
import matplotlib.pyplot as plt

ModelTimeConstant = np.loadtxt("TxtFigS1/taumsRomanUniform.dat")
ExpTimeConstant = np.loadtxt("TxtFigS1/taums_unclasspv.dat")

ModelInputResistance = 1.e+3/np.loadtxt("TxtFigS1/gLsRoman.dat")
ExpInputResistance = np.loadtxt("TxtFigS1/Rins_unclasspv.dat")

ModelRMP = np.loadtxt("TxtFigS1/RMPModel.txt")
ExpRMP = np.loadtxt("TxtFigS1/vrests_unclasspv.dat")

ModelCutoff = np.loadtxt("TxtFigS1/CutoffFreqModel.txt")
ExpCutoff = np.loadtxt("TxtFigS1/CutOffFreqExp.dat")

ModelRheobase = np.loadtxt("TxtFigS1/RheobaseModel.txt")
ExpRheobase = np.loadtxt("TxtFigS1/RheobaseExp.dat")

plt.rcParams["font.family"] = "Arial"
limbinsTimeExp = 61#np.linspace(0,10,26)
limbinsTimeModel = 10
limbinsRin = 26#np.linspace(30,160,26)
limbinsRMPExp = 26#np.linspace(-90,-60,26)
limbinsRMPModel = 12#np.linspace(-90,-60,41)
limbinsCutoff = 20#np.linspace(0,200,26)
limbinsRheobase = 15#np.linspace(0,10,26)

width = 5.2 # in inches
height = 7.# in inches
tmp = plt.figure(figsize=[width,height]);

ax = plt.subplot(5,2,1)
plt.title("Experimental Data")
plt.hist(ExpTimeConstant,limbinsTimeExp)
plt.ylim((0,25)); plt.xlim((1,11));
plt.xlabel("Time Constant (ms)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,2)
plt.title("Model Neurons")
plt.hist(ModelTimeConstant,limbinsTimeModel)
plt.xlim((1,11)); plt.ylim((0,25))
plt.xlabel("Time Constant (ms)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,3)
plt.hist(ExpInputResistance,limbinsRin)
plt.xlim((30,160)); plt.ylim((0,20))
plt.xlabel("Input Resistance (M$\Omega$)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,4)
plt.hist(ModelInputResistance,limbinsRin)
plt.xlim((30,160)); plt.ylim((0,20))
plt.xlabel("Input Resistance (M$\Omega$)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,5)
plt.hist(ExpRMP,limbinsRMPExp)
plt.xlim((-90,-50)); plt.ylim((0,25))
plt.xlabel("Resting Potential (mV)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,6)
plt.hist(ModelRMP,limbinsRMPModel)
plt.xlim((-90,-50)); plt.ylim((0,25))
plt.xlabel("Resting Potential (mV)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,7)
plt.hist(ExpCutoff,limbinsCutoff)
plt.xlim((0,200)); plt.ylim((0,30))
plt.xlabel("Cutoff Frequency (Hz)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,8)
plt.hist(ModelCutoff,limbinsCutoff)
plt.xlim((0,200)); plt.ylim((0,30))
plt.xlabel("Cutoff Frequency (Hz)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,9)
plt.hist(ExpRheobase,limbinsRheobase)
plt.xlim((0,600)); plt.ylim((0,20))
plt.xlabel("Rheobase (pA)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

ax = plt.subplot(5,2,10)
plt.hist(ModelRheobase,limbinsRheobase)
plt.xlim((0,600)); plt.ylim((0,20))
plt.xlabel("Rheobase (pA)")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

plt.tight_layout()
plt.savefig("Figures/FigS1.eps")
plt.savefig("Figures/FigS1.png")

