import numpy as np

phminsYSTD = np.loadtxt("phmins_Esm75_YSTD.dat")
phmaxsYSTD = np.loadtxt("phmaxs_Esm75_YSTD.dat")

print(np.mean(phminsYSTD),np.std(phminsYSTD))
print(np.mean(phmaxsYSTD),np.std(phmaxsYSTD))
