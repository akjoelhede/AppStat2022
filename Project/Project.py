#%%
import pandas as pd  
import numpy as np                                     # Matlab like syntax for linear algebra and functions
import matplotlib as mpl
import matplotlib.pyplot as plt                        # Plots and figures like you know them from Matlab
from iminuit import Minuit                             # The actual fitting tool, better than scipy's
import sys     
from scipy import stats
from scipy.stats import binom, poisson, norm 

# %%

pend1 = pd.read_csv("data_Alice_pendulum14m_25measurements.txt",sep='\t', names = ["n", "m"])
pend2 = pd.read_csv("data_Bob_pendulum14m_25measurements.txt",sep='\t', names = ["n", "m"])
inc1 = pd.read_csv("data_NormDir_MedBall1.txt",sep='\t', names = ["n", "m"])
inc2 = pd.read_csv("data_RevDir_MedBall1.txt",sep='\t', names = ["n", "m"])
# %%

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(pend1["n"], pend1["m"], 'x', color = 'red')
axs[0, 0].set_title("Alice Pendulum")
axs[0, 0].set_xlabel("Experiments, [n]")
axs[0, 0].set_ylabel("Time, [t]")
axs[1, 0].plot(pend2["n"], pend2["m"], 'x', color = 'black')
axs[1, 0].set_title("Bob Pendulum")
axs[1, 0].set_xlabel("Experiments, [n]")
axs[1, 0].set_ylabel("Time, [t]")
axs[1, 0].sharex(axs[0, 0])
axs[0, 1].plot(inc1["n"], inc1["m"])
axs[0, 1].set_title("Normdir Incline")
axs[0, 1].set_xlabel("Time, [t]")
axs[0, 1].set_ylabel("Voltage dif, [dV]")
axs[1, 1].plot(inc2["n"], inc2["m"])
axs[1, 1].set_title("Revdir Incline")
axs[1, 1].set_xlabel("Time, [t]")
axs[1, 1].set_ylabel("Voltage dif, [dV]")
fig.tight_layout()


# %%

def sigma_pendulum(T,L,sigL,sigT):
	return np.sqrt((2*np.pi/T)**2*sigL - L*8*np.pi**2/T**3 * sigT)

def sigma_incline(a, theta, D, d, sigA, sigTheta, sigD, sigd):
	np.csc(theta)(1 + 2/5 * D**2/(D**2-d**2)) * sigA -a (1 + 2/5 * D**2/(D**2-d**2))* np.cot(theta)*np.csc(theta)*sigTheta - a*np.csc(theta)*5*d**2*D/(d**2-D**2)**2 *sigD + a*np.csc(theta)*5*d*D**2/(d**2-D**2)**2 * sigd
# %%
