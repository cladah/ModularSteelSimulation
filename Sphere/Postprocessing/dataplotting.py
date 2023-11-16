import matplotlib.pyplot as plt
import numpy as np
from HelpFile import readresultfile
from scipy import interpolate
def plotJMAK(phasename,filename):
    T = readresultfile(filename, phasename+"/JMAK/T")
    tau = readresultfile(filename, phasename + "/JMAK/tau")
    n = readresultfile(filename, phasename + "/JMAK/n")
    xpoints = -tau * (-np.log(0.98)) ** (1 / n)
    ypoints = T
    plt.plot(xpoints, ypoints, label =phasename + ' start (2%)')
    xpoints = -tau * (-np.log(0.02)) ** (1 / n)
    ypoints = T
    plt.plot(xpoints, ypoints, label =phasename + ' finish (98%)')
    indx = np.logical_not(np.isnan(tau))
    T = T[indx]
    tau = tau[indx]
    n = n[indx]

    # Polytest
    testtau = np.polyfit(T,tau,5)
    testn = np.polyfit(T, n, 1)
    taunew = np.poly1d(testtau)
    nnew = np.poly1d(testn)
    #tauinter = taunew(T)
    #ninter = nnew(T)

    # Spline test
    taufunc = interpolate.splrep(T, tau, s=0)
    nfunc = interpolate.splrep(T, n, s=0)
    tauinter = interpolate.splev(T, taufunc, der=0)
    ninter = interpolate.splev(T, nfunc, der=0)
    xpoints = -tauinter * (-np.log(0.50)) ** (1 / ninter)
    ypoints = T

    plt.plot(xpoints, ypoints, label=phasename + ' test (50%)')

    plt.xscale("log")
def plotKM(phasename,filename):
    Ms = readresultfile(filename, phasename+"/KM/Ms")
    beta = readresultfile(filename, phasename + "/KM/beta")
    xpoints = [1,1E12]
    M02 = Ms + np.log(0.98)/beta
    M50 = Ms + np.log(0.5) / beta
    M98 = Ms + np.log(0.02) / beta
    plt.plot(xpoints, [M02,M02], label =phasename + ' start (2%)')
    plt.plot(xpoints, [M50, M50], label =phasename + ' half (50%)')
    plt.plot(xpoints, [M98, M98], label =phasename + ' start (98%)')
    plt.xscale("log")