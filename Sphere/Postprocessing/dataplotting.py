import matplotlib.pyplot as plt
import numpy as np
from HelpFile import readresultfile
def plotJMAK(phasename):
    T = readresultfile("TTT.hdf5", phasename+"/JMAK/T")
    tau = readresultfile("TTT.hdf5", phasename + "/JMAK/tau")
    n = readresultfile("TTT.hdf5", phasename + "/JMAK/n")
    xpoints = -tau * (-np.log(0.98)) ** (1 / n)
    ypoints = T
    plt.plot(xpoints, ypoints)
    xpoints = -tau * (-np.log(0.02)) ** (1 / n)
    ypoints = T
    plt.plot(xpoints, ypoints)
    plt.xscale("log")