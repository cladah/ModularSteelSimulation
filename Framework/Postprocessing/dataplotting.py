import matplotlib.pyplot as plt
import numpy as np
from .HelpFile import readresultfile, getTTTdata
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
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
    indx = xpoints<1E12*np.logical_not(np.isnan(tau))
    T = T[indx]
    tau = tau[indx]

    print("adding point to " + phasename)
    n = n[indx]
    # Adding a point far to the right to avoid phasetransformation at high T

    tau = np.append(tau, tau[-1]*1E12)
    n = np.append(n, n[-1])
    T = np.append(T, 2*T[-1]-T[-2])


    testT = np.linspace(270, 1000, 74)
    # Polytest
    testtau = np.polyfit(T,tau,1)
    testn = np.polyfit(T, n, 1)
    taunew = np.poly1d(testtau)
    nnew = np.poly1d(testn)
    #tauinter = taunew(testT)
    #ninter = nnew(testT)

    # Spline test
    taufunc = interpolate.splrep(T, tau, s=0)
    nfunc = interpolate.splrep(T, n, s=0)
    tauinter = interpolate.splev(testT, taufunc, der=0)
    ninter = interpolate.splev(testT, nfunc, der=0)
    xpoints = -tauinter * (-np.log(0.50)) ** (1 / ninter)
    ypoints = testT
    #print(xpoints)
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

def showMeshPlot(nodes, elements):

    x = nodes[:, 0]
    y = nodes[:, 1]

    #https://stackoverflow.com/questions/49640311/matplotlib-unstructered-quadrilaterals-instead-of-triangles
    def quatplot(y, z, quatrangles, ax=None, **kwargs):

        if not ax: ax=plt.gca()
        yz = np.c_[y, z]
        verts= yz[quatrangles]
        pc = mpl.collections.PolyCollection(verts, **kwargs)
        ax.add_collection(pc)
        ax.autoscale()

    plt.figure()
    plt.gca().set_aspect('equal')

    quatplot(x, y, np.asarray(elements), color="crimson", facecolor="None")
    # if nodes:
    #     plt.plot(x, y, marker="o", ls="", color="crimson")

    plt.title('This is the plot for: quad')
    plt.xlabel('Y Axis')
    plt.ylabel('Z Axis')


    plt.show()