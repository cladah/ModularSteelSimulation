import matplotlib.pyplot as plt

from Postprocessing.dataextraction import *
from Postprocessing.dataplotting import *

def plotTTTdiagram():
    #fig = plt.subplots(nrows=1, ncols=2, gridspec_kw=dict(hspace=0.3), figsize=(20, 15), sharex=True, sharey=True)
    plt.subplot(121)
    plotJMAK("Bainite", "TTT_center.hdf5")
    plotJMAK("Perlite", "TTT_center.hdf5")
    plotKM("Martensite", "TTT_center.hdf5")
    plt.xlim(1,1E12)
    plt.legend()
    plt.ylabel('Temperature [K]')
    plt.xlabel("Time [s]")

    plt.subplot(122)
    plotJMAK("Bainite", "TTT_surface.hdf5")
    plotJMAK("Perlite", "TTT_surface.hdf5")
    plotKM("Martensite", "TTT_surface.hdf5")
    plt.xlim(1, 1E12)
    plt.legend()
    plt.ylabel('Temperature [K]')
    plt.xlabel("Time [s]")
    plt.show()