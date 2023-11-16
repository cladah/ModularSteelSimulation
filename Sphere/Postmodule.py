import matplotlib.pyplot as plt

from Postprocessing.dataextraction import *
from Postprocessing.dataplotting import *

def plotTTTdiagram():
    plotJMAK("Bainite", "TTT_center.hdf5")
    plotJMAK("Perlite", "TTT_center.hdf5")
    plotKM("Martensite", "TTT_center.hdf5")
    plt.xlim(1,1E12)
    plt.legend()
    plt.show()