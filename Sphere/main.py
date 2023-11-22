import os

import numpy as np
import h5py
import matplotlib.pyplot as plt

from HelpFile import *
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitridingmodule
from TTTmodule import runTTTmodule, runTTTfitmodule
from Postmodule import *
from Quenchingmodule import runquenchingmodule



def start():
    if not os.path.exists('Cachefiles/InputCache.json'):
        createinputcache()
    createMesh()
    runcarbonitridingmodule()
    runTTTmodule()
    runTTTfitmodule()
    #runquenchingmodule()
    #plotTTTdiagram()
    #x,y = getaxisvalues('C', 0)
    #plt.plot(x,y)
    #plt.show()
    #plt.close()
if __name__ == "__main__":
    start()