from Postprocessing.dataextraction import *
from Postprocessing.dataplotting import *

def plotTTTdiagram():
    plotJMAK("Bainite")
    plotJMAK("Perlite")
    #plotKM("Martensite")
    plt.show()