from .ModuleStructure_file import CalcModule
import numpy as np
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream, getaxisvalues
from Framework.HelpFile import read_input
import pandas as pd
class Optimizationmodule(CalcModule):

    def __init__(self):
        super().__init__("Carburization")

    def run(self):
        dataname_Xi = "Martensite"
        dataname_C = "Composition/C"
        Xi_M = getaxisvalues(dataname_Xi, time=-1)
        wC = getaxisvalues(dataname_C)
        hardness_data = pd.read_csv("DataFiles/Hardness_measurement.csv")
        neutron_data = pd.read_csv("DataFiles/Neutron_Fit_Results.csv")

        measuredData = np.array(yvalues['int1'])
        calibrationData = np.array(yvalues['int0'])

        A = np.vstack([measuredData, np.ones(len(measuredData))]).T
        gain, offset = np.linalg.lstsq(A, calibrationData)[0]
