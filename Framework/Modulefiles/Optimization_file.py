from .ModuleStructure_file import CalcModule
import numpy as np
from Framework.Datastream_file import readdatastreamcache, readdatastream, adjustdatastream, getaxisvalues
from Framework.HelpFile import read_input
from Meshing_file import Meshingmodule
from Carbonitriding_file import Carbonitridingmodule, Carbonizationmodule
from TTTdiagram_file import TTTdiagrammodule
from Transformationmodel_file import Transformationmodelmodule
from Quenching_file import Quenchingmodule
from ..Datastream_file import createdatastreamcache, removedatastreamcache, savedatastream
from ..HelpFile import read_input, setupSimulation, createinputcache, change_input, reset_input, analyseTTTdatabase, get_plotlbls


class Optimizationmodule(CalcModule):

    def __init__(self):
        super().__init__("Opitimization check")

    def run(self):
        datanames = ["Martensite", "Composition/C", "HV", "Strain"]

        error_est = [0.01, 0.001,]
        for dataname in datanames:
            ssd = self.run_ssd(dataname)
            if ssd < 1e-6:
                pass
            else:
                raise False
    def run_ssd(self, dataname):
        cal_data = getaxisvalues(dataname, time=-1)
        ssd = np.sum((np.array(cal_data, dtype=np.float32) - np.array(meas_data, dtype=np.float32)) ** 2)
        ssd_norm = ssd/np.mean(cal_data)
        return ssd_norm

class Optimization_loop():
    def __init__(self, exp_data):
        self.exp_data = exp_data
        self.nrofloops = 0
    def run_cycle(self):
        for i in range(20):
            self.nrofloops = self.nrofloops + 1
            # Using the previous loop as input for cached data
            change_input("Datastream", "Cachedirect", str(i-1) + ".xdfm")

            self.run_one_loop()
            self.calculate_res_grad()
            # Saving data as new xdmf file
            savedatastream(str(i)+".xdmf")

    def run_one_loop(self):
        modules = self.getmodules()
        for module in modules:
            module.run()
        datanames = ["Martensite", "Composition/C"]
        for dataname in datanames:
            datafile = "Datastream.xdmf" # Needs to be adjusted depending on the operating stsyem that is in use.

            ssd = self.check_cal_exp(datafile, dataname)


            # The martensite fraction is the most affected by the cooling. Optimization should be ran against the martensite Ms tem
            # This assumes the cooling is quick enough in the center so no bainite structure is forming.
            # Don't forget to add the cementite fraction that is spit out of Thermocalc
            # Ask how cementite affects the neutron imaging. Do I need to take this into account in the analysis?
            #


    def calculate_res_grad(self):
        for i in range(self.nrofloops):
            dataname = "Martensite"
            ssd = self.check_cal_exp(str(i) + ".xdmf", dataname)


    def getmodules(self):
        modules = list()
        modules.append(Meshingmodule())
        modules.append(Carbonizationmodule())
        modules.append(TTTdiagrammodule())
        modules.append(Transformationmodelmodule())
        modules.append(Quenchingmodule())
        return modules
    def check_cal_exp(self, datafile, dataname):
        cal_data = getaxisvalues(dataname, time=-1)
        ssd = np.sum((np.array(cal_data, dtype=np.float32) - np.array(self.exp_data, dtype=np.float32)) ** 2)
        ssd_norm = ssd / np.mean(cal_data)
        return ssd_norm