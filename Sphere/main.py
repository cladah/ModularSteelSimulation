import time
import threading
import meshio
import vtk
import pyvista
import numpy as np
from vtkmodules import vtkIOXdmf2

from CGUImodule import MainApp
from Datastream_file import createdatastream, createdatastreamcache, removedatastreamcache, savedatastream
from HelpFile import read_input, setupSimulation, createinputcache
import customtkinter as ctk
from Modulefiles.Meshing_file import Meshingmodule
from Modulefiles.Carbonitriding_file import Carbonitridingmodule
from Modulefiles.TTTdiagram_file import TTTdiagrammodule
from Modulefiles.Transformationmodel_file import Transformationmodelmodule
from Modulefiles.Quenching_file import Quenchingmodule

import Testfile

def GUI():
    ctk.set_appearance_mode("dark")
    app = MainApp()
    app.mainloop()
    removedatastreamcache()


def modelling():
    setupSimulation()

    data = read_input()

    createdatastreamcache(data["Datastream"]["Cachedirect"])
    # resetdatastream()
    createinputcache()

    modules = list()
    modules.append(Meshingmodule())
    modules.append(Carbonitridingmodule())
    modules.append(TTTdiagrammodule())
    modules.append(Transformationmodelmodule())
    modules.append(Quenchingmodule())

    for currentmodule in modules:
        if currentmodule.modulename() != "Meshing":
            tid = threading.Thread(target=run_single_module, args=(currentmodule,))
            tid.start()
            progressmonitor(tid, currentmodule)  # Making sure thread is done
        else:
            run_single_module(currentmodule)
    removedatastreamcache()
    savedatastream(data["Datastream"]["Savedirect"])

def run_single_module(module):
    module.run()

def progressmonitor(tid, module):
    if tid.is_alive():
        time.sleep(1)
        progressmonitor(tid, module)

def xmdftesting():

    # modules = list()
    # modules.append(Meshingmodule())
    # for currentmodule in modules:
    #     if currentmodule.modulename() != "Meshing":
    #         tid = threading.Thread(target=run_single_module, args=(currentmodule,))
    #         tid.start()
    #         progressmonitor(tid, currentmodule)  # Making sure thread is done
    #     else:
    #         run_single_module(currentmodule)
    vtk_array = vtk.vtkDoubleArray()
    vtk_array.SetArray(np.linspace(0, 1, 2029), 2029, 0)
    vtk_array.SetObjectName("Testdata")

    print("time " + str(vtk_array.GetMTime()))
    new_grid = vtk.vtkDataObject()
    vtk.vtkXdmfDataArray()

    writer2 = vtk.vtkXMLUnstructuredGridWriter()
    reader = vtk.vtkXdmfReader()
    reader.SetFileName("Datastream.xdmf")
    reader.Update(0)
    info = reader.GetOutputInformation(0)
    print(info)
    #info = reader.GetOutputInformation(0)
    #print(reader.GetPointArrayName(0))
    #reader = xdmf.XdmfReader.New()
    #print(info)
    grid = reader.GetOutput()
    # grid.GetBlock(1).GetPointData().AddArray(vtk_array)
    # print(grid)
    nrnodes = grid.GetNumberOfPoints()

    print(grid.GetBlock(1).GetUpdateTime())
    data = grid.GetBlock(1).GetPointData().GetNumberOfArrays()
    # data = [grid.GetPointData().GetArray("Composition/C").GetValue(i) for i in range(nrnodes)]
    print(nrnodes)
    print(data)
    key = vtk.vtkInformationIntegerKey.MakeKey("time", "time")
    timestep = vtk.vtkInformation()
    timestep.SetObjectName("Time")
    timestep.Set(key, 1)
    testobj = vtk.vtkFieldData()
    testobj.SetObjectName("Time2")
    testann = vtk.vtkAnnotation()
    testann.SetInformation(timestep)
    # testobj.SetTuple(np.array([1,2,3]))

    #testobj.SetTuple(1)
    grid.GetBlock(1).SetInformation(timestep)
    grid.GetBlock(1).SetFieldData(testobj)
    writer = vtk.vtkXdmfWriter()
    writer.SetFileName("Datastream.xdmf")
    writer.SetInputData(grid)
    writer.Write()
    # Write the updated XDMF file



    #reader.UpdateTimeStep(0)


if __name__ == "__main__":
    #Testfile.read_data_from_xdmf("Resultfiles/230124_2.xdmf", 0)
    #Testfile.read_data_from_xdmf("Resultfiles/Test.xdmf", 0)
    #Testfile.add_data_to_xdmf("Resultfiles/Datastream.xdmf", [], 0)


    modelling()
    #xmdftesting()




    #print(np.max(data))
    # data = read_input()
    # createdatastreamcache(data["Datastream"]["Cachedirect"])
    # Meshingmodule().run()
    # Carbonitridingmodule().run()
    # TTTdiagrammodule().run()
    # removedatastreamcache()
    # savedatastream(data["Datastream"]["Savedirect"])


    # Bug in meshio TimeSeriesWriter
    # self.h5_filename = self.filename.stem + ".h5"
    # self.h5_filename = self.filename.with_suffix(".h5")