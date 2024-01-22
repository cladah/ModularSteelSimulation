from Solvers.ComsolSolver import setupComsol, runComsol, adjustComsol
from HelpFile import read_input, checkinput, adjustinputcache
import os
def runquenchingmodule(parent):
    if checkinput('Quenching'):
        print('Using old quenching simulation')
        return
    print('\nQuenching module')
    data = read_input()

    if data['Programs']['FEM'] == 'FCSx':
        print('FeniCSx not implemented')
        #rundocker()
    elif data['Programs']['FEM'] == 'Comsol':
        print('Using COMSOL for FEM calculation')
        runComsol(parent)
        return
        if not os.path.exists('Resultfiles/Comsolmodel.mph'):
            setupComsol()
        else:
            print("Using previous setup model")
            runComsol()
    else:
        raise KeyError('Program not implemented')
    adjustinputcache('Mesh')