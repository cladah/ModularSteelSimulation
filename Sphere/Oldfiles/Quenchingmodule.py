from Sphere.Modulefiles.Solvers.ComsolSolver import setupComsol, runComsol
from HelpFile import read_input, checkruncondition, adjustinputcache
import os
def runquenchingmodule(parent):
    if not checkruncondition('Quenching'):
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
        if not os.path.exists('../Resultfiles/Comsolmodel.mph'):
            setupComsol()
        else:
            print("Using previous setup model")
            runComsol()
    else:
        raise KeyError('Program not implemented')
    adjustinputcache('Mesh')