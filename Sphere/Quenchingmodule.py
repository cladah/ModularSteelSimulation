from Solvers.ComsolSolver import setupComsol
from HelpFile import read_input, checkinput, adjustinputcache
import os
def runquenchingmodule():
    if checkinput('Quenching'):
        print('Using old quenching simulation')
        return
    print('Quenching module')
    data = read_input()
    if data['Programs']['FEM'] == 'FCSx':
        print('FeniCSx not implemented')
        #rundocker()
    elif data['Programs']['FEM'] == 'Comsol':
        print('Using COMSOL for FEM calculation')
        if not os.path.exists('Resultfiles/Comsolmodel.mph'):
            setupComsol()
        else:
            print("comsolmodel exists")
            pass
    else:
        raise KeyError('Program not implemented')
    adjustinputcache('Mesh')