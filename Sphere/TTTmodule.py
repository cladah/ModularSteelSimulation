from Sphere.Solvers.Thermocalc import *
from Sphere.Solvers.TTTmodelfit import *
from HelpFile import *

def runTTTmodule():

    from Solvers.Thermocalc import calculateCCT
    import h5py
    import matplotlib.pyplot as plt
    import numpy as np
    from HelpFile import read_input, checkinput
    if checkinput('TTT'):
        print('Using precalculated TTT simulation')
        return
    print('CCT module')
    data = read_input()
    with h5py.File("Resultfiles/Carbonitriding.hdf5", "r") as f:
        x = np.array(f.get("CNcurves/Position"))
        a = dict()
        for element in data['Material']['Composition'].keys():
            #plt.plot(x, 100 * np.array(f.get("CNcurves/"+element)))
            a[element] = 100 * np.array(f.get("CNcurves/Elements/" + element)) # Composition curves at all points along x

    composition = data['Material']['Composition']
    #composition['C'] = a['C'][-1]
    #composition['N'] = a['N'][-1]
    Tsteps = np.linspace(300,1000,10)

    start, half, finish = calculateBainite(Tsteps, composition)

    with h5py.File("Resultfiles/TTT.hdf5", "w") as f:
        f.create_dataset("Bainite/Tsteps", data=Tsteps)
        f.create_dataset("Bainite/start", data=start)
        f.create_dataset("Bainite/half", data=half)
        f.create_dataset("Bainite/finish", data=finish)

    return
    #saveresult("TTT/Surface/Bainite/Start", start)
    #saveresult("TTT/Surface/Bainite/Half", half)
    #saveresult("TTT/Surface/Bainite/Finish", finish)

    #start, half, finish = calculatePerlite()
    #saveresult("TTT/Surface/Perlite/Start", start)
    #saveresult("TTT/Surface/Perlite/Half", half)
    #saveresult("TTT/Surface/Perlite/Finish", finish)

    #calculateCCT(composition)
    #Ccurve = np.array(f.get("CNcurves/C"))
    #Ncurve = np.array(f.get("CNcurves/N"))
    #plt.plot(Ccurve)
    #plt.plot(Ncurve)
    #plt.legend(data['Material']['Composition'].keys(), loc='upper left')
    #plt.show()
    #x = [239.7, 126.4, 73.5, 53.0, 55.6, 67.6, 68.1, 90.6, 185.0, 838.0]
    #y = [780, 800, 820, 840, 860, 880, 900, 920, 940, 960]
    #Perlite = modelfitting('JMAK', x, y)
    #print(Perlite)
    #T_new = np.linspace(700,1000,50)
    #t_new = Perlite(T_new).T

    #plt.plot(x, y)
    #plt.plot(t_new, T_new)
    #plt.show()
    #Martensite = modelfitting(data['Material']['Martensite']['model'], [1, 1], 2)
    #Perlite = modelfitting(data['Material']['Perlite']['model'], [1, 1], [1, 1])
    #Bainite = modelfitting(data['Material']['Bainite']['model'], [1, 1], [1, 1])

    with h5py.File("Resultfiles/ThermoResult.hdf5", "r+") as f:
        try:
            f.create_dataset('JMAK/Tau', data=1)
            f.create_dataset('JMAK/n', data=1)
        except ValueError:
            del f['JMAK/Tau']
            del f['JMAK/n']
            f.create_dataset('JMAK/Tau', data=1)
            f.create_dataset('JMAK/n', data=1)
def runTTTfitmodule():
    import h5py
    with h5py.File("Resultfiles/TTT.hdf5", "r") as f:
        Tsteps = np.array(f.get("Bainite/Tsteps"))
        start = np.array(f.get("Bainite/start"))
        half = np.array(f.get("Bainite/half"))
        finish = np.array(f.get("Bainite/finish"))
    b_Tlist, b_n, b_tau = JMAKfit([start,Tsteps], [half,Tsteps], [finish,Tsteps])
    saveresult("TTT.hdf5","Bainite/JMAK/T",b_Tlist)
    saveresult("TTT.hdf5", "Bainite/JMAK/n",b_n)
    saveresult("TTT.hdf5", "Bainite/JMAK/tau",b_tau)
    #with h5py.File("Resultfiles/TTT.hdf5", "r+") as f:
    #    f.create_dataset('Bainite/JMAK/Tau', data=b_Tlist)
    #   f.create_dataset('Bainite/JMAK/Tau', data=b_n)
    #    f.create_dataset('Bainite/JMAK/Tau', data=b_tau)