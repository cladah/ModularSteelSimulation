from HelpFile import *
def JMAKfit(composition, phase):
    import numpy as np
    TTTdata = getTTTdata(composition,"TTTdata")
    print(TTTdata[phase]["start"])
    input("Pause")
    data1 = [start, Tsteps]
    data2 = [half, Tsteps]
    data3 = [finish, Tsteps]
    import numpy as np
    tau = np.array([])
    n = np.array([])
    Tlist = np.array([])
    for x in data1[1]:
        i = np.where(data3[1] == x)[0]
        k = np.where(data2[1] == x)[0]
        j = np.where(data1[1] == x)[0]
        if i.size==0 or j.size==0 or k.size==0:
            pass
        elif data1[0][j[0]] == data3[0][i[0]]: # Take away data point with the same values
            pass
        else:
            # add if to take 98% or 50% lines as references.
            i = i[0]
            j = j[0]
            tmpn = np.log(np.log(0.98)/np.log(0.02))/np.log(data1[0][j]/data3[0][i])
            tmptau = - data1[0][j]/(-np.log(0.98))**(1/tmpn)
            n = np.append(n, tmpn)
            tau = np.append(tau, tmptau)
            Tlist = np.append(Tlist, x)
    return Tlist, n, tau

def KMfit(phasename,filename): # Koistinen marburger fitting process
    import numpy as np
    start = readresultfile(filename, phasename + "/start")
    half = readresultfile(filename, phasename + "/half")
    finish = readresultfile(filename, phasename + "/finish")

    # 0.02 = 1- exp(-beta * (Ms - start))
    # 0.98 = 1- exp(-beta * (Ms - finish))
    #np.log(0.98)*(Ms-finish) = np.log(0.02)*(Ms-start)
    Ms = (np.log(0.98)*finish-np.log(0.02)*start)/(np.log(0.98)-np.log(0.02))
    beta = -np.log(0.98)/(Ms - start)

    return Ms,beta