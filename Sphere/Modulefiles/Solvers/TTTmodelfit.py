from Sphere.HelpFile import getTTTdata
import numpy as np
def JMAKfit(composition, phase):
    TTTdata = getTTTdata(composition,"TTTdata")
    if TTTdata == None:
        raise KeyError("Can't find " + str(composition) + " in database")
    try:
        start = TTTdata[phase]["start"]
        half = TTTdata[phase]["half"]
        finish = TTTdata[phase]["finish"]
    except:
        print("TTTdata for " + str(composition) + " doesn't exist for " + phase + " transformation")
        return None, None, None
    tau = np.array([])
    n = np.array([])
    Tlist = np.array([])
    for x in start[0]:
        i = np.where(start[0] == x)[0][0]
        j = np.where(half[0] == x)[0][0]
        k = np.where(finish[0] == x)[0][0]
        if i.size==0 or j.size==0 or k.size==0:
            pass
        elif start[1][j] == finish[1][i]: # ignore data point with the same values
            pass
        else:
            # add if to take 98% or 50% lines as references.
            if np.isnan(finish[1][k]):
                tmpn = np.log(np.log(0.98) / np.log(0.50)) / np.log(start[1][i] / half[1][j])

                tmpn = 3.0

                tmptau = start[1][i] / ((-np.log(0.98)) ** (1 / tmpn))
            else:
                tmpn = np.log(np.log(0.98) / np.log(0.02)) / np.log(start[1][i] / finish[1][k])
                tmpn = 3.0

                tmptau = start[1][i] / ((-np.log(0.98)) ** (1 / tmpn))


            n = np.append(n, tmpn)
            tau = np.append(tau, tmptau)
            Tlist = np.append(Tlist, x)

    # print(phase)
    # print(n)

    return Tlist, tau, n

def KMfit(composition, phase): # Koistinen marburger fitting process
    import numpy as np
    TTTdata = getTTTdata(composition, "TTTdata")

    start = TTTdata[phase]["start"][0][0]
    half = TTTdata[phase]["half"][0][0]
    finish = TTTdata[phase]["finish"][0][0]
    # 0.02 = 1- exp(-beta * (Ms - start))
    # 0.98 = 1- exp(-beta * (Ms - finish))
    #np.log(0.98)*(Ms-finish) = np.log(0.02)*(Ms-start)
    Ms = (np.log(0.98)*finish-np.log(0.02)*start)/(np.log(0.98)-np.log(0.02))
    beta = -np.log(0.98)/(Ms - start)

    return Ms, beta