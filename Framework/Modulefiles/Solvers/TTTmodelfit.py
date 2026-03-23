from Framework.HelpFile import getTTTdb
import numpy as np

def JMAKfit(rawdata, minput):
    try:
        start = rawdata["start"]
        half = rawdata["half"]
        finish = rawdata["finish"]
    except:
        print("TTTdata for doesn't exist")
        return None, None, None
    tau = np.array([])
    n = np.array([])
    Tlist = np.array([])
    for x in start[0]:
        # Explain!
        i = np.where(start[0] == x)[0][0]
        j = np.where(half[0] == x)[0][0]
        k = np.where(finish[0] == x)[0][0]
        if i.size == 0 or j.size == 0 or k.size == 0:
            pass
        elif start[1][j] == finish[1][i]: # ignore data point with the same values
            pass
        else:
            tmpn = 3.0


            tmptau = start[1][i] / ((-np.log(0.98)) ** (1 / tmpn))
            n = np.append(n, tmpn)
            tau = np.append(tau, tmptau)
            Tlist = np.append(Tlist, x)


    return Tlist, tau, n

def KMfit(raw_data, minput): # Koistinen marburger fitting process
    start = raw_data["start"][0][0]
    half = raw_data["half"][0][0]
    finish = raw_data["finish"][0][0]
    # 0.02 = 1- exp(-beta * (Ms - start))
    # 0.98 = 1- exp(-beta * (Ms - finish))
    if np.isnan(start):
        Ms = np.nan
        beta = np.nan
    elif np.isnan(half):
        Ms = start
        beta = 0.01
    elif np.isnan(finish):
        Ms = (np.log(0.98) * half - np.log(0.5) * start) / (np.log(0.98) - np.log(0.5))
        beta = -np.log(0.98) / (Ms - start)
    else:
        Ms = (np.log(0.98) * finish - np.log(0.02) * start) / (np.log(0.98) - np.log(0.02))
        beta = -np.log(0.98)/(Ms - start)
    return Ms, beta