from HelpFile import readdatastream
import numpy as np
def getaxisvalues(dataname, axis):

    node_y = readdatastream('nodes')[:, 1]
    indx = np.where(node_y == 0)
    y = readdatastream(dataname)[indx]
    x = readdatastream('nodes')[:, 0][indx]
    indx = np.argsort(x)
    x = np.array(x)[indx]
    y = np.array(y)[indx]
    return x, y