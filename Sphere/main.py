from HelpFile import *
from Meshmodule import createMesh
from Carbonitridingmodule import runcarbonitriding

def start():
    try:
        createinputcache()
    except:
        pass
    try:
        createdatastream()
    except:
        pass
    createMesh()
    runcarbonitriding()

if __name__ == "__main__":
    start()