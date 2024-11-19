import docker
import os
import time

import subprocess
def rundocker():
    directory = os.getcwd()
    print(directory)
    dockervolume = directory + ':/root/shared'
    dockervolume = dockervolume.replace('\\', '/')

    if 1 == 1:
        client = docker.from_env()
        print('Creating docker container\n') # dokken92/dolfinx_mpc # dolfinx/dolfinx:v0.8.0 # ghcr.io/jorgensd/dolfinx_mpc:v0.8.1
        container = client.containers.run('dolfinx/dolfinx:v0.9.1', ["python3", "Modulefiles/Solvers/Fenicsx_test.py"],
                                          detach=True,
                                          auto_remove=True,
                                          #tty=True,
                                          #stdin_open=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          environment=['PYTHONPATH=/usr/local/lib/python3/dist-packages:/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:'],
                                          name='fenicscxcont')
        print('Running FeniCSx\n ')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print(" \n")
        print('FeniCSx calculation done!')
        print('Removing docker container')

    else:
        raise KeyError('Solver not implemented')
        client = docker.from_env()
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Solvers/CoupledSolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('Fenicsx calculation done!')


def rundocker_1D():
    directory = os.getcwd()
    print(directory)
    dockervolume = directory + ':/root/shared'
    dockervolume = dockervolume.replace('\\', '/')

    if 1 == 1:
        client = docker.from_env()
        print('Creating docker container') # dokken92/dolfinx_mpc # dolfinx/dolfinx:v0.8.0 # dolfinx/dolfinx:stable
        container = client.containers.run('ghcr.io/jorgensd/dolfinx_mpc:v0.8.1', ["python3", "Modulefiles/Solvers/FenicsxSolvers.py"],
                                          detach=True,
                                          auto_remove=True,
                                          #tty=True,
                                          #stdin_open=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          environment=['PYTHONPATH=/usr/local/lib/python3/dist-packages:/usr/local/dolfinx-real/lib/python3.10/dist-packages:/usr/local/lib:'],
                                          name='fenicscxcont')
        print('Running FeniCSx\n')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('\nFeniCSx calculation done!')
        print('Removing docker container')

    else:
        raise KeyError('Solver not implemented')
        client = docker.from_env()
        container = client.containers.run('dolfinx/dolfinx:stable', ["python3", "Solvers/CoupledSolver.py"],
                                          detach=True,
                                          auto_remove=True,
                                          volumes=[dockervolume],
                                          working_dir='/root/shared',
                                          name='fenicscxcont')
        for log in container.logs(stream=True, stdout=True, stderr=True):
            print(log)
        print('Fenicsx calculation done!')