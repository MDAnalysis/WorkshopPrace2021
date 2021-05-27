#!/usr/bin/env python
import MDAnalysis as mda
import numpy as np
from mpi4py import MPI
import sys
import platform

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_workers = comm.Get_size()

# We print some simple reporting messages so that the parallelism can be followed in the log.
print(f'Rank {rank} reporting for duty on node {platform.node()}, in a pool of {n_workers} workers.')

basedir = '/project/jhlsrf005/JHL_data/Day3-Session1-Practical/'
u = mda.Universe(basedir + 'md.tpr',
                 basedir + 'md.xtc',
                 basedir + 'md.part0002.xtc',
                 basedir + 'md.part0003.xtc',
                 basedir + 'md.part0004.xtc',
                 basedir + 'md.part0005.xtc',
                 basedir + 'md.part0006.xtc',
                 basedir + 'md.part0007.xtc',)
print(f'Rank {rank} has loaded the Universe. Now iterating...')

areas = []
for frame in u.trajectory[rank::n_workers]:
    areas.append(u.dimensions[0] * u.dimensions[1])  # in Angstrom squared

print(f'Rank {rank} has finished iterating')
areas = comm.gather(areas, root=0)  # rank 0 gets all areas from everyone.

if rank:  # Everyone but rank 0 quits now
    print(f'Rank {rank} says goodbye!')
    sys.exit()

all_areas = np.concatenate(areas)

print(f'Rank {rank} is writing the output')
np.savetxt('areas_MPI.dat', all_areas)  # we save the array of values for later comparison
print(f'Rank {rank} leaves and closes the door. Bye!')
