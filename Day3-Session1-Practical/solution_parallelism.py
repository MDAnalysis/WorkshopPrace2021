#!/usr/bin/env python
# coding: utf-8
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np
from multiprocessing import Pool
import itertools
import sys
import os

basedir = '/project/jhlsrf005/JHL_data/Day3-Session1-Practical/'
u = mda.Universe(basedir + 'md.tpr',
                 basedir + 'md.xtc',
                 basedir + 'md.part0002.xtc',
                 basedir + 'md.part0003.xtc',
                 basedir + 'md.part0004.xtc',
                 basedir + 'md.part0005.xtc',
                 basedir + 'md.part0006.xtc',
                 basedir + 'md.part0007.xtc',)


peptides = u.select_atoms('protein')
popes = u.select_atoms('resname POPE')

class HeadgroupZAngle(AnalysisBase):
    def __init__(self, peptides, lipids, binwidth=5.0, **kwargs):  # If you don't know what **kwargs does, it's ok to ignore it for now
        super().__init__(peptides.universe.trajectory, **kwargs)
        self.peptides = peptides
        self.u = self.peptides.universe
        self.lipids = lipids
        
        # Histogram stuff
        self.binwidth = binwidth # The binwidth of the angle distributions
        # The distributions will be histograms between 0 and 180 degrees
        self.bins = np.arange(0, 180 + self.binwidth, self.binwidth)

    def _prepare(self):
        self.results.close_angles = np.zeros(len(self.bins)-1, dtype=int)
        self.results.far_angles = np.zeros(len(self.bins)-1, dtype=int)

        # We need to split the lipid headgroups into top and bottom groups
        #  so that we know for which the angle distribution must be reversed.
        # At the first frame the membrane is quite flat. We can assign
        #  top/bottom leaflets simply by comparing the phosphate position
        #  relative to the membrane's center-of-geometry in z at that frame.
        self.u.trajectory[0]
        membrane_zcog = self.lipids.center_of_geometry()[2]
        self.top_heads = self.lipids.select_atoms('name PO4 NH3 and same residue as '
                                                  f'name PO4 and prop z >= {membrane_zcog}')
        self.bottom_heads = self.lipids.select_atoms('name PO4 NH3 and same residue as '
                                                     f'name PO4 and prop z < {membrane_zcog}')

    def compute_angle(self, headgroups, sign=1):
        # In this topology, NH3 atoms immediately precede PO4 atoms
        pos_nh3 = headgroups.positions[::2]
        pos_po4 = headgroups.positions[1::2]
        vecs = pos_nh3 - pos_po4
        # a correction for vectors crossing the PBC (assumes an orthogonal box)
        vecs += self.u.dimensions[:3]/2
        vecs = mda.lib.distances.apply_PBC(vecs, self.u.dimensions)
        vecs -= self.u.dimensions[:3]/2

        norms = np.linalg.norm(vecs, axis=1)
        angles = sign * np.rad2deg(np.arccos(vecs[:,2]/norms))
        
        # Histogramming the angles over the bins
        dist = np.histogram(angles, bins=self.bins)[0]
        return dist
    
    def _single_frame(self):
        close_heads = self.lipids.select_atoms('name PO4 NH3 and '
                                               'same residue as around 10 global group peptides',
                                                peptides=self.peptides)
        close_heads_top = self.top_heads & close_heads
        close_heads_bottom = self.bottom_heads & close_heads
        far_heads_top = self.top_heads - close_heads_top
        far_heads_bottom = self.bottom_heads - close_heads_bottom
        
        # For each frame we just add up the histograms. We normalize at the end, in _conclude.
        self.results.close_angles += self.compute_angle(close_heads_top)
        self.results.close_angles += self.compute_angle(close_heads_bottom, sign=-1)
        self.results.far_angles += self.compute_angle(far_heads_top)
        self.results.far_angles += self.compute_angle(far_heads_bottom, sign=-1)

    def _conclude(self):
        # Normalization to a probability density
        self.results.close_angles = self.results.close_angles / (self.results.close_angles.sum() * self.binwidth)
        self.results.far_angles = self.results.far_angles / (self.results.far_angles.sum() * self.binwidth)

angle_analysis = HeadgroupZAngle(peptides, popes)

xs = angle_analysis.bins[:-1]

# The 'parallelize_run' wrapper function will be executed by
#  each child process.
# It takes as arguments the HeagroupZAngle analysis object,
#  the number of workers, and the id of the current worker, which
#  it uses to decide which frames to work on.
def parallelize_run(analysis, n_workers, worker_id):
    # To avoid too many progress bars we only switch verbosity on for
    #  the parallel worker with worker_id == 0
    analysis.run(start=worker_id, step=n_workers, verbose=not worker_id)
    return analysis

# How many subprocesses do we want to spawn?
# In a non-HPC use, we can typically let multiprocessing use all
#  the cores in a machine. But here we may be limited by cluster
#  constraints. os.sched_getaffinity(0) tells us which computer
#  cores are available to us.
n_workers = len(os.sched_getaffinity(0))

# We need to build a generator of the set of arguments per worker,
#  meaning we must repeat the analysis and the number of workers.
params = zip(itertools.repeat(angle_analysis),
             itertools.repeat(n_workers),
             range(n_workers))

pool = Pool(processes=n_workers)

# starmap is the call that makes each child run the parallelize_run function.
# starmap will take one tuple at a time from the 'params' generator and pass it to
#  one of the children, which then uses the tuple's elements as the three separate
#  arguments to 'parallelize_run'. starmap won't wait for a child to finish before
#  assigning the next tuple, but after all have been assigned, it waits for
#  every assigned child to finish working. Because we made 'params' to have as
#  many entries as there are children, everyone should have exactly one tuple to work on.
# Read more in the documentation for multiprocessng:
#  https://docs.python.org/3/library/multiprocessing.html
analyses = pool.starmap(parallelize_run, params)

# Free the children's resources
pool.close()

n_frames = [partial_analysis.n_frames for partial_analysis in analyses]

close_angles = np.average([partial_analysis.results.close_angles
                           for partial_analysis in analyses],
                          weights=n_frames,
                          axis=0)

far_angles = np.average([partial_analysis.results.far_angles
                         for partial_analysis in analyses],
                        weights=n_frames,
                        axis=0)

# We save this as a text file for later comparison
np.savetxt('output_multiprocessing_SLURM.dat', np.column_stack((xs, close_angles, far_angles)))

