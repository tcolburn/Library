#!/usr/bin/env python
# Welcome! This script is for the purpose of DIMs trajectory analysis.
# At the moment, it can calculate and plot Delta RMSD and Gate order
# parameter changes with respect to simulation time, as well as conatenate trajectory files. Be sure to include
# the .crd files for your initial/target structures in the working directory,
# and name your dcds in the form of "dims_mhp1_oi_*"

import MDAnalysis as mda
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis.analysis.rms as rms
import MDAnalysis.analysis.align as align 
from MDAnalysis import Universe
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.analysis.psa import PSAnalysis
from pair_id import PairID
import scipy.cluster.hierarchy as hier
from scipy.spatial.distance import directed_hausdorff
import matplotlib.pyplot as plt
from pylab import *
import MDAnalysis as mda
import numpy as np
import matplotlib as mpl
import MDAnalysis.analysis.rms as rms
import matplotlib.pyplot as plt
import seaborn as sns
import pylab as pl
from matplotlib import collections  as mc
import pandas as pd
import glob
import os




def get_sims(counts, method_names=['EX_10', 'IM_25', 'ABMD_10', 
                'ABMD_25', 'TMD_fast_10', 'TMD_fast_25'], exclude=[], 
                 ext='.dcd', methods_basedir = 'paper_methods'):
     """Fetch simulation filepaths by method, with the possibility of exclusion.

    Data are pulled from some specified directory, with the subdirectories named 
    via `method_names`. A list of excluded trajectories can be specified to eliminate 
    bad runs from the batch, or else specify some other subset of the data. Labels are 
    generated using the `method_names` supplied.
    
    :Arguments:
    
       counts
          a list of the number of trajectories per method
          
       method_names
          a list of the names of each method

       exclude
          a list of excluded trajectory numbers, where numbering corresponds 
          to the order of the `simulations` list returned by the function
          
       methods_basedir
          filepath to the data 
          
       ext
          file extension

    :Returns:
    
       simulations, labels
          function and edges (``midpoints = 0.5*(edges[:-1]+edges[1:])``)
    (This function originated as
    :func:`recsql.sqlfunctions.regularized_function`.)"""

    labels = []  # Heat map labels
    simulations = []  # List of simulation topology/trajectory filename pairs
    universes = []  # List of MDAnalysis Universes representing simulation\
    exclude = exclude
    
    for method in method_names:
        # Note: DIMS uses the PSF topology format
        topname = 'top/init.psf'
        pathname1 = 'path' + ext
        pathname2 = pathname1
        method_dir = '{}/{}'.format(methods_basedir, method)
        print method + '|',   
        
        if method is 'EX_10':
            for run in xrange(1,counts[0]+1): # Number of runs per method
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname1)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is 'IM_25' :
            for run in xrange(1,counts[1]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname1)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is 'ABMD_10':
            for run in xrange(1,counts[2]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))
                
        if method is 'ABMD_25':
            for run in xrange(1,counts[3]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))   
                
        if method is 'TMD_fast_10':
            for run in xrange(1,counts[4]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))  
                
        if method is 'TMD_fast_25':
            for run in xrange(1,counts[5]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))
                
        if method is 'TMD_slow_10':
            for run in xrange(1,counts[6]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory)) 
                
        if method is 'TMD_slow_25':
            for run in xrange(1,counts[7]+1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

    simulations = [simulations[i] for i in xrange(len(simulations)) if i not in exclude]
    labels = [labels[i] for i in xrange(len(labels)) if i not in exclude]
    
    return simulations, labels


