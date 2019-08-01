#!usr/bin/env python
# Welcome! This script is for the purpose of DIMs trajectory analysis,
# designed in particular for the transporter Mhp1.

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis.analysis.rms as rms
import MDAnalysis
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.analysis.psa import PSAnalysis
import pandas as pd
import scipy.cluster.hierarchy as hier
from scipy.spatial.distance import directed_hausdorff
import matplotlib.pyplot as plt
from pylab import *
from matplotlib import collections  as mc
import glob
import os
import shutil
import sys
from tqdm import tqdm_notebook as tqdm

# Endpoint structures

target = mda.Universe("/nfs/homes3/tcolburn/Projects/Beckstein/Mhp1/analysis/mk_figs/inter_DIMS/TC/PSA/PSAnalysisTutorial/paper_methods/top/init.psf",
                   "/nfs/homes3/tcolburn/Projects/Beckstein/Mhp1/analysis/mk_figs/inter_DIMS/TC/PSA/PSAnalysisTutorial/paper_methods/EX_10/targ.crd")

initial = mda.Universe("/nfs/homes3/tcolburn/Projects/Beckstein/Mhp1/analysis/mk_figs/inter_DIMS/TC/PSA/PSAnalysisTutorial/paper_methods/top/init.psf",
                   "/nfs/homes3/tcolburn/Projects/Beckstein/Mhp1/analysis/mk_figs/inter_DIMS/TC/PSA/PSAnalysisTutorial/paper_methods/EX_10/init.crd")

# RESIDs adjusted ****(confirm selection rigorously!)****

mhp1_bun_resids = "(resid 32:64 or resid 67:102 or resid 203:300)"

# Reference structure

candidate = '/nfs/homes3/tcolburn/Projects/Beckstein/Mhp1/analysis/mk_figs/inter_DIMS/TC/PSA/PSAnalysisTutorial/IF_OCC_2jlo_avg_candidate.pdb'
ref_IF_OCC_avg = mda.Universe(candidate)

# Selections

rmsd_alignment_selection   =  mhp1_bun_resids + ' and name CA'
rmsd_measurement_selection =  'not name H*'

# Gate info

ECthin_resid_pair = ("resid 47:49", "resid 360:362")
ICthin_resid_pair = ("resid 229:231", "resid 161:163")
thick_resid_pair = ("resid 38 or resid 41", "resid 309 or resid 312:313")

gate_resid_pairs = [ECthin_resid_pair, ICthin_resid_pair, thick_resid_pair]
gate_names = ['EC thin', 'IC thin', 'Thick']


# Misc.

delta_label="$\Delta \Delta RMSD$"

gate_labels = ["Extracellular thin", "Intracellular thin", "Thick gate"]

df_labels = ["EC Thin/Bundle distance", "IC Thin/Bundle distance", "Thick/Bundle distance", "$\Delta \Delta RMSD$"]


###############################################################################




def get_sims(counts, method_names=[], exclude=[],
                 ext='.dcd', methods_basedir = 'paper_methods', topname="", prefix=""):
    """Load data into MDA universes.

    Data are pulled from the filepaths provided and then loaded into their own
    MDA Universe

    :Arguments:

       counts
          The number of sims from each method
       method_names
          list of method names
       exclude
          a list of sims to exclude
       ext
          the type of trajectory file being located
       methods_basedir
          the directory containing the data

    :Returns:

       universes
           a list of universes"""
    labels = []  # Heat map labels
    simulations = []  # List of simulation topology/trajectory filename pairs
    universes = []  # List of MDAnalysis Universes representing simulation\
    exclude = exclude

    for method in tqdm(method_names, position=0):
        # Note: DIMS uses the PSF topology format
        topname = topname
        pathname1 = prefix + ext
        pathname2 = pathname1
        method_dir = '{}/{}'.format(methods_basedir, method)
        sys.stdout.write(method + '|')

        if method is method_names[0]:
            for run in tqdm(range(1,counts[0]+1), position=1): # Number of runs per method
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname1)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is method_names[1] :
            for run in tqdm(range(1,counts[1]+1), position=1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname1)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is method_names[2]:
            for run in tqdm(range(1,counts[2]+1), position=1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is method_names[3]:
            for run in tqdm(range(1,counts[3]+1), position=1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is method_names[4]:
            for run in tqdm(range(1,counts[4]+1), position=1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        if method is method_names[5]:
            for run in tqdm(range(1,counts[5]+1), position=1):
                run_dir = '{}/{:03n}'.format(method_dir, run)
                topology = '{}/{}'.format(methods_basedir, topname)
                trajectory = '{}/{}'.format(run_dir, pathname2)
                labels.append(method + '(' + str(run) + ')')
                simulations.append((topology, trajectory))

        # if method is method_names[6]:
            # for run in range(1,counts[6]+1):
                # run_dir = '{}/{:03n}'.format(method_dir, run)
                # topology = '{}/{}'.format(methods_basedir, topname)
                # trajectory = '{}/{}'.format(run_dir, pathname2)
                # labels.append(method + '(' + str(run) + ')')
                # simulations.append((topology, trajectory))

        # if method is method_names[7]:
            # for run in range(1,counts[7]+1):
                # run_dir = '{}/{:03n}'.format(method_dir, run)
                # topology = '{}/{}'.format(methods_basedir, topname)
                # trajectory = '{}/{}'.format(run_dir, pathname2)
                # labels.append(method + '(' + str(run) + ')')
                # simulations.append((topology, trajectory))

    simulations = [simulations[i] for i in range(len(simulations)) if i not in exclude]
    labels = [labels[i] for i in range(len(labels)) if i not in exclude]

    return simulations, labels


from MDAnalysis import Universe

def get_unis(simulations):
    """Load data into MDA universes.

    Data are pulled from the filepaths provided and then loaded into their own
    MDAnalysis

    :Arguments:

       simulations
          a list of filepaths to the simulation data

    :Returns:

       universes
           a list of universes"""
    universes = []
    for i,sim in enumerate(simulations):
        sys.stdout.write("\rUniverse #" + str(i+1))
        universes.append(Universe(*sim))
    return universes


def separate_sims(simulations, counts, exclude=[]):
    """Separate sims by methods, generating a list of lists of filepaths.

    ** WORK IN PROGRESS **

    :Arguments:

       simulations
          a list of filepaths to the simulation data
       counts
          a list of the number of trajectories per method

    :Returns:

       universes
           a list of universes"""
    # There's definitely a better way to do this. . .
    I = counts[0]
    J = I + counts[1]
    K = J + counts[2]
    L = K + counts[3]
    M = L + counts[4]
    N = M + counts[5]

    EX_10 = [simulations[i] for i in range(0,I) if i not in exclude]
    IM_25  = [simulations[i] for i in range(I,J) if i not in exclude]
    ABMD_10 = [simulations[i] for i in range(J,K) if i not in exclude]
    ABMD_25 = [simulations[i] for i in range(K,L) if i not in exclude]
    TMD_fast_10 = [simulations[i] for i in range(L,M) if i not in exclude]
    TMD_fast_25 = [simulations[i] for i in range(M,N) if i not in exclude]

#     data = [EX_10, IM_25, TMD_fast_10, TMD_fast_25]
    data = [EX_10, IM_25, ABMD_10, ABMD_25, TMD_fast_10, TMD_fast_25]
    return data

def compact_labels(counts):
    """Replace method labels with smaller symbols -- for tight label-fitting.

    ** WORK IN PROGRESS **

    :Arguments:

       counts
          a list of the number of trajectories per method

    :Returns:

       labels
           a list of labels"""

    ex_lab = ["*" for i in range(counts[0])]
    im_lab = ["~" for i in range(counts[1])]
    abmd_10_lab = ["_" for i in range(counts[2])]
    abmd_25_lab = ["__" for i in range(counts[3])]
    tmd_10_lab = ["+" for i in range(counts[4])]
    tmd_25_lab = ["++" for i in range(counts[5])]

    labels = ex_lab + im_lab + abmd_10_lab + abmd_25_lab + tmd_10_lab + tmd_25_lab
    return labels




def rmsd_to_init(atomgroup,
                   initial_atomgroup = initial.atoms,
                   selection = rmsd_measurement_selection,
                   center = True, superposition = True):
    """Calculate RMSD to the initial structure

    :Arguments:

       atomgroup
           the MDA atomgroup associated with a trajectory

       initial_atomgroup
           atomgroup of the initial structure
       selection
           atom selection used for fitting
       center, superposition
           conditions for fitting by MDA (should be Boolean)

    :Returns:

       rmsd
           a scalar RMSD-to-initial value"""
    ag        = atomgroup.select_atoms(selection)
    initial_ag = initial_atomgroup.select_atoms(selection)
    return rmsd(ag.positions, initial_ag.positions, center=center, superposition=superposition)

def rmsd_to_target(atomgroup,
                   target_atomgroup = target.atoms,
                   selection = rmsd_measurement_selection,
                   center = True, superposition = True):
    """Calculate RMSD to the target structure

    :Arguments:

       atomgroup
           the MDA atomgroup associated with a trajectory
       initial_atomgroup
           atomgroup of the target structure
       selection
           atom selection used for fitting

       center, superposition
           conditions for fitting by MDA (should be Boolean)

    :Returns:

       rmsd
           a scalar RMSD-to-target value"""
    ag        = atomgroup.select_atoms(selection)
    target_ag = target_atomgroup.select_atoms(selection)
    return rmsd(ag.positions, target_ag.positions, center=center, superposition=superposition)

def dd_rmsd(atomgroup, initial_atomgroup, target_atomgroup,
                   selection = rmsd_measurement_selection,
                   center = True, superposition = True):

    """Calculate the delta-delta RMSD of the atomgroup. This is a simple
    difference between two previously defined functions.

    :Arguments:

       atomgroup
           the MDA atomgroup associated with a trajectory
       initial_atomgroup
           atomgroup of the initial structure
       selection
           atom selection used for fitting
       center, superposition
           conditions for fitting by MDA (should be Boolean)

    :Returns:

       ddrmsd
           delta-delta RMSD"""
    return rmsd_to_init(atomgroup, initial_atomgroup, selection,
        center, superposition)-rmsd_to_target(atomgroup, target_atomgroup, selection,
        center, superposition)

def ic_thin(atomgroup, ref_select= "resid 229:231", ic_select="resid 161:163"):
    """Calculate the intracellular (IC) gate distance using MDA

    :Arguments:

       atomgroup
           the MDA atomgroup associated with a trajectory
       ref_select
           atom selection for an immobile domain
       ic_select
           atom selection for the IC gate

    :Returns:

       gate_d
           IC gate distance"""
    b_ic = atomgroup.select_atoms(ref_select) # Bundle
    ic_thin = atomgroup.select_atoms(ic_select) # Select intracellular thin gate
    gate_d = np.linalg.norm(b_ic.center_of_mass() - ic_thin.center_of_mass())
    return gate_d


def ec_thin(atomgroup, ref_select="resid 47:49", ec_select="resid 360:362"):
    """Calculate the extracellular (EC) gate distance using MDA

    :Arguments:

       atomgroup
           the MDA atomgroup associated with a trajectory
       ref_select
           atom selection for an immobile domain
       ec_select
           atom selection for the EC gate

    :Returns:

       gate_d
           EC gate distance"""
    b_ec = atomgroup.select_atoms(ref_select) # Bundle
    ec_thin = atomgroup.select_atoms(ec_select) # Select intracellular thin gate
    gate_d = np.linalg.norm(b_ec.center_of_mass() - ec_thin.center_of_mass())
    return gate_d


def thick(atomgroup, bun_select="resid 28 or resid 41", ic_select="resid 309 or resid 312:313"):
    b_thick = atomgroup.select_atoms(bun_select) # Bundle
    thick = atomgroup.select_atoms(ic_select) # Select intracellular thin gate
    gate_d = np.linalg.norm(b_thick.center_of_mass() - thick.center_of_mass())
    return gate_d

def compute_timeseries(universe, frame_function, delta_t=.001, **kwargs):
    """Calculate a time series for a trajectory

    :Arguments:

       universe
           the MDA universe of the trajectory
       frame_function
           some input function that takes in single trajectory frames
       delta_t
           timestep, in picoseconds (ps)

    :Returns:

       pd.Series(timeseries), pd.Series(t_DIMS)
           time series and "physical" time as pandas Series"""
    atoms = universe.atoms
    timeseries = np.zeros(universe.trajectory.n_frames)
    t_DIMS = np.zeros(universe.trajectory.n_frames)
    for i, ts in enumerate(universe.trajectory):
        timeseries[i] = frame_function(atoms, **kwargs)
        t_DIMS[i] = int(i)*delta_t        # Instantaneous "physical" time
        # sys.stdout.write("\rTimestep: [" + str(i+1).zfill(3) + "/" + str(len(universe.trajectory)).zfill(3) + "] ")
    return pd.Series(timeseries), pd.Series(t_DIMS)


def hdistance_calcs(psa_object):
    """Calculates the Hausdorff distance, and generates important
    book-keeping information associated with the Hausdorff pairs.

    ** WORK IN PROGRESS **

    :Arguments:

       numpaths
          the total number of simulations in the ensemble

       psa_object
          the Path Similarity Analysis (PSA) object associated with your analysis

    :Returns:

       D, pair_nums, HP_list, tot_frames

           D: the distance matrix
    #         print "j = " + str(j)

           path_nums: pair numbering corresponding to the list of sim
                      universes fed to PSA during setup.

           HP_list: Hausdorff pair frame numbers, for pulling from the
                    the list of universes

           tot_frames: total number of frames for each trajectory """

    HP_list = [] # Keep track of HPs outside PSA object (redundant)
    tot_frames = [] # Book-keeping of total frames
    path_nums = [] # Simulation numbers for pulling frames

    numpaths = psa_object.npaths
    D = np.zeros((numpaths,numpaths)) # Hausdorff distance matrix
    stride = 1
    N = len(psa_object.paths[0][0]) // 3
    sqrtN = N**0.5

    n = 0
    for i in tqdm(range(0, numpaths-1), position=0):
        for j in tqdm(range(i+1, numpaths), position=1):
            P = psa_object.paths[i][::stride] # Pull array of path info
            Q = psa_object.paths[j][::stride]
            d_pq, a_pq, b_pq = directed_hausdorff(P, Q)
            d_qp, a_qp, b_qp = directed_hausdorff(Q, P)
            if d_pq > d_qp:
                stuff = d_pq, a_pq, b_pq
                path_nums.append((i,j))
                n+=1
            else:
                stuff = d_qp, a_qp, b_qp
                path_nums.append((j,i))
                n+=1
            d, index1, index2 = stuff # (un-normalized) distance, frame1, frame2
            D[i,j] = d / sqrtN
            D[j,i]= D[i,j] # Matrix should be symmetric
            HP_list.append((index1, index2)) # Hpair frames
            tot_frames.append((len(P), len(Q)))
            # sys.stdout.write("\rSim #" + str(i + 1) + " [" + str(j+1) + "/" + str(numpaths) + "]")
    return D, path_nums, HP_list, tot_frames

def ensemble_analysis_df(universe_list, frame_functions, func_names, traj_dir = "df_pkls/", delta_label = r"$\Delta \Delta RMSD$",  prefix="sim_", **kwargs):
    """		Calculate a time series for a trajectory

    :Arguments:

       universe_list
           the list of MDA universes associated with an ensemble
       frame_functions
           some set of input functions for operating on frame data
       delta_t
           timestep, in picoseconds (ps)

    :Returns:

       pd.Series(timeseries), pd.Series(t_DIMS)
           time series and "physical" time as pandas Series"""
    for i, universe in tqdm(enumerate(universe_list), position=0):
        traj_num = str((i+1)).zfill(3)
        df_list = []
        for j, function in tqdm(enumerate(frame_functions), position=1):
            # sys.stdout.write("                          Processing Trajectory " + traj_num + " using function [" + str(j+1) + "/" + str(len(frame_functions)) + "]" + "\r")
            ff_data, physical_time = compute_timeseries(universe, function)
            df_list.append(ff_data.to_frame(name=func_names[j]))
        df_list.append(physical_time.to_frame(name=func_names[-1]))
        Data = pd.concat(df_list, axis=1)
        ddrmsd = Data[func_names[0]] - Data[func_names[1]]
        df_list.append(ddrmsd.to_frame(name=delta_label))
        Data = pd.concat(df_list, axis=1)
        Data = Data.set_index(delta_label)
        pickle_name = traj_dir + prefix + str(i+1).zfill(3) + ".pkl"
        Data.to_pickle(pickle_name)



def regularized_function(x, y, func, bins=100, brange=None):
    """Compute *func()* over data aggregated in bins.
    ``(x,y) --> (x', func(Y'))``  with ``Y' = {y: y(x) where x in x' bin}``
    First the data is collected in bins x' along x and then *func* is
    applied to all data points Y' that have been collected in the bin.
    .. function:: func(y) -> float
       *func* takes exactly one argument, a np 1D array *y* (the
       values in a single bin of the histogram), and reduces it to one
       scalar float.
    .. Note:: *x* and *y* must be 1D arrays.
    :Arguments:
       x
          abscissa values (for binning)
       y
          ordinate values (func is applied)
       func
          a np ufunc that takes one argument, func(Y')
       bins
          number or array
       brange
          limits (used with number of bins)
    :Returns:
       F,edges
          function and edges (``midpoints = 0.5*(edges[:-1]+edges[1:])``)
    (This function originated as
    :func:`recsql.sqlfunctions.regularized_function`.)"""
    _x = np.asarray(x)
    _y = np.asarray(y)

    if len(_x.shape) != 1 or len(_y.shape) != 1:
        raise TypeError("Can only deal with 1D arrays.")

    # setup of bins (taken from np.histogram)
    if (brange is not None):
        mn, mx = brange
        if (mn > mx):
                    raise AttributeError('max must be larger than min in range parameter.')

    if not np.iterable(bins):
        if brange is None:
            brange = (_x.min(), _x.max())
        mn, mx = [float(mi) for mi in brange]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = np.linspace(mn, mx, bins+1, endpoint=True)
    else:
        bins = np.asarray(bins)
        if (np.diff(bins) < 0).any():
            raise ValueError('bins must increase monotonically.')

    sorting_index = np.argsort(_x)
    sx = _x[sorting_index]
    sy = _y[sorting_index]

    # boundaries in SORTED data that demarcate bins; position in bin_index is the bin number
    bin_index = np.r_[sx.searchsorted(bins[:-1], 'left'),
                         sx.searchsorted(bins[-1], 'right')]

    # naive implementation: apply operator to each chunk = sy[start:stop] separately
    #
    # It's not clear to me how one could effectively block this procedure (cf
    # block = 65536 in np.histogram) because there does not seem to be a
    # general way to combine the chunks for different blocks, just think of
    # func=median
    F = np.zeros(len(bins)-1)  # final function
    F[:] = [func(sy[start:stop]) for start,stop in zip(bin_index[:-1],bin_index[1:])]
    return F,bins



def gate_regularize(delta_delta, gate_dists, vmin=-2.85, vmax=2.85,  brange=(-2.85,2.85), func=np.std, bins=50):
        """Calculate a time series for a trajectory

        :Arguments:

           delta_delta
               array-like object of delta-delta data (abscissa)
           gate_dists
               array-like object of gate distances (ordinate)
           func
               a np ufunc that takes one argument, func(Y')
           bins
               number or array
           brange
               limits (used with number of bins)

        :Returns:

           mean, std, zeta
               time series and "physical" time as pandas Series"""
        mean, q1 = regularized_function(delta_delta, gate_dists, np.std, bins=bins, brange=brange)
        std, q1 = regularized_function(delta_delta, gate_dists, np.mean, bins=bins, brange=brange)
        zeta = 0.5*(q1[1:] + q1[:-1])
        return mean, std, zeta

def gates_DIM_proc_adjusted(traj, top = "/nfs/homes/tcolburn/Projects/Beckstein/Mhp1/top/2x79_r10_g470.psf"):
    u = mda.Universe(top, traj)
    traj_len = len(u.trajectory)

    gate1 = []
    gate2 = []
    gate3 = []

    for ts in tqdm(u.trajectory):
        b_ec = u.select_atoms("resid 47:49") #select bundle for ec_thin comparison
        ec_thin = u.select_atoms("resid 360:362") #select extracellular thin gate
        b_ic = u.select_atoms("resid 161:163") #sic
        ic_thin = u.select_atoms("resid 229:231") #select intracellular thin gate
        b_tg = u.select_atoms("resid 28 or resid 41") #sic
        thick = u.select_atoms("resid 309 or resid 312:313") #select thick gate

        gate1_a = abs(b_ec.center_of_mass() - ec_thin.center_of_mass())
        gate2_a = abs(b_ic.center_of_mass() - ic_thin.center_of_mass())
        gate3_a = abs(b_tg.center_of_mass() - thick.center_of_mass())
        gate1_b = np.linalg.norm(gate1_a)
        gate2_b = np.linalg.norm(gate2_a)
        gate3_b = np.linalg.norm(gate3_a)

        gate1.append(gate1_b) #Ec_thin
        gate2.append(gate2_b) #Ic_thin
        gate3.append(gate3_b) #Thick Gate

    return gate1, gate2, gate3, traj_len

def get_axis_limits(ax, scale=.8):
    return ax.get_xlim()[1]*scale - .7, ax.get_ylim()[1]*scale


def gspace_plot(pickle_lists, bbox = (1.1, 0.1), cbar_loc=[0.5, 0.53, 0.03, 0.3],
                df_labels=[], x_freq=0.5,  y_freq=0.5, vmin=-2.98, vmax=2.98,  gnums=[], gate_labels=[],
                method_labels = [], target="explicit/1fs/targ.crd", top="paper_methods/top/init.psf",
                window = [6, 13, 6, 12.5], savename = "full_gspace.pdf", colors=[], fontsize = "11",
                figsize=(10, 7)):
    """Calculate a time series for a trajectory

        :Arguments:

           pickle_lists
               list of pickled dataframes from sims
           bbox
               location of legend
           cbar_loc
               location and width/height of color bar
           df_labels
               dataframe headers
           x_freq
               number of ticks on the unit interval in x
           x_freq
               number of ticks on the unit interval in y
           vmin
               minimum delta-delta RMSD
           vmax
               maximum delta-delta RMSD
           gnums
               the choice of Mhp1 gates to compare;
               numbering corresponds to `df_labels`
           gate_labels
               . . . these are gate labels.
           method_labels
               . . . method labels. . .
           target
               filepath to target structure
           top
               filepath to topology
           window
               x- and y- ranges; list of length 4
           savename
               the name of your plot (include your own extension)
           colors
               limits (used with number of bins)
           fontsize
               limits (used with number of bins)
           figsize
               limits (used with number of bins)

        :Returns:

           mean, std, zeta
               time series and "physical" time as pandas Series"""
    # Plot preparations
    sns.set_style("ticks")
    fig = plt.figure(figsize=figsize)
    method_count = len(pickle_lists)
    sub_plt_nums = [method_count*100 + 20 + (i+1) for i in range(method_count)]
    t1, t2, t3, traj_len = gates_DIM_proc_adjusted(target, top)       # Target processing
    tgate_data = [t1, t2, t3]
    axes = []
    for i in range(method_count):
        ax = plt.subplot(sub_plt_nums[i])
        axes.append(ax)
        ax.axis(window)

        ax.grid(False)
#         ax.set_rasterized(True)

    # Plotting

    for i, p_list in tqdm(enumerate(pickle_lists), position=0):
        # sys.stdout.write("Method #" + str(i+1) + "\r")
        ax = axes[i]
        for j, pickle in enumerate(p_list):
            df = pd.read_pickle(pickle)
            xgate = df[df_labels[gnums[1]]].values
            ygate = df[df_labels[gnums[0]]].values
            ax.plot(xgate, ygate, c='k', zorder=-1, lw=0.75, rasterized=True)
            c =  df.index.values
            if i % 2==0:
                plot = ax.scatter(x=xgate, y=ygate, marker='.', c=c,vmin=vmin, vmax=vmax,  alpha=.21, s=20, cmap="cividis_r", rasterized=True)
            else:
                ax.scatter(x=xgate, y=ygate, marker='.', c=c,vmin=vmin, vmax=vmax,  alpha=.21, s=20, cmap="cividis_r", rasterized=True)
        if i % 2==0:
            ax.set_ylabel(df_labels[gnums[0]] + r" ($\AA$)", size = 10)
        else:
            pass
        ax.scatter(x=xgate[0], y=ygate[0], marker='>', color=sns.xkcd_rgb["black"], s=120, label="Initial", rasterized=True)
        ax.scatter(tgate_data[gnums[1]], tgate_data[gnums[0]], marker='*', color=sns.xkcd_rgb["gold"], s=120, label="Target", rasterized=True)
        xycoords = np.asarray(get_axis_limits(ax)) - np.array([0,-2])
        ax.annotate(method_labels[i], xy=tuple(xycoords))
        ax.xaxis.set_ticks(np.arange(window[0], window[1]+1, x_freq))
        ax.yaxis.set_ticks(np.arange(window[2], window[3]+1, y_freq))

    n_p = len(pickle_lists)

    if n_p == 1:
        axes[0].legend(loc = 'center left', bbox_to_anchor=bbox, prop = {'size':13})
        axes[0].set_xlabel(df_labels[gnums[1]]+ r" ($\AA$)", size = fontsize)
    if n_p > 1 and not n_p % 2 == 0:
        axes[1].legend(loc = 'center left', bbox_to_anchor=(1.1, 0.5), prop = {'size':13})
        axes[-1].set_xlabel(df_labels[gnums[1]]+ r" ($\AA$)", size = fontsize)
    if n_p > 1 and n_p % 2 == 0:
        axes[1].legend(loc = 'center left', bbox_to_anchor=(1.1, 0.5), prop = {'size':13})
        axes[-1].set_xlabel(df_labels[gnums[1]]+ r" ($\AA$)", size = fontsize)
        axes[-2].set_xlabel(df_labels[gnums[1]]+ r" ($\AA$)", size = fontsize)
    else:
        pass

    cbar_ax = fig.add_axes(cbar_loc)
    sns.despine(fig, right=True, top=True)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cb = mpl.colorbar.ColorbarBase(cbar_ax, cmap="cividis_r", norm=norm)
    cb.set_label(r"$\Delta \Delta RMSD$")
    plt.tight_layout()
    plt.savefig(savename, bbox_inches="tight")

def full_delta_plot(pickle_lists, gate_labels=gate_labels, df_labels=df_labels,
                    method_labels = [], init="explicit/1fs/init.crd",
                    target="explicit/1fs/targ.crd", top="paper_methods/top/init.psf",
                    gate_colors = ["blue","red", "black"], window = [-3.2, 3.2, 0, 30],
                    fontsize = "11", figsize=(7, 10), bins=200,
                    savename = "full_delta_delta.pdf"):

    """Calculate a time series for a trajectory

    :Arguments:

       pickles_lists
           a list of lists of the pickle filepaths for each method
       gate_lables
           gate labels
       df_labels
           dataframe column labels
       method_labels
           annotated text for each subplot (i.e. method)
       init
           path to initial structure
       target
           path to target structure
       top
           path to topology
       gate_colors
           gate colors
       window
           plotting window [x_min, x_max, y_min, y_max]
       fontsize
           font size
       figsize
           figure dimensions [(x, y)]
       bins
           number of bins
       savename
           filename for the plot


    :Returns:

       mean, std, zeta
           time series and "physical" time as pandas Series"""

    # Plot preparations
    sns.set_style("ticks")
    fig = plt.figure(figsize=figsize)
    # print(df_labels)
    method_num = len(pickle_lists)
    sub_plt_nums = [method_num*100 + 10 + (i+1) for i in range(method_num)]
    t1, t2, t3, traj_len = gates_DIM_proc_adjusted(target, top)       # Target processing
    tgate_data = [t1, t2, t3]
    i1, i2, i3, traj_len = gates_DIM_proc_adjusted(init, top)         # Initial processing
    igate_data = [i1, i2, i3]
    ep_x, time = delta_proc(init, target, init, top)     # Absolute endpoint delta-delta value
    axes = []
    for i in range(method_num):
        ax = plt.subplot(sub_plt_nums[i])
        axes.append(ax)
        ax.axis(window)
        pax = ax.twinx()
        pax.set_ylim(ax.get_ylim())
        ax.grid(False)

    # Data regularization and organization
    for j, gate in tqdm(enumerate(gate_labels)):
        # sys.stdout.write("Processing data for " + df_labels[j] + ":")
        data_list = []
        for k, p_list in enumerate(pickle_lists):
            # sys.stdout.write("Method " + str(k+1) + "\r")
            method_dfs = []
            for l, pkl in enumerate(p_list):
                df = pd.read_pickle(pkl)
                gate_data = df[df_labels[j]]
                method_dfs.append(gate_data)
            method_tot = pd.concat(method_dfs)
            data_list.append(method_tot)
        for n, ax in tqdm(enumerate(axes), position=2):
            delta_delta = data_list[n].index.values
            gate_dist = data_list[n]
            g_std, q1 = regularized_function(delta_delta, gate_dist, np.std, bins=bins)
            g_mean, q1 = regularized_function(delta_delta, gate_dist, np.mean, bins=bins)
            zeta = 0.5*(q1[1:] + q1[:-1])

            # Plotting
            if n % 2 == 0:
                ax.set_ylabel(r" Gate distance ($\AA$)", size = 10)

            ax.fill_between(zeta, g_mean - g_std, g_mean + g_std,
                            color=sns.xkcd_rgb[gate_colors[j]], alpha=0.35)
            ax.plot(zeta, g_mean, c=sns.xkcd_rgb[gate_colors[j]],
                             label = gate_labels[j], zorder = 0, rasterized=True)
            ax.axhline(tgate_data[0], -2.927, 0, c=sns.xkcd_rgb["black"], linestyle=":", markersize=10)
            ax.axhline(tgate_data[1], -2.927, 0, c=sns.xkcd_rgb["black"], linestyle=":")
            ax.axhline(tgate_data[2], -2.927, 0, c=sns.xkcd_rgb["black"], linestyle=":")
            xycoords = np.asarray(get_axis_limits(ax)) - np.array([0.6,0])
            ax.annotate(method_labels[n], xy=tuple(xycoords), fontsize = fontsize)
            if j==0:
                ax.plot(ep_x[0], igate_data[j], c=sns.xkcd_rgb["black"],
                        zorder = 0, marker=">", label = "Initial",
                        markersize=10, linestyle=":", rasterized=True) # Initial

                ax.plot(ep_x[0]*(-1), tgate_data[j], c=sns.xkcd_rgb["gold"],
                        label = "Target", zorder = 0, marker="*",
                        markersize=10, linestyle=":", rasterized=True) # Target
            else:
                ax.plot(ep_x[0], igate_data[j], c=sns.xkcd_rgb["black"],
                        zorder = 0, marker=">", markersize=10, linestyle=":", ) # Initial
                ax.plot(ep_x[0]*(-1), tgate_data[j], c=sns.xkcd_rgb["gold"],
                        zorder = 0, marker="*", markersize=10, linestyle=":", rasterized=True)
#             ax.set_rasterized(True)

    # Misc.
    axes[0].legend(loc = 'center left', bbox_to_anchor=(1.1, 0.5), prop = {'size':10})
    axes[-1].set_xlabel(r"$\Delta\Delta$ - RMSD", size = 10)
    plt.tight_layout()
    plt.savefig(savename, bbox_inches="tight")

class PDBToBinaryTraj(object):

    def __init__(self, universe, output_type='.dcd', infix=''):
        self.universe = universe
        self.universe.atoms.write('new_top.pdb') # write first frame as topology

        self.frames = self.universe.trajectory
        base, ext = os.path.splitext(self.frames.filename)
        path, name = os.path.split(base)
        self.newname = name + infix + output_type

    def convert(self):
        w = mda.Writer(self.newname, self.frames.n_atoms)
        for ts in self.frames:
            w.write(ts)
        w.close()

def coor_assign(traj, traj_top="../mhp1.psf", traj_sel="protein and not resid 1 and not (atom PROA 2 HN)",
             ref="../ref.crd", ref_top="../step5_assembly.psf",
             ref_sel="protein and not (atom PROA 11 HT1 or atom PROA 11 HT2 or atom PROA 11 HT3)",
             outname="newpath.dcd"):
             # Load universes

            ref = mda.Universe(ref_top, ref)
            data = mda.Universe(traj_top, traj)

            # Make appropriate selections

            ref_select = ref.select_atoms(ref_sel)
            data_select = data.select_atoms(traj_sel)
            N = ref.select_atoms('all').n_atoms

            # Write out a new trajectory
            n=1
            with mda.Writer(outname, N) as W:
                for ts in tqdm(data.trajectory):
                    ref_select.positions = data_select.positions
                    W.write(ref.select_atoms("all"))
                    n+=1

def generate_top_from_trj(u, select='all', outname=None, ext='pdb', frame=0):
    selection = u.select_atoms(select)
    if outname:
        u.trajectory[frame]  # Start writing from first frame
        if type(outname) is not str:
            base, ext = os.path.splitext(u.trajectory.filename)
            outname = base + '_' + select + '.' + topext
        with mda.Writer(outname, selection.n_atoms) as W:
            W.write(selection)
            # sys.stdout.write("Wrote topology {} using the selection: {}".format(outname, select))

def traj_pare_select(u, outname=None, topoutname=False, topext='.pdb', select="all", b=0, e=None, skip=None):
    selection = u.select_atoms(select)

         # Optionally generate topology from inputted universe
    if topoutname:
        generate_top_from_trj(u, select=select, outname=topoutname, ext=topext)

         # Generate trajectory using "select" as the atom selectioni
    if outname is None:
        base, ext = os.path.splitext(u.trajectory.filename)
        outname = base + '_' + select + ext
    with mda.Writer(outname, selection.n_atoms) as W:
        for i, ts in tqdm(enumerate(u.trajectory[b:e:skip])):
            W.write(selection)
            # if i % 100 == 0:
                # sys.stdout.write("Writing {} using the selection: {} ; Frame {}\r".format(outname, select, str(i+1).zfill(3)))

def traj_pare(traj, out="path.dcd", top="step5_assembly.psf", select="all", b=0, e=None, skip=10):
    u = mda.Universe(top, traj)
    system = u.select_atoms(select)
    with mda.Writer(out, system.n_atoms) as W:
        for i, ts in tqdm(enumerate(u.trajectory[b:e:skip])):
            W.write(system)
            # sys.stdout.write("Frame #" + str(i+1))

def gates_DIM_proc(traj, top = "/nfs/homes/tcolburn/Projects/Beckstein/Mhp1/top/2x79_r10_g470.psf"):
    u = mda.Universe(top, traj)
    traj_len = len(u.trajectory)

    gate1 = []
    gate2 = []
    gate3 = []

    for ts in tqdm(u.trajectory):
        b_ec = u.select_atoms("resid 38:40") #select bundle for ec_thin comparison
        ec_thin = u.select_atoms("resid 351:353") #select extracellular thin gate
        b_ic = u.select_atoms("resid 152:154") #sic
        ic_thin = u.select_atoms("resid 220:222") #select intracellular thin gate
        b_tg = u.select_atoms("resid 29 or resid 32") #sic
        thick = u.select_atoms("resid 300 or resid 303:304") #select thick gate

        gate1_a = abs(b_ec.center_of_mass() - ec_thin.center_of_mass())
        gate2_a = abs(b_ic.center_of_mass() - ic_thin.center_of_mass())
        gate3_a = abs(b_tg.center_of_mass() - thick.center_of_mass())
        gate1_b = np.linalg.norm(gate1_a)
        gate2_b = np.linalg.norm(gate2_a)
        gate3_b = np.linalg.norm(gate3_a)

        gate1.append(gate1_b) #Ec_thin
        gate2.append(gate2_b) #Ic_thin
        gate3.append(gate3_b) #Thick Gate

    return gate1, gate2, gate3, traj_len

# This adjusted form of the gates function is an extremely sloppy fix due to a hard deadline -- DELETE WHEN FINISHED

def gates_DIM_proc_adjusted(traj, top = "/nfs/homes/tcolburn/Projects/Beckstein/Mhp1/top/2x79_r10_g470.psf"):
    u = mda.Universe(top, traj)
    traj_len = len(u.trajectory)

    gate1 = []
    gate2 = []
    gate3 = []

    for ts in u.trajectory:
        b_ec = u.select_atoms("resid 47:49") #select bundle for ec_thin comparison
        ec_thin = u.select_atoms("resid 360:362") #select extracellular thin gate
        b_ic = u.select_atoms("resid 161:163") #sic
        ic_thin = u.select_atoms("resid 229:231") #select intracellular thin gate
        b_tg = u.select_atoms("resid 28 or resid 41") #sic
        thick = u.select_atoms("resid 309 or resid 312:313") #select thick gate

        gate1_a = abs(b_ec.center_of_mass() - ec_thin.center_of_mass())
        gate2_a = abs(b_ic.center_of_mass() - ic_thin.center_of_mass())
        gate3_a = abs(b_tg.center_of_mass() - thick.center_of_mass())
        gate1_b = np.linalg.norm(gate1_a)
        gate2_b = np.linalg.norm(gate2_a)
        gate3_b = np.linalg.norm(gate3_a)

        gate1.append(gate1_b) #Ec_thin
        gate2.append(gate2_b) #Ic_thin
        gate3.append(gate3_b) #Thick Gate

    return gate1, gate2, gate3, traj_len

def delta_proc(init, targ, traj, top = "/nfs/homes/tcolburn/Projects/Beckstein/Mhp1/top/2x79_r10_g470.psf",
               direction = "i2occ2o", label = "TC"):
    initial = mda.Universe(top, init)
    target = mda.Universe(top, targ)
    if type(traj) is str:
        trajectory = mda.Universe(top, traj)
    else:
        trajectory = traj

    r_init =  rms.RMSD(trajectory, initial, select='name CA and protein')
    r_init.run()
    r_init.save("r_init_" + label)
    r_targ =  rms.RMSD(trajectory, target, select='name CA and protein')
    r_targ.run()
    r_targ.save("r_targ" + label)

    rmsd_init = r_init.rmsd.T
    rmsd_targ = r_targ.rmsd.T
    del_rmsd = rmsd_init[2] - rmsd_targ[2]
    time = rmsd_init[1]

    return del_rmsd, time

def OP_proc(traj, list_name = "d", top = "2jln_r10_g470_c22.psf"):
    u = mda.Universe(top, traj)

    d1 = []
    d2 = []
    d3 = []
    sys.stdout.write("Tabulating gate distances...")
    for ts in u.trajectory:
        b_ec = u.select_atoms("resid 38:40") #select bundle for ec_thin comparison
        ec_thin = u.select_atoms("resid 351:353") #select extracellular thin gate
        b_ic = u.select_atoms("resid 152:154") #sic
        ic_thin = u.select_atoms("resid 220:222") #select intracellular thin gate
        b_tg = u.select_atoms("resid 29 or resid 32") #sic
        thick = u.select_atoms("resid 300 or resid 303:304") #select thick gate

        d1_a = abs(b_ec.center_of_mass() - ec_thin.center_of_mass())
        d2_a = abs(b_ic.center_of_mass() - ic_thin.center_of_mass())
        d3_a = abs(b_tg.center_of_mass() - thick.center_of_mass())
        d1_b = np.linalg.norm(d1_a)
        d2_b = np.linalg.norm(d2_a)
        d3_b = np.linalg.norm(d3_a)

        d1.append(d1_b) #Ec_thin
        d2.append(d2_b) #Ic_thin
        d3.append(d3_b) #Thick Gate

    sys.stdout.write("Lists populated")
    return d1, d2, d3

def cat(basename = "dims_mhp1_i2occ2o", trjdir = ".", psf="2jln_r10_g470_c22.psf", n_trj=100, select="protein", outname="full_traj.dcd", ext=".dcd"):
    out = outname
    trj_name_list = []

    # Append all trajectory (filesnames) to complete list of parts
    for i in tqdm(range(n_trj)):
        # sys.stdout.write("writing out frames for trajectory %i" % (i+1))
        # trj_name = basename + str(i+1) + ".dcd"
        name_string='{}/{}_{:n}' + ext
        trj_name = name_string.format(trjdir, basename, i)
        # trj_name = '%s_%i.dcd' % (basename, i+1)
        if not os.path.exists(trj_name) or os.path.getsize(trj_name) == 0:
            break
        trj_name_list.append(trj_name)

    # Use MDAnalysis (ChainReader) to read in constructed list of trajectory parts
    mhp1 = mda.Universe(psf, trj_name_list)
    ag = mhp1.select_atoms(select)
    # sys.stdout.write("Universe loaded; writing out trajectory")
    with mda.Writer(out, ag.n_atoms) as W:
        for ts in mhp1.trajectory:
            W.write(ag)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Analyze Mhp1 DIMs trajectories.",
                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    parser_delta = subparsers.add_parser('delta', help = "Plots the change in RMSD between the current"
                                                     " configuration and the initial/target configurations")
    parser_delta.add_argument('init', help = 'Initial protein configuration')
    parser_delta.add_argument('targ', help = 'Target protein configuration')
    parser_delta.add_argument('traj', help = 'Trajectory')
    parser_delta.set_defaults(func=delta)

    parser_opar = subparsers.add_parser('gates', help = "Plots gate-specific order parameters:"
     " distance between centers of mass of"
     " a particular residue selections on"
     " both the bundle domain and the gates"
     " of interest.")
    parser_opar.add_argument('traj', help = 'Trajectory')
    parser_opar.add_argument('top', help = 'Topology file')
    parser_delta.set_defaults(func=gates)

    parser_opar = subparsers.add_parser('cat', help = "Concatenates trajectory files")
    parser_opar.add_argument('first', help = 'First dcd')
    parser_opar.add_argument('n', help = 'Total number of non-zero dcds')
    parser_delta.set_defaults(func=cat)

    args = parser.parse_args()
