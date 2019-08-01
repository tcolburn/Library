#!/usr/bin/env python
# Welcome! This script is for the purpose of DIMs trajectory analysis.
# At the moment, it can calculate and plot Delta RMSD and Gate order
# parameter changes with respect to simulation time, as well as conatenate trajectory files. Be sure to include
# the .crd files for your initial/target structures in the working directory,
# and name your dcds in the form of "dims_mhp1_oi_*"

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import MDAnalysis.analysis.rms as rms
import numpy as np
from itertools import izip
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import os


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
    F[:] = [func(sy[start:stop]) for start,stop in izip(bin_index[:-1],bin_index[1:])]
    return F,bins

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
                for ts in data.trajectory:
                    print "Time step #" + str(n)
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
            print("Wrote topology {} using the selection: {}".format(outname, select)
r

def traj_pare_select(u, outname=None, topoutname=False, topext='.pdb',
                     select="all", b=0, e=None, skip=None):
    selection = u.select_atoms(select)

    # Optionally generate topology from inputted universe
    if topoutname:
        generate_top_from_trj(u, select=select, outname=topoutname, ext=topext)

    # Generate trajectory using "select" as the atom selection
    if outname is None:
        base, ext = os.path.splitext(u.trajectory.filename)
        outname = base + '_' + select + ext
    with mda.Writer(outname, selection.n_atoms) as W:
        print("Writing {} using the selection: {}".format(outname, select))
        for i, ts in enumerate(u.trajectory[b:e:skip]):
            W.write(selection)
            if i % 100 == 0:
                print "Frame {}\r".format(i),
        print "\n"


def traj_pare(traj, out="path.dcd", top="step5_assembly.psf", select="all", b=0, e=None, skip=10):
    u = mda.Universe(top, traj)
    system = u.select_atoms(select)
    with mda.Writer(out, system.n_atoms) as W:
        for i, ts in enumerate(u.trajectory[b:e:skip]):
            W.write(system)
            print "Frame #" + str(i+1)

def gates_DIM_proc(traj, top = "/nfs/homes/tcolburn/Projects/Beckstein/Mhp1/top/2x79_r10_g470.psf"):
    u = mda.Universe(top, traj)
    traj_len = len(u.trajectory)

    gate1 = []
    gate2 = []
    gate3 = []

    for ts in u.trajectory:
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
    print "Tabulating gate distances..."
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

    print "Lists populated"
    return d1, d2, d3

def cat(basename = "dims_mhp1_i2occ2o", trjdir = ".", psf="2jln_r10_g470_c22.psf", n_trj=100, select="protein"):
    out = "full_traj.dcd"
    trj_name_list = []

    # Append all trajectory (filesnames) to complete list of parts
    for i in xrange(n_trj):
        print "writing out frames for trajectory %i" % (i+1)
        # trj_name = basename + str(i+1) + ".dcd"
        trj_name = '{}/{}_{:n}.dcd'.format(trjdir, basename, i+1)
        # trj_name = '%s_%i.dcd' % (basename, i+1)
        if not os.path.exists(trj_name) or os.path.getsize(trj_name) == 0:
            break
        trj_name_list.append(trj_name)

    # Use MDAnalysis (ChainReader) to read in constructed list of trajectory parts
    mhp1 = mda.Universe(psf, trj_name_list)
    ag = mhp1.select_atoms(select)
    print "Universe loaded; writing out trajectory"
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
