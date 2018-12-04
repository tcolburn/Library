#!/usr/bin/env python
# Welcome! This script is for the purpose of DIMs trajectory analysis.
# At the moment, it can calculate and plot Delta RMSD and Gate order
# parameter changes with respect to simulation time, as well as conatenate trajectory files. Be sure to include
# the .crd files for your initial/target structures in the working directory,
# and name your dcds in the form of "dims_mhp1_oi_*"

import numpy
from itertools import izip
import os
import MDAnalysis as mda
import matplotlib as mpl
import MDAnalysis.analysis.rms as rms
import matplotlib.pyplot as plt
import seaborn as sns

def regularized_function(x, y, func, bins=100, brange=None):
    """Compute *func()* over data aggregated in bins.
    ``(x,y) --> (x', func(Y'))``  with ``Y' = {y: y(x) where x in x' bin}``
    First the data is collected in bins x' along x and then *func* is
    applied to all data points Y' that have been collected in the bin.
    .. function:: func(y) -> float
       *func* takes exactly one argument, a numpy 1D array *y* (the
       values in a single bin of the histogram), and reduces it to one
       scalar float.
    .. Note:: *x* and *y* must be 1D arrays.
    :Arguments:
       x
          abscissa values (for binning)
       y
          ordinate values (func is applied)
       func
          a numpy ufunc that takes one argument, func(Y')
       bins
          number or array
       brange
          limits (used with number of bins)
    :Returns:
       F,edges
          function and edges (``midpoints = 0.5*(edges[:-1]+edges[1:])``)
    (This function originated as
    :func:`recsql.sqlfunctions.regularized_function`.)"""
    _x = numpy.asarray(x)
    _y = numpy.asarray(y)

    if len(_x.shape) != 1 or len(_y.shape) != 1:
        raise TypeError("Can only deal with 1D arrays.")

    # setup of bins (taken from numpy.histogram)
    if (brange is not None):
        mn, mx = brange
        if (mn > mx):
                    raise AttributeError('max must be larger than min in range parameter.')

    if not numpy.iterable(bins):
        if brange is None:
            brange = (_x.min(), _x.max())
        mn, mx = [float(mi) for mi in brange]
        if mn == mx:
            mn -= 0.5
            mx += 0.5
        bins = numpy.linspace(mn, mx, bins+1, endpoint=True)
    else:
        bins = numpy.asarray(bins)
        if (numpy.diff(bins) < 0).any():
            raise ValueError('bins must increase monotonically.')

    sorting_index = numpy.argsort(_x)
    sx = _x[sorting_index]
    sy = _y[sorting_index]

    # boundaries in SORTED data that demarcate bins; position in bin_index is the bin number
    bin_index = numpy.r_[sx.searchsorted(bins[:-1], 'left'),
                         sx.searchsorted(bins[-1], 'right')]

    # naive implementation: apply operator to each chunk = sy[start:stop] separately
    #
    # It's not clear to me how one could effectively block this procedure (cf
    # block = 65536 in numpy.histogram) because there does not seem to be a
    # general way to combine the chunks for different blocks, just think of
    # func=median
    F = numpy.zeros(len(bins)-1)  # final function
    F[:] = [func(sy[start:stop]) for start,stop in izip(bin_index[:-1],bin_index[1:])]
    return F,bins



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
        gate1_b = numpy.linalg.norm(gate1_a)
        gate2_b = numpy.linalg.norm(gate2_a)
        gate3_b = numpy.linalg.norm(gate3_a)

        gate1.append(gate1_b) #Ec_thin
        gate2.append(gate2_b) #Ic_thin
        gate3.append(gate3_b) #Thick Gate

    return gate1, gate2, gate3, traj_len



def delta_proc(init, targ, traj, top = "/nfs/homes/tcolburn/Projects/Beckstein/Mhp1/top/2x79_r10_g470.psf",
               direction = "i2occ2o", label = "TC"):
    initial = mda.Universe(top, init)
    target = mda.Universe(top, targ)
    trajectory = mda.Universe(top, traj)

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
        d1_b = numpy.linalg.norm(d1_a)
        d2_b = numpy.linalg.norm(d2_a)
        d3_b = numpy.linalg.norm(d3_a)

        d1.append(d1_b) #Ec_thin
        d2.append(d2_b) #Ic_thin
        d3.append(d3_b) #Thick Gate

    print "Lists populated"
    return d1, d2, d3

def cat(direction = "i2occ2o", trjdir = ".", psf="2jln_r10_g470_c22.psf", n_trj=100, select="protein"):
    out = "full_dims_mhp1_" + direction + ".dcd"
    trj_name_list = []
    dcd_basename = "dims_mhp1_" + direction

    # Append all trajectory (filesnames) to complete list of parts
    for i in xrange(n_trj):
        print "writing out frames for trajectory %i" % (i+1)
        # trj_name = dcd_basename + str(i+1) + ".dcd"
        trj_name = '{}/{}_{:n}.dcd'.format(trjdir, dcd_basename, i+1)
        # trj_name = '%s_%i.dcd' % (dcd_basename, i+1)
        if not os.path.exists(trj_name) or os.path.getsize(trj_name) == 0:
            break
        trj_name_list.append(trj_name)

    # Use MDAnalysis (ChainReader) to read in constructed list of trajectory parts
    mhp1 = mda.Universe(psf, trj_name_list)
    ag = mhp1.select_atoms(select)
    with mda.Writer(trjdir + out, ag.n_atoms) as W:
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

