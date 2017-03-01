#!/usr/bin/env python
# Welcome! This script is for the purpose of DIMs trajectory analysis.
# At the moment, it can calculate and plot Delta RMSD and Gate order
# parameter changes with respect to simulation time, as well as conatenate trajectory files. Be sure to include
# the .crd files for your initial/target structures in the working directory,
# and name your dcds in the form of "dims_mhp1_oi_*"

import os
import MDAnalysis as mda
import numpy as np
import matplotlib as mpl
import MDAnalysis.analysis.rms as rms
import matplotlib.pyplot as plt
import seaborn as sns

def delta(init, targ, traj,top = "2jln_r10_g470_c22.psf", direction = "i2occ2o", show = False):
    initial = mda.Universe(top, init)
    target = mda.Universe(top, targ)
    trajectory = mda.Universe(top, traj)

    r_init =  rms.RMSD(trajectory, initial, select='all')
    r_init.run()
    r_init.save("r_init")
    r_targ =  rms.RMSD(trajectory, target, select='all')
    r_targ.run()
    r_targ.save("r_targ")

    rmsd_init = r_init.rmsd.T
    rmsd_targ = r_targ.rmsd.T
    del_rmsd = rmsd_init[2] - rmsd_targ[2]
    time = rmsd_init[1]
    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(111)
    ax.set_title('$\Delta \Delta$ RMSD: OB Gates')
    ax.plot(time, del_rmsd, 'k--')
    ax.set_xlabel("Time Step")
    ax.set_ylabel(r"Delta RMSD ($\AA$)")
    fig.savefig("del_rmsd_" + direction + "_mhp1.pdf")
    if show:
        plt.show()

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

def plt_gates(path_ID = "i2occ2o"):
    fig = plt.figure(figsize = (8,5))
    ax = fig.add_subplot(111)
    ax.plot(d1, 'k-', label = "Extracellular")
    ax.plot(d2, 'b-', label = "Intracellular")
    ax.plot(d3, 'r-', label = "Thick")
    box = ax.get_position()
    ax.set_title('Gate Order Parameters')
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylim([0, 25])
    ax.set_xlabel("timestep")
    ax.set_ylabel(r"Distance to Bundle ($\AA$)")
    fig.savefig("mhp1_gate_tseries_" + path_ID + "_mhp1.pdf")

    print "Gate lists are now populated"

    fig = plt.figure(figsize = (8,5))
    ax = fig.add_subplot(111)
    ax.plot(d1, 'k-', label = "Extracellular")
    ax.plot(d2, 'b-', label = "Intracellular")
    ax.plot(d3, 'r-', label = "Thick")
    box = ax.get_position()
    ax.set_title('Gate Order Parameters')
    ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylim([0, 25])
    ax.set_xlabel("timestep")
    ax.set_ylabel(r"Distance to Bundle ($\AA$)")
    fig.savefig("mhp1_gate_tseries_" + path_ID + "_mhp1.pdf")

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
