#!/usr/bin/env python

usage=(
"""Create a single DCD trajectory from any number of input trajectories
in which the NapA dimer is fit to the dimerization domain.""")

import sys
sys.path.append('/nfs/homes3/dldotson/Projects/Transporters/lib/python/')
from aaTransport.Definitions import *
import argparse
import MDAnalysis as md
import MDAnalysis.analysis.align as aln
from MDAnalysis.core.log import ProgressMeter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=usage,
           # formatter_class=argparse.RawDescriptionHelpFormatter,
            )
    parser.add_argument("-r", "--ref", metavar="REF")
    parser.add_argument("top", metavar="TOP")
    parser.add_argument("traj", metavar="TRAJ", nargs='+')
    parser.add_argument("-o", "--outtraj", default="transition.dcd")
    args = parser.parse_args()

    # create universe with trajectory elements; make reference universe
    u = md.Universe(args.top, args.traj)
    uref = md.Universe(args.ref)

    # establish resnums
    A = u.selectAtoms('segid PROA')
    B = u.selectAtoms('segid PROB')
    A.residues.set_resnum(A.residues.resids() + 2)
    B.residues.set_resnum(A.residues.resids() + 2)

    Aref = uref.selectAtoms('segid PROA')
    Bref = uref.selectAtoms('segid PROB')
    Aref.residues.set_resnum(Aref.residues.resids() + 2)
    Bref.residues.set_resnum(Aref.residues.resids() + 2)

    # get selections of dimerization domain
    dimersel = [ NapA_ss[x] for x in NapA_dimer ]
    dimer = u.selectAtoms(*dimersel)
    dimerref = uref.selectAtoms(*dimersel)

    # loop through trajectory fitting to ref, write out
    w = md.Writer(args.outtraj, u.trajectory.numatoms)
    pm = ProgressMeter(u.trajectory.numframes, interval=100)
    for ts in u.trajectory:
        pm.echo(ts.frame)
        aln.alignto(dimer, dimerref)
        w.write(ts)

    w.close_trajectory()
