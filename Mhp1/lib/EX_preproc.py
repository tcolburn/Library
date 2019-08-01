#!/usr/bin/env python
##!/nfs/homes3/sseyler/.virtualenvs/research_stable/bin/python2.7

import MDAnalysis as mda
import os
import shutil
import sys

top = sys.argv[1]
traj = sys.argv[2]
out = sys.argv[3]

print "Processing path"
u = mda.Universe(top, traj)
protein = u.select_atoms("protein")
with mda.Writer(out, protein.n_atoms) as W:
    for ts in u.trajectory:
        W.write(protein)
print "Protein-only path written to current directory"
with open("extract_protein_traj.out", "a") as text_file:
    text_file.write('Loaded {}  |  writing to:\n  >>> {}\n'.format(traj, out))
