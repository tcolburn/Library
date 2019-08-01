from CodeTiming import Timer
import numpy as np
from MDAnalysis import Universe
import MDAnalysis.analysis.contacts
import os

struct = '4ake'
#################################################
# DIMS
#################################################
basedir = '/nfs/homes/sseyler/Simulations/adk/charmm/dims/c_o'
# basedir = '/nfs/homes/sseyler/Simulations/adk/charmm/emin/' + struct

# top_dims = basedir + '/struct/' + 'step1_pdbreader.psf'
top_dims = basedir + '/struct/' + 'open_final.psf'
trj_dims = basedir + '/trj/raw' + '/'

dims_tr = []
for i in xrange(1, 201):
    name = 'dims%04i.dcd' % i
    dims_tr.append(trj_dims + name)

# for i in xrange(1, 200):
#     name = 'dims%04i.dcd' % i
#     if os.path.exists(basedir + name): dims_tr.append(basedir + name)
#     else: pass
#     dims_tr.append(basedir + name)
# with Timer(verbose=True) as t:
#     pgdims1 = PG.PathGroup(dims_tr, topology=top_dims, selection="name CA", \
#         ref_struct=refavgCORE_ca, processes=3)

#################################################
# GP
#################################################
# nt = 40
# basedir = '/nfs/homes/sseyler/Simulations/adk/pathways/forward'
# top_gp = basedir + '1/run41/adk_a_top.pdb'
# filename = 'path_a.dcd'
# gp_tr = []
# for step in xrange(4,5):
#     begin = 40*step+1
#     for run in xrange(begin, begin+nt):
#         filepath = "%s%i/run%i/%s" % (basedir, step, run, filename)
#         gp_tr.append(filepath)


#################################################
# Minimize RMSD DIMS trajs w.r.t. reference
#################################################
if True:
    # from CodeTiming import Timer
    from MDAnalysis import Universe
    from MDAnalysis.analysis.align import *
    import multiprocessing as mp

    # Get raw trajectories
    basedir = '/nfs/homes/sseyler/Simulations/adk/charmm/dims/c_o'
    top_dims = basedir + '/struct/' + 'open_final.psf'
    trj_dims = basedir + '/trj/raw' + '/'
    dims_tr = []
    for i in xrange(1, 201):
        name = 'dims%04i.dcd' % i
        dims_tr.append(trj_dims + name)

    # Input trajectories
    trjs = dims_tr
    top = top_dims
    ref_struct = '4ake'
    outfp = '/nfs/homes/sseyler/Simulations/adk/charmm/dims/c_o/trj/fit' +    \
        '/' + ref_struct + '/'

    # Get reference structure
    prefix = '/nfs/homes/sseyler/Dropbox/Beckstein/data/AdK'
    reffpath = prefix + '/boundary_conformations'
    ref_c = Universe(reffpath+'/1AKE_A.pdb')
    ref_o = Universe(reffpath+'/4AKE_A.pdb')
    ref_c.atoms.translate(-ref_c.atoms.centerOfMass())
    ref_o.atoms.translate(-ref_o.atoms.centerOfMass())

    def rmsd_fit(top, filenames):
        for filename in filenames:
            fp, sep, fn = filename.rpartition("/")
            lhs, sep, rhs = fn.rpartition(".")
            newname = lhs + "_fit-" + ref_struct + sep + rhs
            print "Writing file " + newname + " ..."
            u = Universe(top, filename)
            rms_fit_trj(u, ref_o, select='name CA', filename=outfp + newname)

    # with Timer(verbose=True) as t:
    ntrajs = len(trjs)
    nproc = mp.cpu_count()
    chunksize = ntrajs/nproc + 1 # integer division
    jobs = []
    s = e = 0
    for i in xrange(nproc):
        s = i*chunksize
        e = s + chunksize if e <= ntrajs else None
        p = mp.Process(target=rmsd_fit, args=(top, trjs[s:e]))
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()

#################################################


#################################################
# Setup NCA
#################################################
import numpy as np
from MDAnalysis.analysis.contacts import ContactAnalysis1
from MDAnalysis.analysis.align import *
from MDAnalysis import Universe

basedir = '/nfs/homes/sseyler/Simulations/adk/charmm/dims/c_o'
outputname = basedir + '/analysis/q1q2' + '/dims-co_nca_rmsd-' + struct       \
    + '.dat'

top_dims = basedir + '/struct/' + 'open_final.psf'
trj_dims = basedir + '/trj/raw' + '/'
dims_tr = []
for i in xrange(1, 201):
    name = 'dims%04i.dcd' % i
    dims_tr.append(trj_dims + name)
trjs = dims_tr
top = top_dims

prefix = '/nfs/homes/sseyler/Dropbox/Beckstein/data'
prot = '/AdK'
ofp = basedir + '/analysis'
ofn1 = ofp + '/' + 'q1-1ake_'
ofn2 = ofp + '/' + 'q2-4ake_'
ext1 = '.dat'
cutoff = 8.0 # in Angstroms

#------------------------------------------------
# Setup reference structures
#------------------------------------------------
reffpath = prefix + prot + '/boundary_conformations'
ref_c = Universe(reffpath+'/1AKE_A.pdb')
ref_o = Universe(reffpath+'/4AKE_A.pdb')
atomsel = 'name CA' # change to 'name CA and segid A' for top files with chains
ca_c = ref_c.selectAtoms(atomsel)
ca_o = ref_o.selectAtoms(atomsel)
# Reference coordinates
ref_o.atoms.translate(-ref_o.atoms.centerOfMass())
ref_o_coor =  ref_o.atoms.CA.coordinates()

# Process a single trajectory for q1q2 analysis
if False:
    u = Universe(top, basedir + 'adk_eq_langevin_'+struct+'.dcd')
    # u = Universe(filename, multiframe=True) # Excruciatingly slow LoL
    CA1 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_c, \
        radius=cutoff,outfile=ofn1+ext1)
    CA2 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_o, \
        radius=cutoff,outfile=ofn2+ext1)
    CA1.run(force=True)
    CA2.run(force=True)
    q1 = CA1.timeseries[1]
    q2 = CA2.timeseries[1]
    np.savetxt(outputname, np.vstack((q1,q2)).T, fmt="%7.5e")

if False:
    u = Universe(top, basedir+struct+'_a_emin.dcd')
    # u = Universe(filename, multiframe=True) # Excruciatingly slow LoL
    CA1 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_c, \
        radius=cutoff,outfile=ofn1+ext1)
    CA2 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_o, \
        radius=cutoff,outfile=ofn2+ext1)
    CA1.run(force=True)
    CA2.run(force=True)
    q1 = CA1.timeseries[1]
    q2 = CA2.timeseries[1]
    np.savetxt(outputname, np.vstack((q1,q2)).T, fmt="%7.5e")


# Process an ensemble of trajectories for q1q2 analysis
if True:
    import multiprocessing as mp

    # Get raw trajectories
    top_dims = basedir + '/struct/' + 'open_final.psf'
    trj_dims = basedir + '/trj/fit' + '/' + struct + '/'
    dims_tr = []
    for i in xrange(1, 201):
        name = 'dims%04i_fit-%s.dcd' % (i, struct)
        dims_tr.append(trj_dims + name)

    # with Timer(verbose=True) as t:

    def gen_ts(top, filenames, ts_q):
        pid = str(os.getpid())
        i = 0
        for filename in filenames:
            print filename

            fp, sep, fn = filename.rpartition("/")
            lhs, sep, rhs = fn.rpartition(".")
            u = Universe(top, filename)

            #-- NCA -----
            CA1 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_c, \
                radius=cutoff,outfile=ofn1+lhs+ext1)
            CA2 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_o, \
                radius=cutoff,outfile=ofn2+lhs+ext1)
            CA1.run(force=True)
            CA2.run(force=True)
            q1 = CA1.timeseries[1]
            q2 = CA2.timeseries[1]

            #-- RMSD -----
            rmsd_trj = np.zeros(len(u.trajectory))
            idx = 0
            u.trajectory.rewind()
            for ts in u.trajectory:
                mobile_coor =  u.atoms.CA.coordinates()
                rmsd_trj[idx] = rmsd(mobile_coor, ref_o_coor)
                idx += 1

            # Cutting out first ts from rmsd_trj b/c ContactAnalysis1
            # lops off a timestep somewhere...
            ts_q.put(np.vstack((q1,q2,rmsd_trj[1:])).T)

    ts_q = mp.Queue()
    ntrajs = len(trjs)
    nproc = mp.cpu_count()
    chunksize = ntrajs/nproc + 1 # integer division
    jobs = []
    s = e = 0
    for i in xrange(nproc):
        s = i*chunksize
        e = s + chunksize if e <= ntrajs else None
        p = mp.Process(target=gen_ts, args=(top, trjs[s:e], ts_q))
        jobs.append(p)
        p.start()

    temp = []
    for i in xrange(ntrajs):
        temp.append(ts_q.get())
    np.savetxt(outputname, np.vstack(tuple(temp)), fmt="%7.5e")

    for j in jobs:
        j.join()
    ts_q.close()