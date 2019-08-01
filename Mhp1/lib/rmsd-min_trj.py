# from CodeTiming import Timer
import numpy as np
from MDAnalysis import Universe
import MDAnalysis.analysis.contacts
from MDAnalysis.analysis.align import *
import MDAnalysis as mda
import os

#################################################
# Minimize RMSD DIMS trajs w.r.t. reference
#################################################

def select_sims(method, protein, direction, pname, N, pf=''):
    if method == 'dims':
        basedir = '/nfs/homes/sseyler/Simulations/%s/charmm/%s/implicit' \
                % (protein, method)
        if protein == 'adk':
            top_file = '%s/top/%s' % (basedir, 'adk4ake.psf')
        elif protein == 'dt':
            top_file = '%s/top/%s' % (basedir, '1mdt_a_core.pdb') #
        else:
            print "Must specify either \"adk\" or \"dt\" to select protein."
        tr_list = ['%s/trj/%s/%s%04i%s.dcd'                                    \
                % (basedir, direction, pname, i+1, pf) for i in xrange(N)]
    elif method == 'gp':
        basedir = '/nfs/homes/sseyler/Simulations/%s/%s' % (protein, method)
        if protein == 'adk':
            top_file = '%s/top/1ake.pdb' % basedir
        elif protein == 'dt':
            top_file = '%s/top/%s' % (basedir, '1mdt_a.pdb') #'1mdt_a.pdb'
        else:
            print "Must specify either \"adk\" or \"dt\" to select protein."
        tr_list = ["%s/trj/%s/%s%i%s.dcd"                                      \
                % (basedir, direction, pname, i+1, pf) for i in xrange(N)]
    else:
        print "Must specify either \"dims\" or \"gp\" to select method."
    return basedir, top_file, tr_list


# def setup_ref_structs(protein, struct_c, struct_o):
#     prefix = '/nfs/homes/sseyler/Simulations/%s' % protein # workstation
#     reffpath = prefix + '/coords'
#     fp_c = '%s/core_norm/%s' % (reffpath, struct_c)
#     fp_o = '%s/core_norm/%s' % (reffpath, struct_o)
#     return mda.Universe(fp_c), mda.Universe(fp_o)

def setup_ref_structs(protein, struct_c, struct_o, pf='/core_norm'):
    #------------------------------------------------
    # Setup reference structures
    #------------------------------------------------
    ref_c, ref_o = get_refs_asUniverses(protein, struct_c, struct_o, pf='/core_norm')
    atomsel = 'name CA' # change to 'name CA and segid A' for top files with chains
    print type(ref_c), type(ref_o)
    ca_c = ref_c.selectAtoms(atomsel)
    ca_o = ref_o.selectAtoms(atomsel)
    ref_c.atoms.translate(-ca_c.centerOfMass())
    ref_o.atoms.translate(-ca_o.centerOfMass())
    ref_c_coor_temp = ca_c.coordinates()
    ref_o_coor = ca_o.coordinates()

    R, rmsd_val = rotation_matrix(ref_c_coor_temp, ref_o_coor)
    ref_c_coor = np.asarray(ref_c_coor_temp * (R.T))
    return ref_c_coor, ref_o_coor

def get_refs_asUniverses(protein, struct_c, struct_o, pf='/core_norm'):
    prefix = '/nfs/homes/sseyler/Simulations/%s/coords' % protein # workstation
    reffpath = prefix + pf
    fp_c = '%s/%s' % (reffpath, struct_c)
    fp_o = '%s/%s' % (reffpath, struct_o)
    return mda.Universe(fp_c), mda.Universe(fp_o)

def rmsd_trajectory_fitting(params, N, adk=('',''), ref_struct='4ake', \
        pname='', write_refs_tofile=False, force=False):
    import multiprocessing as mp

    method, protein, direction = params
    struct_c, struct_o = adk
    #######################################################
    prefix = '/nfs/homes/sseyler/Simulations/%s' % protein # workstation
    reffpath = prefix + '/coords'
    fp_c = '%s/core_norm/%s' % (reffpath, struct_c)
    fp_o = '%s/core_norm/%s' % (reffpath, struct_o)
    ref_c = mda.Universe(fp_c)
    ref_o = mda.Universe(fp_o)

    atoms = 'name CA'
    refcoords_c = ref_c.selectAtoms(atoms).coordinates() \
            - ref_c.selectAtoms(atoms).centerOfMass()
    refcoords_o = ref_o.selectAtoms(atoms).coordinates() \
            - ref_o.selectAtoms(atoms).centerOfMass()

    R, rmsd_val = rotation_matrix(refcoords_c, refcoords_o)
    # refcoords_c_rot = np.asarray(refcoords_c[:,:] * (R.T))

    #######################################################
    # Get reference structures
    ref_c, ref_o = get_refs_asUniverses(protein, struct_c, struct_o)
    u_ref = ref_o if ref_struct=='4ake' else ref_c

    if write_refs_tofile:
        ref_o.atoms.translate(-ref_o.atoms.CA.centerOfMass())
        ref_c.atoms.translate(-ref_c.atoms.CA.centerOfMass())
        ref_c.atoms.rotate(R)
        ref_c.atoms.write(fp_c)
        ref_o.atoms.write(fp_o)
    #######################################################

    basedir, top_file, tr_list =                                               \
            select_sims(method, protein, direction, pname, N)

    # Input trajectories
    trjs = tr_list
    top = top_file

    if method == 'dims':  charmm, implicit = '/charmm', '/implicit'
    else:                 charmm, implicit = '', ''
    outfp = '/nfs/homes/sseyler/Simulations/%s/%s/%s%s/trj/%s/fit/%s/' \
            % (protein, charmm, method, implicit, direction, ref_struct)

    # Get reference structures
    # ref_c, ref_o = setup_ref_structs(protein, struct_c, struct_o)

    def rmsd_fit(top, filenames):
        for filename in filenames:
            fp, sep, fn = filename.rpartition("/")
            lhs, sep, rhs = fn.rpartition(".")
            newname = lhs + "_fit-" + ref_struct + sep + rhs
            outputfp = outfp + newname
            if os.path.isfile(outputfp) and not force:
                print 'The file %s already exists. Skipping...' % newname
            else:
                print "Writing file " + newname + " ..."
                u = Universe(top, filename)
                rms_fit_trj(u, ref_o, select='name CA', filename=outfp + newname)

    # with Timer(verbose=True) as t:
    ntraj = len(trjs)
    nproc = mp.cpu_count()
    chunksize = 1.0/nproc*ntraj
    jobs = []
    for i in xrange(nproc):
        s = int(round(i*chunksize))
        e = int(round((i+1)*chunksize))
        p = mp.Process(target=rmsd_fit, args=(top, trjs[s:e]))
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()
#################################################