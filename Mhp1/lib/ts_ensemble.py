# from CodeTiming import Timer
import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.contacts
from MDAnalysis.analysis.align import *
import os

#------------------------------------------------
# Options
#------------------------------------------------
method = 'dims'
protein = 'dt'
direction = 'co'
N = 6

nca, rmsda, pla, angles, dmcd = False, True, False, False, False

pname_pf = 'dims' if method=='dims' else 'pathway'
ref_struct = '4ake' if protein=='adk' else '1ddt'
# ref_struct = '1ake' if protein=='adk' else '1mdt'

if protein == 'adk':
    struct_c, struct_o = '1ake_a_core.pdb', '4ake_a_core.pdb'
    name_c, name_o = '1ake', '4ake'
elif protein == 'dt':
    struct_c, struct_o = '1mdt_a_fit.pdb', '1ddt_a_fit.pdb'
    name_c, name_o = '1mdt', '1ddt'
if pla or angles or rmsda:
    pathname = 'fit/%s/%s' % (ref_struct, pname_pf)
    postfix = '_fit-%s' % ref_struct
else: pathname = pname_pf


def select_sims(method, protein, direction, name, pf=''):
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
                % (basedir, direction, name, i+1, pf) for i in xrange(N)]
    elif method == 'gp':
        basedir = '/nfs/homes/sseyler/Simulations/%s/%s' % (protein, method)
        if protein == 'adk':
            top_file = '%s/top/1ake.pdb' % basedir
        elif protein == 'dt':
            top_file = '%s/top/%s' % (basedir, '1mdt_a.pdb') #'1mdt_a.pdb'
        else:
            print "Must specify either \"adk\" or \"dt\" to select protein."
        tr_list = ["%s/trj/%s/%s%i%s.dcd"                                      \
                % (basedir, direction, name, i+1, pf) for i in xrange(N)]
    else:
        print "Must specify either \"dims\" or \"gp\" to select method."
    return basedir, top_file, tr_list


def get_refs_asUniverses(protein, struct_c, struct_o, pf=''):
    prefix = '/nfs/homes/sseyler/Simulations/%s/coords' % protein # workstation
    reffpath = prefix + pf
    fp_c = '%s/%s' % (reffpath, struct_c)
    fp_o = '%s/%s' % (reffpath, struct_o)
    return mda.Universe(fp_c), mda.Universe(fp_o)

def get_ref_coords(protein, struct_c, struct_o, selection='name CA'):
    u_ref_c, u_ref_o = get_refs_asUniverses(protein, struct_c, struct_o)
    ref_c_coor = u_ref_c.selectAtoms(selection).coordinates()
    ref_o_coor = u_ref_o.selectAtoms(selection).coordinates()
    return ref_c_coor, ref_o_coor

def rmsd_trajectory_fitting(method, protein, direction, pathname,              \
        struct_c, struct_o, selection='name CA', write_refs_tofile=False):
    import multiprocessing as mp

    #######################################################
    prefix = '/nfs/homes/sseyler/Simulations/%s' % protein # workstation
    reffpath = prefix + '/coords/boundary_conformations'
    fp_c = '%s/core_norm/%s' % (reffpath, struct_c)
    fp_o = '%s/core_norm/%s' % (reffpath, struct_o)
    ref_c = mda.Universe(fp_c)
    ref_o = mda.Universe(fp_o)

    refcoords_c, refcoords_o = get_ref_coords(protein, struct_c, struct_o)
    # atoms = 'name CA'
    # refcoords_c = ref_c.selectAtoms(atoms).coordinates() \
    #         - ref_c.selectAtoms(atoms).centerOfMass()
    # refcoords_o = ref_o.selectAtoms(atoms).coordinates() \
    #         - ref_o.selectAtoms(atoms).centerOfMass()

    R, rmsd_val = rotation_matrix(refcoords_c, refcoords_o)
    # refcoords_c_rot = np.asarray(refcoords_c[:,:] * (R.T))

    if write_refs_tofile:
        ref_o.atoms.translate(-ref_o.atoms.CA.centerOfMass())
        ref_c.atoms.translate(-ref_c.atoms.CA.centerOfMass())
        ref_c.atoms.rotate(R)
        ref_c.atoms.write(reffpath + "/core_norm/" + struct_c)
        ref_o.atoms.write(reffpath + "/core_norm/" + struct_o)
    #######################################################

    basedir, top_file, tr_list = select_sims(method, protein, direction, pathname)

    # Input trajectories
    trjs = tr_list
    top = top_file

    if method == 'dims':  charmm, implicit = '/charmm', '/implicit'
    else:                 charmm, implicit = '', ''
    outfp = '/nfs/homes/sseyler/Simulations/%s/%s/%s%s/trj/%s/fit/%s/' \
            % (protein, charmm, method, implicit, direction, ref_struct)

    # Get reference structures
    ref_c, ref_o = get_refs_asUniverses(protein, struct_c, struct_o)

    u_ref = ref_o if ref_struct=='4ake' else ref_c

    def rmsd_fit(top, filenames):
        for filename in filenames:
            fp, sep, fn = filename.rpartition("/")
            lhs, sep, rhs = fn.rpartition(".")
            newname = lhs + "_fit-" + ref_struct + sep + rhs
            outputfp = outfp + newname
            if os.path.isfile(outputfp):
                print 'The file %s already exists. Skipping...' % newname
            else:
                print "Writing file " + newname + " ..."
                u = mda.Universe(top, filename)
                rms_fit_trj(u, u_ref, select='name CA', filename=outfp + newname)

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


#################################################
# Setup NCA
#################################################
if nca:
    from MDAnalysis.analysis.contacts import ContactAnalysis1

    basedir, top_file, tr_list =                                                   \
            select_sims(method, protein, direction, pathname, pf=postfix)
    outputname = '%s/data/nca/%s-%s-%s_nca-rmsd-%s.dat'                            \
            % (basedir, method, protein, direction, ref_struct)

    ofp = basedir + '/data'
    ofn1 = '%s/q1-%s_' % (ofp, name_c)
    ofn2 = '%s/q2-%s_' % (ofp, name_o)
    ext1 = '.dat'
    cutoff = 8.0 # in Angstroms

#------------------------------------------------
# Setup reference structures
#------------------------------------------------
ref_c, ref_o = get_refs_asUniverses(protein, struct_c, struct_o, pf='/core_norm')
atomsel = 'name CA' # change to 'name CA and segid A' for top files with chains

ca_c = ref_c.selectAtoms(atomsel)
ca_o = ref_o.selectAtoms(atomsel)
ref_c.atoms.translate(-ca_c.centerOfMass())
ref_o.atoms.translate(-ca_o.centerOfMass())
ref_c_coor_temp = ca_c.coordinates()
ref_o_coor = ca_o.coordinates()

# ref_c_coor_temp, ref_o_coor = get_ref_coords(protein, struct_c, struct_o,      \
#         selection='name CA')
# ref_c_coor_temp = ref_c.atoms.CA.coordinates() - ref_c.atoms.CA.centerOfMass()
# ref_o_coor = ref_o.atoms.CA.coordinates() - ref_o.atoms.CA.centerOfMass()

R, rmsd_val = rotation_matrix(ref_c_coor_temp, ref_o_coor)
ref_c_coor = np.asarray(ref_c_coor_temp * (R.T))
# print rmsd_val
# print ref_c.atoms.centerOfMass()
# print ref_o.atoms.centerOfMass()
# print R

# Process a single trajectory for q1q2 analysis
if False: # equilibrium langevin dynamics
    u = mda.Universe(top, basedir + 'adk_eq_langevin_'+adk+'.dcd')
    # u = mda.Universe(filename, multiframe=True) # Excruciatingly slow LoL
    CA1 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_c, \
        radius=cutoff,outfile=ofn1+ext1)
    CA2 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_o, \
        radius=cutoff,outfile=ofn2+ext1)
    CA1.run(force=True)
    CA2.run(force=True)
    q1 = CA1.timeseries[1]
    q2 = CA2.timeseries[1]
    np.savetxt(outputname, np.vstack((q1,q2)).T, fmt="%7.5e")

if False: # energy minimization
    u = mda.Universe(top, basedir+adk+'_a_emin.dcd')
    # u = mda.Universe(filename, multiframe=True) # Excruciatingly slow LoL
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

    # NMP: centers of geometry of the backbone and C-beta atoms in residues 115-125
    # (CORE-LID), 90-100 (CORE), and 35-55 (NMP)
    res_nmp_hinge = 'resid 90-100'
    res_nmp_in = 'resid 35-55'
    res_nmp_out = 'resid 115-125'

    # LID: 179-185 (CORE), 115-125 (CORE-hinge-LID), and 125-153 (LID)
    res_lid_hinge = res_nmp_out
    res_lid_in = 'resid 125-153'
    res_lid_out = 'resid 179-185'

    res_lid = "resid 122-159"
    res_core = "resid 1-29 or resid 60-121 or resid 160-214"
    res_nmp = "resid 30-59"

    def vsqnorm(v, axis=0):
        return np.sum(v*v, axis=axis)

    toDegrees = 180.0/np.pi
    def computeangle(v1, v2):
        return np.arccos( np.dot(v1, v2)/((vsqnorm(v1)*vsqnorm(v2))**0.5) )*toDegrees


    def gen_ts(top, filenames, ts_q):
        pid = str(os.getpid())
        i = 0
        for filename in filenames:
            print filename

            fp, sep, fn = filename.rpartition("/")
            lhs, sep, rhs = fn.rpartition(".")
            u = mda.Universe(top, filename)

            #-- NCA ---------
            if nca:
                CA1 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_c, \
                    radius=cutoff,outfile=ofn1+lhs+ext1)
                CA2 = ContactAnalysis1(u, selection=atomsel,refgroup=ca_o, \
                    radius=cutoff,outfile=ofn2+lhs+ext1)
                CA1.run(force=True)
                CA2.run(force=True)
                q1 = CA1.timeseries[1]
                q2 = CA2.timeseries[1]

            #-- DXDY --------
            if pla:
                dxdy_trj = np.zeros((len(u.trajectory), 2))

                C, O = ref_c_coor, ref_o_coor
                OC = np.ravel(C - O)
                isq_natoms = (3./len(OC))**0.5
                uOC = OC/vsqnorm(OC)**0.5

                idx = 0
                u.trajectory.rewind()
                for ts in u.trajectory:
                    P = u.atoms.CA.coordinates()
                    OP = np.ravel(P - O)
                    dx = np.dot(OP, uOC)
                    OM = dx*uOC
                    dy = vsqnorm((OP - OM))**0.5
                    dxdy_trj[idx,:] = np.array([dx, dy])*isq_natoms
                    idx += 1
                    # or can get dy vector and divide by dy unit vector
                    # computed by rotating in xy-plane
                ts_q.put(dxdy_trj)


            #-- RMSD --------
            if rmsda:
                rmsd_trj = np.zeros(len(u.trajectory))
                idx = 0
                u.trajectory.rewind()
                for ts in u.trajectory:
                    mobile_coor = u.atoms.CA.coordinates() - u.atoms.CA.centerOfMass()
                    rmsd_trj[idx] = rmsd(mobile_coor, ref_o_coor)
                    idx += 1
                ts_q.put(rmsd_trj)

            #-- theta -------
            if angles:
                aa_trj = np.zeros((len(u.trajectory), 2))
                # theta_nmp = theta_lid = np.zeros(len(u.trajectory))
                bb = u.selectAtoms("backbone")

                atoms_nmp_out = bb.selectAtoms(res_nmp_out)
                atoms_nmp_hinge = bb.selectAtoms(res_nmp_hinge)
                atoms_nmp_in = bb.selectAtoms(res_nmp_in)
                atoms_lid_hinge = bb.selectAtoms(res_lid_hinge)
                atoms_lid_out = bb.selectAtoms(res_lid_out)
                atoms_lid_in = bb.selectAtoms(res_lid_in)
                idx = 0
                u.trajectory.rewind()
                for ts in u.trajectory:

                    nmp_hinge_centroid = atoms_nmp_hinge.centroid()
                    v1 = atoms_nmp_out.centroid() - nmp_hinge_centroid
                    v2 = atoms_nmp_in.centroid() - nmp_hinge_centroid
                    aa_trj[idx,0] = computeangle(v1, v2)

                    lid_hinge_centroid = atoms_lid_hinge.centroid()
                    v1 = atoms_lid_in.centroid() - lid_hinge_centroid
                    v2 = atoms_lid_out.centroid() - lid_hinge_centroid
                    aa_trj[idx,1] = computeangle(v1, v2)
                    idx += 1
                ts_q.put(aa_trj)

            #-- DMCD --------
            if dmcd:
                ca = u.selectAtoms("name CA")
                mobile_lid = ca.selectAtoms(res_lid)
                mobile_core = ca.selectAtoms(res_core)
                mobile_nmp = ca.selectAtoms(res_nmp)
                dmcd_nmp = dmcd_lid = np.zeros(len(u.trajectory))
                idx = 0
                for ts in u.trajectory:
                    lid_com = mobile_lid.centerOfMass()
                    core_com = mobile_core.centerOfMass()
                    nmp_com = mobile_nmp.centerOfMass()
                    dmcd_nmp[idx] = vsqnorm(nmp_com-core_com)
                    dmcd_lid[idx] = vsqnorm(lid_com-core_com)
                    idx += 1
                dmcd_nmp = dmcd_nmp**0.5
                dmcd_lid = dmcd_lid**0.5

            # Cutting out first ts from rmsd_trj b/c ContactAnalysis1
            # lops off a timestep somewhere...
            # ts_q.put(np.vstack((q1,q2,rmsd_trj[1:],theta_nmp[1:],             \
            #     theta_lid[1:],dmcd_nmp[1:],dmcd_lid[1:])).T)
            if nca and rmsda:
                ts_q.put(np.vstack((q1,q2,rmsd_trj[1:])).T)


    # direction = direction + '/fit/1ddt'
    basedir, top_file, tr_list =                                               \
        select_sims(method, protein, direction, pathname, pf=postfix)
    datadir = '/nfs/homes/sseyler/Projects/PathwayAnalysis/Phase1'

    if pla:
        name = 'pla'
    elif angles:
        name = 'aa'
    elif nca:
        name = 'nca'
    elif rmsda:
        name = 'rmsd'

    fname = '%s_%s-%s-%s-%s.dat' % (name, method, protein, direction, ref_struct)
    outputname = '%s/data/%s/%s' % (basedir, name, fname)
    if nca:
        outputname2 = '%s/data/ensemble/%s/%s/%s/%s/%s'                   \
            % (datadir, name, method, protein, direction, fname)
    else:
        outputname2 = '%s/data/ensemble/%s/%s' % (datadir, name, fname)


    top = top_file
    trjs = tr_list

    ts_q = mp.Queue()
    ntraj = len(trjs)
    nproc = mp.cpu_count()
    chunksize = 1.0/nproc*ntraj
    jobs = []
    for i in xrange(nproc):
        s = int(round(i*chunksize))
        e = int(round((i+1)*chunksize))
        p = mp.Process(target=gen_ts, args=(top, trjs[s:e], ts_q))
        jobs.append(p)
        p.start()

    temp = []
    for i in xrange(ntraj):
        temp_trj = ts_q.get()
        temp.append(temp_trj)

        fname = '%s_%s-%s-%s-%s_%i.dat' % (name, method, protein, direction, ref_struct, i)
        outputname = '%s/data/%s/%s' % (basedir, name, fname)
        outputname2 = '%s/data/ensemble/%s/%s' % (datadir, name, fname)
        np.savetxt(outputname, temp_trj, fmt="%7.5e")
        np.savetxt(outputname2, temp_trj, fmt="%7.5e")
    datadir = None
    if nca and rmsda:
        np.savetxt(outputname, np.vstack(tuple(temp)), fmt="%7.5e")

    print 'Wrote data to: \n  --> %s' % outputname
    if datadir is not None:
        np.savetxt(outputname2, np.vstack(tuple(temp)), fmt="%7.5e")
        print 'Wrote data to: \n  --> %s' % outputname2

    for j in jobs:
        j.join()
    ts_q.close()