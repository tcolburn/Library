from MDAnalysis.analysis.align import *
import MDAnalysis as mda
import multiprocessing as mp


def rmsd_trajectory_fitting(prot, method, ntrj, ref, sel='name CA', postf='core'):

    if method == 'dims':
        basedir = workstation + '/Simulations/%s/charmm/dims/implicit' % prot
        top = basedir + '/top/adk4ake.psf'
        trjpref = 'dims'
    elif method == 'gp':
        basedir = workstation + '/Simulations/%s/gp' % prot
        top = basedir + '/top/1ake.pdb'
        trjpref = 'pathway'
    trjdir = basedir + '/trj/co/'
    
    # outfp = '/nfs/homes/sseyler/Simulations/%s/%s/%s%s/trj/%s/fit/%s/' \
    #          % (protein, charmm, method, implicit, direction, ref_struct)

    outfp = '%s/trj/co/fit/%s/' % (basedir, postf)

    trjs = []
    for i in xrange(1, ntrj+1):
        if method == 'dims':
            name = '%s%04i.dcd' % (trjpref, i)
        elif method == 'gp':
            name = '%s%i.dcd' % (trjpref, i)
        trjs.append(trjdir + name)

    def rmsd_fit(top, filenames):
        for filename in filenames:
            fp, sep, fn = filename.rpartition("/")
            lhs, sep, rhs = fn.rpartition(".")
            newname = lhs + "_fit-" + postf + sep + rhs
            outputfp = outfp + newname
            if os.path.isfile(outputfp):
                print 'The file %s already exists. Skipping...' % newname
            else:
                print "Writing file " + newname + " ..."
                u = mda.Universe(top, filename)
                rms_fit_trj(u, ref, select=sel, filename=outfp + newname)

    # with Timer(verbose=True) as t:
    ntrj = len(trjs)
    nproc = mp.cpu_count()
    chunksize = 1.0/nproc*ntrj
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


if __name__ == '__main__':

    prot = 'adk'
    co = True
    oc = False

    home = '/home/sseyler/Remote/spudda'
    office = '/nfs/homes/sseyler'
    workstation = office
    if 'Remote' in workstation:
        print 'Using home workstation...'
    else:
        print 'Using office workstation...'

    pf = workstation + '/Simulations/' # workstation

    # Get reference structure
    reffpath = '%s/%s/coords/core_norm' % (pf, prot)

    reffpath1AKE = '%s/1ake_a_core.pdb' % reffpath
    reffpath4AKE = '%s/4ake_a_core.pdb' % reffpath
    reffpath1AKE_ca = '%s/1ake_a_ca_core.pdb' % reffpath
    reffpath4AKE_ca = '%s/4ake_a_ca_core.pdb' % reffpath

    ref1AKE = mda.Universe(reffpath1AKE)
    ref4AKE = mda.Universe(reffpath4AKE)

    ref1AKE_ca = ref1AKE.selectAtoms('name CA')
    ref4AKE_ca = ref4AKE.selectAtoms('name CA')

    COREstr = "(resid 1:29 or resid 60:121 or resid 160:214)"
    ref1akeCORE_ca = ref1AKE_ca.selectAtoms(COREstr).coordinates() \
            - ref1AKE_ca.selectAtoms(COREstr).centerOfMass()
    ref4akeCORE_ca = ref4AKE_ca.selectAtoms(COREstr).coordinates() \
            - ref4AKE_ca.selectAtoms(COREstr).centerOfMass()
    ref_coords = 0.5*(ref1akeCORE_ca + ref4akeCORE_ca)

    selection = COREstr + " and name CA"

    u_ref = mda.Universe(reffpath1AKE)
    u_ref_ca = mda.Universe(reffpath1AKE_ca)

    ref_ca_core = u_ref_ca.selectAtoms(COREstr)

    # print 'coords before: ', len(ref_ca_core.coordinates()), ref_ca_core.coordinates()
    # print 'all coords before: ', len(u_ref_ca.atoms.coordinates()), u_ref_ca.atoms.coordinates()
    ref_ca_core.set_positions(ref_coords)
    # print 'coords after: ', len(ref_ca_core.coordinates()), ref_ca_core.coordinates()
    # print 'all coords after: ', len(u_ref_ca.atoms.coordinates()), u_ref_ca.atoms.coordinates()

    rmsd_trajectory_fitting('adk','gp',200, u_ref_ca, sel=selection)