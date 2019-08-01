import os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import matplotlib.patches as mpatches

import prettyplotlib as ppl
import MDAnalysis as mda

mpl.rc('font',**{'family':'serif','serif':['Helvetica']})
mpl.rc('text', usetex=False)

# Wrap around code using "with suppress_stdout()" to suppress output
from contextlib import contextmanager
@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

#------------------------------------------------
# Helper functions
#------------------------------------------------
def vsqnorm(v, axis=0):
        return np.sum(v*v, axis=axis)

deg_per_rad = 180/(np.pi)
def compute_angle(v1, v2):
    return np.arccos( np.dot(v1, v2)/((vsqnorm(v1)*vsqnorm(v2))**0.5) )*deg_per_rad

def select_sims(method, protein, direction, name, N, pf=''):
    if method == 'dims':
        bd = '/nfs/homes/sseyler/Simulations/%s/charmm/%s/implicit' \
                % (protein, method)
        if protein == 'adk':
            top_dims = '%s/top/%s' % (bd, 'adk4ake.psf')
        elif protein == 'dt':
            top_dims = '%s/top/%s' % (bd, '1mdt_a_core.pdb') #
        else:
            print "Must specify either \"adk\" or \"dt\" to select protein."
        tr_list = ['%s/trj/%s/%s%04i%s.dcd' % (bd, direction, name, i+1, pf) for i in xrange(N)]
    elif method == 'gp':
        bd = '/nfs/homes/sseyler/Simulations/%s/%s' % (protein, method)
        if protein == 'adk':
            top_gp = '%s/top/1ake.pdb' % bd
        elif protein == 'dt':
            top_gp = '%s/top/%s' % (bd, '1mdt_a.pdb') #'1mdt_a.pdb'
        else:
            print "Must specify either \"adk\" or \"dt\" to select protein."
        tr_list = ["%s/trj/%s/%s%i%s.dcd" % (bd, direction, name, i+1, pf) for i in xrange(N)]
    else:
        print "Must specify either \"dims\" or \"gp\" to select method."
    return bd, top_file, tr_list


def get_adk_angle_defs():
    # NMP: centers of geometry of the backbone and C-beta atoms in residues 115-125
    # (CORE-LID), 90-100 (CORE), and 35-55 (NMP)
    res_nmp_hinge = 'resid 90-100'
    res_nmp_in = 'resid 35-55'
    res_nmp_out = 'resid 115-125'
    nmp_angle_def = res_nmp_hinge, res_nmp_in, res_nmp_out

    # LID: 179-185 (CORE), 115-125 (CORE-hinge-LID), and 125-153 (LID)
    res_lid_hinge = res_nmp_out
    res_lid_in = 'resid 125-153'
    res_lid_out = 'resid 179-185'
    lid_angle_def = res_lid_hinge, res_lid_in, res_lid_out

    return nmp_angle_def, lid_angle_def

def compute_ref_angles(ref_c, ref_o, angle_defs, of_refnames, sel='backbone'):
    '''
        Get LID-CORE and NMP-CORE angles for 1ake and 4ake
        and store as data files for future retrieval.
    '''
    c_bb = ref_c.selectAtoms(sel)
    o_bb = ref_o.selectAtoms(sel)
    nmp_angle_def, lid_angle_def = angle_defs

    #=== Compute NMP angles ========================
    res_nmp_hinge, res_nmp_in, res_nmp_out = nmp_angle_def
    # 1ake
    c_nmp_hinge = c_bb.selectAtoms(res_nmp_hinge)
    c_nmp_out = c_bb.selectAtoms(res_nmp_out)
    c_nmp_in = c_bb.selectAtoms(res_nmp_in)
    a1 = c_nmp_out.centroid() - c_nmp_hinge.centroid()
    a2 = c_nmp_in.centroid() - c_nmp_hinge.centroid()
    c_theta_nmp = compute_angle(a1, a2)
    # 4ake
    o_nmp_hinge = o_bb.selectAtoms(res_nmp_hinge)
    o_nmp_out = o_bb.selectAtoms(res_nmp_out)
    o_nmp_in = o_bb.selectAtoms(res_nmp_in)
    a1 = o_nmp_out.centroid() - o_nmp_hinge.centroid()
    a2 = o_nmp_in.centroid() - o_nmp_hinge.centroid()
    o_theta_nmp = compute_angle(a1, a2)

    #=== Compute LID angles ========================
    res_lid_hinge, res_lid_in, res_lid_out = lid_angle_def
    # 1ake
    c_lid_hinge = c_bb.selectAtoms(res_lid_hinge)
    c_lid_out = c_bb.selectAtoms(res_lid_out)
    c_lid_in = c_bb.selectAtoms(res_lid_in)
    a1 = c_lid_in.centroid() - c_lid_hinge.centroid()
    a2 = c_lid_out.centroid() - c_lid_hinge.centroid()
    c_theta_lid = compute_angle(a1, a2)
    # 4ake
    o_lid_hinge = o_bb.selectAtoms(res_lid_hinge)
    o_lid_out = o_bb.selectAtoms(res_lid_out)
    o_lid_in = o_bb.selectAtoms(res_lid_in)
    a1 = o_lid_in.centroid() - o_lid_hinge.centroid()
    a2 = o_lid_out.centroid() - o_lid_hinge.centroid()
    o_theta_lid = compute_angle(a1, a2)

    ref_angles_nmp = np.array([c_theta_nmp, o_theta_nmp])
    ref_angles_lid = np.array([c_theta_lid, o_theta_lid])
    if sel == 'name CA':
        postf = '_ca'
    elif sel == 'backbone':
        postf = '_bb'
    np.savetxt(of_refnames[0]+postf+'.dat', ref_angles_nmp, fmt="%7.5f")
    np.savetxt(of_refnames[1]+postf+'.dat', ref_angles_lid, fmt="%7.5f")

def gen_angle_trajectories(analysis_dir, params, N,                            \
        runs=None, adk=('',''), sel='backbone', pathname=None, force=False):
    ###############################################################################
    # AA - angle-angle coordinates
    ###############################################################################
    methods, protein, directions = params
    adk_c, adk_o = adk

    tr_list_co, tr_list_oc = [], []
    for m in methods:
        method, postfix = m[0], m[2]
        for d in directions:
            if method == 'dims':
                bd = '/nfs/homes/sseyler/Simulations/%s' % protein
                dimsdir = '%s/charmm/%s/implicit' % (bd, method)
                top_dir = '%s/top' % dimsdir
                top_dims = '%s/%s' % (top_dir, 'adk4ake.psf')
                tl = ['%s/trj/%s/dims%04i.dcd' % (dimsdir, d, i+1) for i in xrange(N)]
            elif method == 'tmdslow':
                bd = '/nfs/homes/sseyler/Simulations/%s' % protein
                tmddir = '%s/charmm/tmd/implicit' % (bd)
                top_dir = '%s/top' % tmddir
                top_tmd = '%s/%s' % (top_dir, 'adk_closed.psf')
                tl = ['%s/trj/%s/slow/tmd%04i.dcd' % (tmddir, d, i+1) for i in xrange(N)]
            elif method == 'tmdfast':
                bd = '/nfs/homes/sseyler/Simulations/%s' % protein
                tmddir = '%s/charmm/tmd/implicit' % (bd)
                top_dir = '%s/top' % tmddir
                top_tmd = '%s/%s' % (top_dir, 'adk_closed.psf')
                tl = ['%s/trj/%s/fast/tmd%04i.dcd' % (tmddir, d, i+1) for i in xrange(N)]
            elif method == 'gp':
                bd = '/nfs/homes/sseyler/Simulations/%s' % protein
                gpdir = '%s/%s' % (bd, method)
                top_dir = '%s/top' % gpdir
                top_gp = '%s/%s' % (top_dir, '1ake.pdb')
                tl = ["%s/trj/%s/pathway%i.dcd" % (gpdir, d, i+1) for i in xrange(N)]
            else:
                bd = '/nfs/homes/sseyler/Simulations/%s' % protein
                top_dir = '%s/coords/core_norm' % bd
                top_file = '%s/%s' % (top_dir, adk_c)
                if method != 'intp':
                    if postfix:
                        tl = ['%s/public_servers/%s/%s/%s/%s_%s.pdb'                 \
                                % (bd, method, d, run, pathname, postfix) for run in runs]
                    else:
                        tl = ['%s/public_servers/%s/%s/%s/%s.pdb'                     \
                            % (bd, method, d, run, pathname) for run in runs]
                else:
                    if d == 'co':
                        if sel == 'name CA':
                            tl = ['%s/public_servers/intp/001/interp_adk_co_ca.pdb' % bd]
                        else:
                            tl = ['%s/public_servers/intp/001/interp_adk_co_full.pdb' % bd]
            if d == 'co': tr_list_co.extend(tl)
            elif d == 'oc': tr_list_oc.extend(tl)

    refdir = '/nfs/homes/sseyler/Simulations/adk/coords'
    ifname_top_c = '%s/core_norm/%s' % (refdir, adk_c)
    ifname_top_o = '%s/core_norm/%s' % (refdir, adk_o)
    ref_c = mda.Universe(ifname_top_c)
    ref_o = mda.Universe(ifname_top_o)
    of_refnames = (analysis_dir+'/adk-nmp_ref', analysis_dir+'/adk-lid_ref')

    angle_defs = get_adk_angle_defs()
    nmp_angle_def, lid_angle_def = angle_defs
    res_nmp_hinge, res_nmp_in, res_nmp_out = nmp_angle_def
    res_lid_hinge, res_lid_in, res_lid_out = lid_angle_def
    compute_ref_angles(ref_c, ref_o, angle_defs, of_refnames, sel)

    Nm, Nd = len(methods), len(directions)
    for i, m in zip(xrange(Nm), methods):
        method, postfix = m[0], m[2]
        for j, d in zip(xrange(Nd), directions):
            print 'CV calculation: angle-angle for %s|%s|%s, run: '                \
                    % (method, protein, d)
            if postfix != None:
                basename = '%s/%s_%s-%s-%s-%s'                                     \
                        % (analysis_dir, 'aa', method, postfix, d, protein)
            else:
                basename = '%s/%s_%s-%s-%s'                                        \
                        % (analysis_dir, 'aa', method, d, protein)
            for run in xrange(N):
                ofname1 = '%s-%s_%03i.dat' % (basename, 'nmp', run+1)
                ofname2 = '%s-%s_%03i.dat' % (basename, 'lid', run+1)
                # Skip iteration if angle-angle trajectory has already been computed
                if os.path.isfile(ofname1) and force==False:
                    fname = ofname1.split('/')[-1]
                    print '  --> %i (AA trajectory file exists, skipping...)' % run
                    pass
                else:
                    idx = i*N + run
                    if d == 'co': trj = tr_list_co[idx]
                    elif d == 'oc': trj = tr_list_oc[idx]

                    print '  --> %i' % (run+1)
                    if method == 'dims':
                        u = mda.Universe(top_dims, trj)
                    elif 'tmd' in method:
                        u = mda.Universe(top_tmd, trj, format="LAMMPS")
                    elif method == 'gp':
                        u = mda.Universe(top_gp, trj)
                    else:
                        u = mda.Universe(trj, multiframe=True)
                    bb = u.selectAtoms(sel)

                    # Theta NMP: 115-125 (CORE-hinge-LID), 90-100 (CORE), and 35-55 (NMP)
                    atoms_nmp_out = bb.selectAtoms(res_nmp_out)
                    atoms_nmp_hinge = bb.selectAtoms(res_nmp_hinge)
                    atoms_nmp_in = bb.selectAtoms(res_nmp_in)

                    theta_nmp_trj = np.zeros(len(u.trajectory))
                    idx = 0
                    for ts in u.trajectory:
                        nmp_hinge_centroid = atoms_nmp_hinge.centroid()
                        v1 = atoms_nmp_out.centroid() - nmp_hinge_centroid
                        v2 = atoms_nmp_in.centroid() - nmp_hinge_centroid
                        theta_nmp_trj[idx] = compute_angle(v1, v2)
                        idx += 1

                    # Theta LID: 179-185 (CORE), 115-125 (CORE-hinge-LID), and 125-153 (LID)
                    atoms_lid_hinge = bb.selectAtoms(res_lid_hinge)
                    atoms_lid_out = bb.selectAtoms(res_lid_out)
                    atoms_lid_in = bb.selectAtoms(res_lid_in)

                    theta_lid_trj = np.zeros(len(u.trajectory))
                    idx = 0
                    for ts in u.trajectory:
                        lid_hinge_centroid = atoms_lid_hinge.centroid()
                        v1 = atoms_lid_in.centroid() - lid_hinge_centroid
                        v2 = atoms_lid_out.centroid() - lid_hinge_centroid
                        theta_lid_trj[idx] = compute_angle(v1, v2)
                        idx += 1
                    ofname1 = '%s-%s_%s' % (basename, 'nmp', runs[run])
                    ofname2 = '%s-%s_%s' % (basename, 'lid', runs[run])
                    # angles = np.transpose(np.vstack((theta_nmp_trj, theta_lid_trj)))
                    if method == 'intp':
                        if sel == 'name CA':
                            postf = '_ca'
                        else:
                            postf = '_full'
                        ofname1 = ofname1 + postf
                        ofname2 = ofname2 + postf
                    np.savetxt(ofname1+'.dat', theta_nmp_trj, fmt="%7.5f")
                    np.savetxt(ofname2+'.dat', theta_lid_trj, fmt="%7.5f")
                if (method == 'intp'): break # only one interp line traj

    return basename, runs, of_refnames


def plot_aa(analysis_dir, params, runs, refnames, plot_type, sel='backbone',   \
        pltparams=(None, None, (5.75, 12)), plotname=None, showplot=True):
    #-----------------------
    # Timeseries angles for 1ake and 4ake trajectories
    #-----------------------
    methods, protein, directions = params
    drmsd, angles, dmcd, idf = False, True, False, False
    st = en = sk = None

    N = len(runs)
    ofname_ref_nmp, ofname_ref_lid = refnames
    timeseries, twod = plot_type

    xlims, ylims, (fsize1, fsize2), lsize = pltparams
    #-- plot settings ------
    xmin, xmax = xlims
    ymin, ymax = ylims
    xdatarange, ydatarange = (xmax-xmin), (ymax-ymin)

    #-----------------------
    # 2D angles for 1ake and 4ake trajectories
    #-----------------------
    if twod:
        print "Plotting angles for 1ake and 4ake trajectories..."

        # Plot colors:
        import brewer2mpl
        import seaborn as sns
        sns.set_style('whitegrid')
        # Save a nice dark grey as a variable
        almost_black = '#262626'
        # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
        dark2 = brewer2mpl.get_map('Dark2', 'qualitative', 8).mpl_colors
        accent = brewer2mpl.get_map('Accent', 'qualitative', 8).mpl_colors
        set1 = brewer2mpl.get_map('Set1', 'qualitative', 9).mpl_colors
        s2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
        s3 = brewer2mpl.get_map('Set3', 'qualitative', 12).mpl_colors
        # colors = [s3[9],s3[5],s2[4],s3[3],s3[7],s2[5],s2[7],s2[6],s2[0],s3[10],s2[1]]
        # colors = [s3[0],s3[11]] + s3[2:8] + s3[9:11] + [s2[6]]
        # colors = [s3[1],s2[6]]+s3[2:8] + [s2[7]] + s3[9:11] +[s3[8]]
        # colors = [s3[11],s2[6],s3[2],s3[3],s2[0],s3[5],s3[7],s3[9],s3[4],s2[4],'gray']
        # colors = [s2[0],s3[11],s3[2],s3[3],s3[4],s3[5],s2[4],s3[7],s2[6],s3[9],'gray']
        # colors = [s2[0],s3[11],s2[3],s3[3],s2[2],s3[7],s2[4],s3[5],s3[9],s2[6],'gray']
        # xkcd_colors = ["pinkish", "brick red",
        #                "goldenrod", "orange",
        #                "buff", "banana yellow",
        #                "pale lime green", "tree green",
        #                "robin egg blue", "windows blue",
        #                "lilac", "purplish",
        #                "tan", "raw umber",
        #                "grey"]
        xkcd_colors = ["pinkish", "brick red", "goldenrod", "orange",
                       "banana yellow", "pale lime green", "tree green",
                       "robin egg blue", "windows blue", "lilac", "purplish",
                       "raw umber", "grey"]
        xkcd_colors = ["cherry", "tree green", "orange", "windows blue",
                   "banana yellow", "bubblegum pink", "goldenrod",
                   "robin egg blue", "pale lime green", "lilac", "purply",
                   "raw umber", "grey"]
        colors = sns.xkcd_palette(xkcd_colors)

        end_colors = sns.xkcd_palette(['vibrant green', 'bright red'])

        # Turn off lines and multiple markers in legend
        mpl.rcParams['legend.handlelength'] = 0
        mpl.rcParams['legend.numpoints'] = 1
        mpl.rcParams['legend.scatterpoints'] = 1
        plt.rcParams['xtick.major.pad']='3'
        plt.rcParams['ytick.major.pad']='3'

        #-- create figure ------
        plt.close("all")
        aspect_ratio = float(xdatarange)/ydatarange
        fig = plt.figure(1, figsize=(fsize1, fsize2))
        ax = fig.add_subplot(111)

        if sel == 'name CA':
            postf = '_ca'
        elif sel == 'backbone':
            postf = '_bb'
        ref_theta_nmp = np.loadtxt(ofname_ref_nmp+postf+'.dat')
        ref_theta_lid = np.loadtxt(ofname_ref_lid+postf+'.dat')
        #-- plot reference structures ----------
        s1 = ax.scatter(ref_theta_nmp[0], ref_theta_lid[0], marker="o", s=100, \
                label=" 1AKE", c=end_colors[0], zorder=15, linewidths=1)
        # s1 = ax.scatter(ref_theta_nmp[0], ref_theta_lid[0], marker="o", s=150, \
        #         c='gray', zorder=14, linewidths=0)
        s2 = ax.scatter(ref_theta_nmp[1], ref_theta_lid[1], marker="D", s=90,  \
                label=" 4AKE", c=end_colors[1], zorder=15, linewidths=1)
        # s2 = ax.scatter(ref_theta_nmp[1], ref_theta_lid[1], marker="D", s=140,  \
        #         c='gray', zorder=14, linewidths=0)

        for i, m in zip(xrange(len(methods)), methods):
            method, methodname, postfix = m[0], m[1], m[2]
            d = 'co'
            theta_nmp_accum = []
            theta_lid_accum = []
            for run in xrange(N):
                mn = ' %s' % methodname
                mtype = 'o' if d=='co' else 'D'
                if run > 0:
                    mn = None
                    if (method == 'intp'): break # only one interp line traj
                if  (method == 'intp' and d == 'oc'): break
                if postfix != None:
                    basename = '%s/%s_%s-%s-%s-%s'                         \
                            % (analysis_dir, 'aa', method, postfix, d, protein)
                else:
                    basename = '%s/%s_%s-%s-%s'                            \
                            % (analysis_dir, 'aa', method, d, protein)
                ### Hack to use all-atom intp trajectories ('_full' vs '_ca')
                if method == 'intp':
                    ofname1 = '%s-%s_%s_full.dat' % (basename, 'nmp', runs[run])
                    ofname2 = '%s-%s_%s_full.dat' % (basename, 'lid', runs[run])
                else:
                    ofname1 = '%s-%s_%s.dat' % (basename, 'nmp', runs[run])
                    ofname2 = '%s-%s_%s.dat' % (basename, 'lid', runs[run])
                sk = 2 if method == 'map' else 1 # reduce # map frames
                theta_trj_nmp = np.loadtxt(ofname1)[st:en:sk]
                theta_trj_lid = np.loadtxt(ofname2)[st:en:sk]
                theta_nmp_accum = np.append(theta_nmp_accum, theta_trj_nmp)
                theta_lid_accum = np.append(theta_lid_accum, theta_trj_lid)
                color = colors[i]
                zorder = 12-i
                (mtype,size) = ('D',15) if i < 6 else ('o',21)
                ax.plot(theta_trj_nmp, theta_trj_lid, lw=1,              \
                    c=color, alpha=1, zorder=zorder-1)
                # zorder = -i #if 'menm' not in method else 9
                ax.scatter(theta_trj_nmp, theta_trj_lid, lw=0.2, s=size,     \
                        c=color, marker=mtype, label=mn, alpha=1, zorder=zorder)


        #----------------------------------------
        # Legends
        #----------------------------------------
        # Thin line around the legend box, and instead fill it with a light grey
        # Also only use one point for the scatterplot legend because the user will
        # get the idea after just one, they don't need three.
        # axl = fig.add_subplot(122)
        # light_grey = np.array([float(248)/float(255)]*3)
        # legend = axl.legend(loc='lower left', bbox_to_anchor=(0.85, 0.025),     \
        #         frameon=True, fancybox=True, scatterpoints=1, borderpad=1,     \
        #         labelspacing=0.5, markerscale=1.0)
        # rect = legend.get_frame()
        # rect.set_facecolor(light_grey)
        # rect.set_linewidth(1.0)
        # Change the legend label colors to almost black, too
        # for t in legend.texts: t.set_color(almost_black)
        # plt.setp(legend.get_texts(), fontsize=(lsize-2))
        #----------------------------------------
        legend = None

        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])

        # Remove top and right axes lines ("spines")
        # spines_to_remove = ['top', 'right']
        # for spine in spines_to_remove:
        #     ax.spines[spine].set_visible(False)

        # Locate major ticks at 5 deg intervals on both axes
        intvl = 10
        x_major_ticklocs = np.arange(40, 75+intvl, intvl)
        y_major_ticklocs = np.arange(100, 155+intvl, intvl)
        locx = plticker.FixedLocator(x_major_ticklocs)
        locy = plticker.FixedLocator(y_major_ticklocs)
        ax.xaxis.set_major_locator(locx)
        ax.yaxis.set_major_locator(locy)

        ax.tick_params(axis='both', which='major', labelsize=lsize, pad=-7.5)

        # major_formatter = plticker.FormatStrFormatter(r'%3i$^\circ$')
        # ax.xaxis.set_major_formatter(major_formatter)
        # ax.yaxis.set_major_formatter(major_formatter)

        xticklabels = np.array(ax.get_xticks())
        yticklabels = np.array(ax.get_yticks())
        # x_major_ticklocs = np.arange(xmin, xmax+intvl, intvl)
        xticklabels = ['' if (t%2).round() != 0 else r'%3i$^\circ$' % t for t in xticklabels] # change fraction to %
        yticklabels = ['' if (t%2).round() != 0 else r'%3i$^\circ$' % t for t in yticklabels] # change fraction to %
        ax.set_xticklabels(xticklabels, minor=False)
        ax.set_yticklabels(yticklabels, minor=False)

        # ax.grid(False)

        # ax.tick_params(axis='both', which='major', labelsize=lsize)

        # ax.grid(True, which='major', ls='-', c=almost_black, alpha=1)

        # rec = mpatches.Rectangle((xmin,ymin), xdatarange, ydatarange,          \
        #         fill=False,lw=1, color='grey', alpha=0.4)
        # rec = ax.add_patch(rec)
        # rec.set_clip_on(False)

        #-- labels ------
        ax.set_xlabel(r'NMP$-$CORE angle', fontsize=lsize)
        ax.set_ylabel(r'LID$-$CORE angle', fontsize=lsize)

        # ax.xaxis.set_ticks_position('none')
        # ax.yaxis.set_ticks_position('none')

        sns.despine(offset=10, bottom=True, left=True)

        plt.tight_layout()

        def multisavefig(name, lgd=None, transparent=True):
            for ext in ['tif', 'pdf', 'png']:
                if lgd is None:
                    plt.savefig(name+'.'+ext, format=ext,                      \
                        transparent=transparent, bbox_inches='tight')
                else:
                    plt.savefig(name+'.'+ext, format=ext,                      \
                        transparent=transparent,                               \
                        bbox_extra_artists=(lgd,), bbox_inches='tight')

        if (plotname is not None):
            if type(plotname) is str:
                pltfp1, pltfp2 = plotname, None
            else: pltfp1, pltfp2 = plotname
            # multisavefig(pltfp1, lgd=legend, transparent=True)
            multisavefig(pltfp1, transparent=True)
            # multisavefig(pltfp1+'-notrans', lgd=legend, transparent=False)
            if pltfp2 is not None:
                import glob, os, shutil, re
                source_dir_split = [x for x in pltfp1.split('/')]
                source_dir = '/'.join(source_dir_split[:-1])
                name = '%s*' % pltfp1
                files = glob.iglob(os.path.join(source_dir, name))
                for f in files:
                    if os.path.isfile(f):
                        ext = f.split('.')[-1]
                        newfile = '%s.%s' % (pltfp2, ext)
                        newfile_t = '%s-notrans.%s' % (pltfp2, ext)
                        shutil.copyfile(f, newfile)
                        shutil.copyfile(f, newfile_t)
        if showplot:
            plt.show()
            plt.close()


if __name__ == '__main__':
    methods = [('dims','DIMS'), ('gp','GP'), ('gomd','GodMD'), ('anmp','ANMP'), ('intp','Interp')]
    protein, directions = 'adk', ('co','oc')
    params = methods, protein, directions
    runs = ['001']#['001','002','003']
    N = len(runs)

    # analysis_dir = '/nfs/homes/sseyler/Projects/PathwayAnalysis/Phase1/data/methods/aa'
    analysis_dir = '/nfs/homes/sseyler/Projects/Methods/PSA/analysis/methods/data/aa'
    pathname = 'pathway_ca.pdb'
    adk_c, adk_o = '1ake_a_core.pdb', '4ake_a_core.pdb'
    plotname = 'aa_methods'
    #-----------------------

    basename, runs, refnames = gen_angle_trajectories(analysis_dir, params, N,
            runs=runs, adk=(adk_c,adk_o), pathname=pathname, force=False)

    timeseries, twod = False, True
    plot_type = timeseries, twod
    plot_aa(analysis_dir, params, runs, refnames, plot_type,
            pltparams=([40,80],[100,160],12,24), plotname=plotname, showplot=True)
