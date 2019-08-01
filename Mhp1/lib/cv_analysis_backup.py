#################################################
# Anton Analysis Script
#################################################

from MDAnalysis import *
from MDAnalysis.analysis.align import *
from MDAnalysis import Universe
import MDAnalysis.analysis.contacts

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#------------------------------------------------
# Options
#------------------------------------------------

showplots = False
saveplots = True

#-----------------------
# 1ake or 4ake analysis
#-----------------------

adk = "1ake"
ADK = adk.upper()
plot_dir = adk + "/plots/"

#-----------------------
# Choose frame subset
#-----------------------

subset = False
multicolor = True

# default values
postfix = "" # default - no postfix
ss = 7

if subset:
    start = 0
    end = 500
    postfix = "_f" + str(start) + "-" + str(end)
    ss = 10


#################################################
# Easier approach for entire trajectory
#################################################
print "CV calculation: delta RMSD"

strct_dir = adk + "/structures/"
trj_dir = adk + "/trj/"
trjraw_dir = trj_dir + "raw/001/"
trjfit_dir = trj_dir + "fit/"

ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb")
ref_clsd.atoms.translate(-ref_clsd.atoms.centerOfMass())
trj = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
rms_fit_trj(trj, ref_clsd, filename=trjfit_dir+adk+"_nw_1-20_all_fit1ake.dcd")

ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")
ref_open.atoms.translate(-ref_open.atoms.centerOfMass())
trj = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
rms_fit_trj(trj, ref_open, filename=trjfit_dir+adk+"_nw_1-20_all_fit4ake.dcd")

u1 = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjfit_dir+adk+"_nw_1-20_all_fit1ake.dcd")
u2 = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjfit_dir+adk+"_nw_1-20_all_fit4ake.dcd")

ref_c_coor =  ref_c.atoms.coordinates()
ref_o_coor =  ref_o.atoms.coordinates()

rmsd_trj = np.zeros((len(u1.trajectory), 2))

idx = 0
for ts in u1.trajectory:
    mobile_coor =  u1.atoms.coordinates()
    rmsd_trj[idx,0] = rmsd(mobile_coor, ref_c_coor)
    idx += 1
idx = 0
for ts in u2.trajectory:
    mobile_coor =  u2.atoms.coordinates()
    rmsd_trj[idx,1] = rmsd(mobile_coor, ref_o_coor)
    idx += 1

np.savetxt(adk+"/"+adk+"_nw_1-20_all_rmsd.txt", rmsd_trj, fmt="%7.5f")

#----------------------------
# Get RMSD b/w 1ake and 4ake structures (after centering and rotation)
#----------------------------
ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb")
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")
ref_clsd_coor = ref_clsd.atoms.coordinates() - ref_clsd.atoms.centerOfMass()
ref_open_coor = ref_open.atoms.coordinates() - ref_open.atoms.centerOfMass()
R, rmsd_val = rotation_matrix(ref_clsd_coor, ref_open_coor)
np.savetxt("rmsd_1ake-4ake.txt", np.array([rmsd_val]), fmt="%7.5f")

#------------------------------------------------
# Plot RMSD w.r.t. boundary conformations for entire length of trajectory
#------------------------------------------------
print "Plotting..."

trj_rmsd = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_rmsd.txt")
if subset:
    trj_rmsd = trj_rmsd[start:end,:]
time = 240*np.arange(0, len(trj_rmsd))/1000.

ref_rmsd = np.loadtxt("rmsd_1ake-4ake.txt")

#-----------------------
# Timeseries plots
#-----------------------
print "\ttimeseries"

fig = plt.figure(figsize=(20,10))
plt.plot(time, ref_rmsd*np.ones_like(time), label="RMSD 1ake/4ake", color='magenta')
plt.plot(time, trj_rmsd[:,0], label="w.r.t. 1ake", color='g')
plt.plot(time, trj_rmsd[:,1], label="w.r.t. 4ake", color='r')
plt.title(r'Evolution of $\Delta$RMSD from 1ake/4ake -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'$\Delta$RMSD ($\AA$)')
plt.legend()
if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_rmsd_TS"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# RMSD w.r.t 1ake vs RMSD w.r.t 4ake plot
#-----------------------
print "\t2D CV"

#-- create figure ------
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
ax.scatter(ref_rmsd, 0, marker="*", s=100, label="4ake", c='cyan')
ax.scatter(0, ref_rmsd, marker="x", s=80, label="1ake", c='magenta')
if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(trj_rmsd[:,0], trj_rmsd[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(trj_rmsd[:,0], trj_rmsd[:,1], c='blue', s=ss, \
	  edgecolor='none')

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_ylabel(r'$\Delta$RMSD from 4ake ($\AA$)', fontsize=20)
ax.set_xlabel(r'$\Delta$RMSD from 1ake ($\AA$)', fontsize=20)
ax.set_title(r'$\Delta$RMSD from End Structures (1ake vs 4ake) -- '+adk+' trj'+postfix)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+'/plots/'+adk+"_nw_1-20_all_rmsd_2D"+postfix+".pdf")
plt.cla()
plt.clf()


#################################################
# Look at angle-angle coordinates
#################################################
print "CV calculation: angle-angle"

# NMP: centers of geometry of the backbone and C-beta atoms in residues 115-125
# (CORE-LID), 90-100 (CORE), and 35-55 (NMP)

# LID: 179-185 (CORE), 115-125 (CORE-hinge-LID), and 125-153 (LID)

strct_dir = adk + "/structures/"
trj_dir = adk + "/trj/"
trjraw_dir = trj_dir + "raw/001/"
trjfit_dir = trj_dir + "fit/"

ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb")
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")

u = Universe(strct_dir+"adk4AKE_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
bb = u.selectAtoms("backbone")

# Theta NMP: 115-125 (CORE-hinge-LID), 90-100 (CORE), and 35-55 (NMP)
core_lid = bb.selectAtoms("resid 115-125")
core1 = bb.selectAtoms("resid 90-100")
nmp = bb.selectAtoms("resid 35-55")

# Theta LID: 179-185 (CORE), 115-125 (CORE-hinge-LID), and 125-153 (LID)
core2 = bb.selectAtoms("resid 179-185")
core_lid = bb.selectAtoms("resid 115-125")
lid = bb.selectAtoms("resid 125-153")

def vsqnorm(v, axis=0):
    return np.sum(v*v, axis=axis)

def computeangle(v1, v2):
    return np.arccos( np.dot(v1, v2)/((vsqnorm(v1)*vsqnorm(v2))**0.5) )*360/(2*np.pi)

theta_nmp_trj = np.zeros(len(u.trajectory))
idx = 0
for ts in u.trajectory:
    core1centroid = core1.centroid()
    v1 = core_lid.centroid() - core1centroid
    v2 = nmp.centroid() - core1centroid
    theta_nmp_trj[idx] = computeangle(v1, v2)
    idx += 1
theta_nmp_trj = theta_nmp_trj

theta_lid_trj = np.zeros(len(u.trajectory))
idx = 0
for ts in u.trajectory:
    core_lid_centroid = core_lid.centroid()
    v1 = lid.centroid() - core_lid_centroid
    v2 = core2.centroid() - core_lid_centroid
    theta_lid_trj[idx] = computeangle(v1, v2)
    idx += 1
theta_lid_trj = theta_lid_trj

angles = np.transpose(np.vstack((theta_nmp_trj, theta_lid_trj)))
np.savetxt(adk+"/"+adk+"_nw_1-20_all_angles.txt", angles, fmt="%7.5f")

#----------------------------
# Get LID-CORE and NMP-CORE angles for 1ake and 4ake
#----------------------------
ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb")
c_bb = ref_clsd.selectAtoms("backbone")
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")
o_bb = ref_open.selectAtoms("backbone")

# Compute NMP-CORE angle for 1ake
c_core1 = c_bb.selectAtoms("resid 90-100")
c_core_lid = c_bb.selectAtoms("resid 115-125")
c_nmp = c_bb.selectAtoms("resid 35-55")
a1 = c_core_lid.centroid() - c_core1.centroid()
a2 = c_nmp.centroid() - c_core1.centroid()
c_theta_nmp = computeangle(a1, a2)

# Compute NMP-CORE angle for 4ake
o_core1 = o_bb.selectAtoms("resid 90-100")
o_core_lid = o_bb.selectAtoms("resid 115-125")
o_nmp = o_bb.selectAtoms("resid 35-55")
a1 = o_core_lid.centroid() - o_core1.centroid()
a2 = o_nmp.centroid() - o_core1.centroid()
o_theta_nmp = computeangle(a1, a2)

# Compute LID-CORE angle for 1ake
c_core2 = c_bb.selectAtoms("resid 179-185")
c_core_lid = c_bb.selectAtoms("resid 115-125")
c_lid = c_bb.selectAtoms("resid 125-153")
a1 = c_lid.centroid() - c_core_lid.centroid()
a2 = c_core2.centroid() - c_core_lid.centroid()
c_theta_lid = computeangle(a1, a2)

# Compute LID-CORE angle for 4ake
o_core2 = o_bb.selectAtoms("resid 179-185")
o_core_lid = o_bb.selectAtoms("resid 115-125")
o_lid = o_bb.selectAtoms("resid 125-153")
a1 = o_lid.centroid() - o_core_lid.centroid()
a2 = o_core2.centroid() - o_core_lid.centroid()
o_theta_lid = computeangle(a1, a2)

ref_angles = np.array([c_theta_nmp, c_theta_lid, o_theta_nmp, o_theta_lid])
np.savetxt("ref_angles.txt", ref_angles, fmt="%7.5f")

#------------------------------------------------
# Plot angle-angle coordinates for entire length of trajectory
#------------------------------------------------
print "Plotting..."

angles = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_angles.txt")
if subset:
    angles = angles[start:end,:]
time = 240*np.arange(0, len(angles))/1000.

ref_angles = np.loadtxt("ref_angles.txt")

#-----------------------
# Timeseries plots
#-----------------------
print "\ttimeseries"

fig = plt.figure(figsize=(20, 10))

plt.plot(time, angles[:,0], label="NMP-CORE angle", color='g')
plt.plot(time, angles[:,1], label="LID-CORE angle", color='r')

plt.plot(time, ref_angles[0]*np.ones_like(time), label="1ake NMP-CORE angle", color='cyan', ls='-')
plt.plot(time, ref_angles[1]*np.ones_like(time), label="1ake LID-CORE angle", color='cyan', ls='-.', linewidth=2)
plt.plot(time, ref_angles[2]*np.ones_like(time), label="4ake NMP-CORE angle", color='magenta', ls='-')
plt.plot(time, ref_angles[3]*np.ones_like(time), label="4ake LID-CORE angle", color='magenta', ls='-.', linewidth=2)

plt.title(r'Evolution of NMP-/LID-CORE Angles -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'$\theta^\circ$')
plt.legend()
if showplots:
    plt.show()

if saveplots:
    fig.savefig(plot_dir+adk+"_nw_1-20_all_angles_TS"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# Angle-Angle plot
#-----------------------
print "\t2D CV"

xmin = 35
xmax = 95
ymin = 90
ymax = 180

#-- create figure ------
fig = plt.figure(figsize=(8, 12))
ax = fig.add_subplot(111)
ax.set_aspect(9/6)
# ax.set_position([0.1,0.1,0.8,0.8])

#-- plot data ----------
if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(angles[:,0], angles[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(angles[:,0], angles[:,1], c='blue', s=ss, \
	  edgecolor='none')

ax.scatter(ref_angles[2], ref_angles[3], marker="*", s=100, label="4ake", c='cyan')
ax.scatter(ref_angles[0], ref_angles[1], marker="x", s=80, label="1ake", c='magenta')

#-- plot settings ------
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_xlabel(r'$\theta_\mathrm{NMP}$', fontsize=20)
ax.set_ylabel(r'$\theta_\mathrm{LID}$', fontsize=20)
ax.set_title(r'NMP-CORE Angle vs LID-CORE Angle -- '+adk+' trj'+postfix)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(plot_dir+adk+"_nw_1-20_all_angles_2D"+postfix+".pdf")
plt.cla()
plt.clf()


#################################################
# Domain CoM distances (LID-CORE and NMP-CORE)
#################################################
print "CV calculation: domain mass-center distances"

strct_dir = adk + "/structures/"
trj_dir = adk + "/trj/"
trjraw_dir = trj_dir + "raw/001/"
trjfit_dir = trj_dir + "fit/"

u = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
ca = u.selectAtoms("name CA")

# Theta NMP: 115-125 (CORE-hinge-LID), 90-100 (CORE), and 35-55 (NMP)
lid = "resid 122-159"
core = "resid 1-29 or resid 60-121 or resid 160-214"
nmp = "resid 30-59"

mobile_lid = ca.selectAtoms(lid)
mobile_core = ca.selectAtoms(core)
mobile_nmp = ca.selectAtoms(nmp)

def vsqnorm(v, axis=0):
    return np.sum(v*v, axis=axis)

dmcd = np.zeros((len(u.trajectory), 2))
idx = 0
for ts in u.trajectory:
    lid_com = mobile_lid.centerOfMass()
    core_com = mobile_core.centerOfMass()
    nmp_com = mobile_nmp.centerOfMass()
    dmcd[idx,:] = np.array([vsqnorm(nmp_com-core_com), \
        vsqnorm(lid_com-core_com)])
    idx += 1
dmcd = dmcd**0.5

np.savetxt(adk+"/"+adk+"_nw_1-20_all_dmcd.txt", dmcd, fmt="%7.5f")

#----------------------------
# Get LID-CORE and NMP-CORE MCDs for 1ake and 4ake
#----------------------------
ref_c = Universe(strct_dir+"adk1AKE_nw.pdb")
ref_o = Universe(strct_dir+"adk4AKE_nw.pdb")

ref_nmp_c = ref_c.selectAtoms(nmp + " and name CA").centerOfMass()
ref_lid_c = ref_c.selectAtoms(lid + " and name CA").centerOfMass()
ref_core_c = ref_c.selectAtoms(core + " and name CA").centerOfMass()
ref_nmp_o = ref_o.selectAtoms(nmp + " and name CA").centerOfMass()
ref_lid_o = ref_o.selectAtoms(lid + " and name CA").centerOfMass()
ref_core_o = ref_o.selectAtoms(core + " and name CA").centerOfMass()

ref_dmcd = np.array([vsqnorm(ref_nmp_c-ref_core_c), \
						vsqnorm(ref_lid_c-ref_core_c), \
						vsqnorm(ref_nmp_o-ref_core_o), \
						vsqnorm(ref_lid_o-ref_core_o)
					])

np.savetxt("dmcd.txt", ref_dmcd**0.5, fmt="%7.5f")

#------------------------------------------------
# Plot domain MCDs for entire length of trajectory
#------------------------------------------------
print "Plotting..."

dmcd = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_dmcd.txt")
if subset:
    dmcd = dmcd[start:end,:]
time = 240*np.arange(0, len(angles))/1000.

ref_dmcd = np.loadtxt("dmcd.txt")

#-----------------------
# Timeseries plots
#-----------------------
print "\ttimeseries"

fig = plt.figure(figsize=(20,10))

plt.plot(time, dmcd[:,0], label="NMP-CORE MCD", color='g')
plt.plot(time, dmcd[:,1], label="LID-CORE MCD", color='r')

plt.plot(time, ref_dmcd[0]*np.ones_like(time), label="1ake NMP-CORE MCD", color='cyan', ls='-')
plt.plot(time, ref_dmcd[1]*np.ones_like(time), label="1ake LID-CORE MCD", color='cyan', ls='-.', linewidth=2)
plt.plot(time, ref_dmcd[2]*np.ones_like(time), label="4ake NMP-CORE MCD", color='magenta', ls='-')
plt.plot(time, ref_dmcd[3]*np.ones_like(time), label="4ake LID-CORE MCD", color='magenta', ls='-.', linewidth=2)

plt.title('Evolution of NMP-/LID-CORE Mass-Center Distance -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'LID/NMP-CORE MCD ($\AA$)')
plt.legend()
if showplots:
    plt.show()

if saveplots:
    fig.savefig(plot_dir+adk+"_nw_1-20_all_dmcd_TS"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# LID-CORE vs NMP-CORE plot
#-----------------------
print "\t2D CV"

xmin = 18
xmax = 26
ymin = 18
ymax = 38

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.set_aspect(0.4)
# ax.set_position([0.1,0.1,0.8,0.8])

#-- plot data ----------
if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(dmcd[:,0], dmcd[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(dmcd[:,0], dmcd[:,1], c='blue', s=ss, \
	  edgecolor='none')

ax.scatter(ref_dmcd[2], ref_dmcd[3], marker="*", s=100, label="4ake", c='cyan')
ax.scatter(ref_dmcd[0], ref_dmcd[1], marker="x", s=80, label="1ake", c='magenta')

#-- plot settings ------

ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_xlabel(r'NMP-CORE MCD ($\AA$)', fontsize=20)
ax.set_ylabel(r'LID-CORE MCD ($\AA$)', fontsize=20)
ax.set_title('NMP-CORE vs LID-CORE Mass-Center Distances -- '+adk+' trj'+postfix)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(plot_dir+adk+"_nw_1-20_all_dmcd_2D"+postfix+".pdf")
plt.cla()
plt.clf()


#################################################
# Native contacts
#################################################
print "CV calculation: native contacts"

strct_dir = adk + "/structures/"
trj_dir = adk + "/trj/"
trjraw_dir = trj_dir + "raw/001/"
trjfit_dir = trj_dir + "fit/"

ref_clsd = Universe([strct_dir+"adk1AKE_nw.pdb")
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")
atomselection = 'name CA' # change this to 'name CA and segid A' for topology files that specify a chain
calp1 = ref_clsd.selectAtoms(atomselection)
calp2 = ref_open.selectAtoms(atomselection)

u = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")

cutoff = 8.0 # in Angstroms
prefix = '/nfs/homes/sseyler/Simulations/adk/anton/'
ofn = prefix + 'q1q2/' + adk + '/'
ofn1 = ofn + '1ake_q1'
ofn2 = ofn + '4ake_q2'
ext1 = '.dat'
outputname = adk+"_nw_1-20_all_q1q2.dat"

CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(u, selection=atomselection,refgroup=calp1, \
    radius=cutoff,outfile=ofn1+ext1)
CA2 = MDAnalysis.analysis.contacts.ContactAnalysis1(u, selection=atomselection,refgroup=calp2, \
    radius=cutoff,outfile=ofn2+ext1)
CA1.run(force=True)
CA2.run(force=True)
q1 = CA1.timeseries[1]
q2 = CA2.timeseries[1]

np.savetxt(outputname, np.vstack((q1,q2)).T, fmt="%7.5e")

#----------------------------
# Get LID-CORE and NMP-CORE MCDs for 1ake and 4ake
#----------------------------
# NEEDS TO BE FIXED
ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb", [strct_dir+"adk1AKE_nw.pdb",strct_dir+"adk1AKE_nw.pdb"])
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb", [strct_dir+"adk4AKE_nw.pdb",strct_dir+"adk4AKE_nw.pdb"])

atomselection = 'name CA' # change this to 'name CA and segid A' for topology files that specify a chain
calp1 = ref_clsd.selectAtoms(atomselection)
calp2 = ref_open.selectAtoms(atomselection)
cutoff = 8.0 # in Angstroms
prefix = '/nfs/homes/sseyler/Simulations/adk/anton/'
ofn = prefix + 'q1q2/' + adk + '/'
ofn1 = ofn + '1ake_q1_0'
ofn2 = ofn + '4ake_q2_0'
ext1 = '.dat'
outputname = adk+"_nw_1-20_all_q1q2.dat"

CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(ref_open, selection=atomselection,refgroup=calp1, \
    radius=cutoff,outfile=ofn1+ext1)
CA2 = MDAnalysis.analysis.contacts.ContactAnalysis1(ref_clsd, selection=atomselection,refgroup=calp2, \
    radius=cutoff,outfile=ofn2+ext1)
CA1.run(force=True)
CA2.run(force=True)
q1 = CA1.timeseries[1]
q2 = CA2.timeseries[1]

np.savetxt("ref_nc.txt", np.array([q1[0], q2[0]]), fmt="%7.5e")

#-----------------------
# Plotting
#-----------------------
print "Plotting..."

from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('medium')

#-----------------------
# Timeseries plots
#-----------------------
print "\ttimeseries"

xmin = ymin = 0.8
xmax = ymax = 1.0

q = np.loadtxt(outputname)
if subset:
    q = q[start:end,:]
time = 240*np.arange(0, len(q))/1000.

ref_nc = np.loadtxt("ref_nc.txt")

#-- plot data ----------
fig = plt.figure(figsize=(20,10))
ax = fig.add_subplot(111)

ax.plot(time, ref_nc[0]*np.ones_like(time), label="Fraction of 1ake contacts in 4ake", color='cyan', ls='-')
ax.plot(time, ref_nc[1]*np.ones_like(time), label="Fraction of 4ake contacts in 1ake", color='cyan', ls='-.')

ax.plot(time, q[:,0], linewidth=0.5, linestyle='-', color='g', label='q_1ake')
ax.plot(time, q[:,1], linewidth=0.5, linestyle='-', color='r', label='q_4ake')

#-- labels/legend ------
ax.set_ylabel('Fraction native contacts', fontsize=20)
ax.set_xlabel('Time (ns)', fontsize=20)
ax.set_title('Evolution of Native Contacts Fraction w.r.t 1ake/4ake -- '+adk+' trj'+postfix)
leg = ax.legend(loc='center left', bbox_to_anchor=(0.95, 0.5), prop=fontP)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(plot_dir+adk+"_nw_1-20_all_q1q2_TS"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# q1-q2 plot
#-----------------------
print "\t2D CV"

q = np.loadtxt(outputname)
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.set_position([0.1,0.1,0.8,0.8])
ax.set_aspect(1)

#-- plot data ----------
if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(q[:,0], q[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(q[:,0], q[:,1], c='blue', s=ss, \
	  edgecolor='none')

#-- plot settings ------
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_ylabel('Fraction of 4ake Native Contacts', fontsize=20)
ax.set_xlabel('Fraction of 1ake Native Contacts', fontsize=20)
ax.set_title('2D Native Contacts Space -- '+adk+' trj'+postfix)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(plot_dir+adk+"_nw_1-20_all_q1q2_2D"+postfix+".pdf")
plt.cla()
plt.clf()


#################################################
# Intradomain flexibility
#################################################
print "CV calculation: delta RMSD LID/NMP (intradomain flexibility)"

strct_dir = adk + "/structures/"
trj_dir = adk + "/trj/"
trjraw_dir = trj_dir + "raw/001/"
trjfit_dir = trj_dir + "fit/"

outprefix = adk + "/" + adk + "_nw_1-20_all_"
ext = ".dcd"

lid = "resid 122-159"
core = "resid 1-29 or resid 60-121 or resid 160-214"
nmp = "resid 30-59"

# Loop through trajectories

# 1ake reference
ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb")

clsd_lid = ref_clsd.selectAtoms(lid)
clsd_core = ref_clsd.selectAtoms(core)
clsd_nmp = ref_clsd.selectAtoms(nmp)

# outname = outprefix + "fit1ake-lid" + ext

ref_lid_com_c = clsd_lid.CA.coordinates() - clsd_lid.CA.centerOfMass()
u1 = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
u1_lid = u1.selectAtoms(lid)

rmsd_trj_lid = np.zeros((len(u1.trajectory), 2))
idx = 0
for ts in u1.trajectory:
    mobile_lid = u1_lid.CA.coordinates() - u1_lid.CA.centerOfMass()
    R, rmsd_val = rotation_matrix(mobile_lid, ref_lid_com_c)
    # u1_lid.atoms.translate(-lid.centerOfMass())
    # u1_lid.atoms.rotate(R)
    rmsd_trj_lid[idx,0] = rmsd_val
    idx += 1

ref_nmp_com_c = clsd_nmp.CA.coordinates() - clsd_nmp.CA.centerOfMass()
u1 = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
u1_nmp = u1.selectAtoms(nmp)

rmsd_trj_nmp = np.zeros((len(u1.trajectory), 2))
idx = 0
for ts in u1.trajectory:
    mobile_nmp = u1_nmp.CA.coordinates() - u1_nmp.CA.centerOfMass()
    R, rmsd_val = rotation_matrix(mobile_nmp, ref_nmp_com_c)
    # u1_nmp.atoms.translate(-nmp.centerOfMass())
    # u1_nmp.atoms.rotate(R)
    rmsd_trj_nmp[idx,0] = rmsd_val
    idx += 1

# 4ake reference
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")

open_lid = ref_open.selectAtoms(lid)
open_core = ref_open.selectAtoms(core)
open_nmp = ref_open.selectAtoms(nmp)

ref_lid_com_o = open_lid.CA.coordinates() - open_lid.CA.centerOfMass()
u2 = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
u2_lid = u2.selectAtoms(lid)

idx = 0
for ts in u2.trajectory:
    mobile_lid = u2_lid.CA.coordinates() - u2_lid.CA.centerOfMass()
    R, rmsd_val = rotation_matrix(mobile_lid, ref_lid_com_o)
    # u2_lid.atoms.translate(-lid.centerOfMass())
    # u2_lid.atoms.rotate(R)
    rmsd_trj_lid[idx,1] = rmsd_val
    idx += 1

ref_nmp_com_o = open_nmp.CA.coordinates() - open_nmp.CA.centerOfMass()
u2 = Universe(strct_dir+"adk"+ADK+"_nw.pdb", trjraw_dir+adk+"_nw_1-20_all.dcd")
u2_nmp = u2.selectAtoms(nmp)

idx = 0
for ts in u2.trajectory:
    mobile_nmp = u2_nmp.CA.coordinates() - u2_nmp.CA.centerOfMass()
    R, rmsd_val = rotation_matrix(mobile_nmp, ref_nmp_com_o)
    # u2_nmp.atoms.translate(-nmp.centerOfMass())
    # u2_nmp.atoms.rotate(R)
    rmsd_trj_nmp[idx,1] = rmsd_val
    idx += 1

np.savetxt(outprefix+"idf-nmp.txt", rmsd_trj_nmp, fmt="%7.5f")
np.savetxt(outprefix+"idf-lid.txt", rmsd_trj_lid, fmt="%7.5f")

#----------------------------
# RMSDs of domains b/w 1ake/4ake
#----------------------------
lid = "resid 122-159"
core = "resid 1-29 or resid 60-121 or resid 160-214"
nmp = "resid 30-59"

ref_clsd = Universe(strct_dir+"adk1AKE_nw.pdb")
clsd_lid = ref_clsd.selectAtoms(lid)
clsd_core = ref_clsd.selectAtoms(core)
clsd_nmp = ref_clsd.selectAtoms(nmp)
ref_open = Universe(strct_dir+"adk4AKE_nw.pdb")
open_lid = ref_open.selectAtoms(lid)
open_core = ref_open.selectAtoms(core)
open_nmp = ref_open.selectAtoms(nmp)

ref_lid_com_c = clsd_lid.CA.coordinates() - clsd_lid.CA.centerOfMass()
ref_nmp_com_c = clsd_nmp.CA.coordinates() - clsd_nmp.CA.centerOfMass()
ref_lid_com_o = open_lid.CA.coordinates() - open_lid.CA.centerOfMass()
ref_nmp_com_o = open_nmp.CA.coordinates() - open_nmp.CA.centerOfMass()

R, rmsd_nmp = rotation_matrix(ref_nmp_com_c, ref_nmp_com_o)
R, rmsd_lid = rotation_matrix(ref_lid_com_c, ref_lid_com_o)

np.savetxt("ref_idf.txt", np.array([rmsd_nmp, rmsd_lid]), fmt="%7.5f")

#------------------------------------------------
# Plot RMSD w.r.t. boundary conformations for entire length of trajectory
#------------------------------------------------
print "Plotting..."

#-----------------------
# Timeseries plots
#-----------------------
print "\ttimeseries"

rmsd_trj_nmp = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_idf-nmp.txt")
rmsd_trj_lid = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_idf-lid.txt")
if subset:
    rmsd_trj_nmp = rmsd_trj_nmp[start:end,:]
    rmsd_trj_lid = rmsd_trj_lid[start:end,:]
time = 240*np.arange(0, len(rmsd_trj_nmp))/1000.

ref_idf = np.loadtxt("ref_idf.txt")

fig = plt.figure(figsize=(16, 8))
plt.plot(time, ref_idf[0]*np.ones_like(time), label="1ake/4ake NMP separation", color='cyan', ls='-')

plt.plot(time, rmsd_trj_nmp[:,0], label="from 1ake", color='g')
plt.plot(time, rmsd_trj_nmp[:,1], label="from 4ake", color='r')
plt.title(r'Evolution of $\Delta$RMSD from 1ake/4ake NMP -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'$\Delta$RMSD ($\AA$)')
plt.legend()
if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_1ake-4ake-nmp_TS"+postfix+".pdf")
plt.cla()
plt.clf()

fig = plt.figure(figsize=(16, 8))
plt.plot(time, ref_idf[1]*np.ones_like(time), label="1ake/4ake LID separation", color='cyan', ls='-')

plt.plot(time, rmsd_trj_lid[:,0], label="from 1ake", color='g')
plt.plot(time, rmsd_trj_lid[:,1], label="from 4ake", color='r')
plt.title(r'Evolution of $\Delta$RMSD from 1ake/4ake LID -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'$\Delta$RMSD ($\AA$)')
plt.legend()
if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_1ake-4ake-lid_TS"+postfix+".pdf")
plt.cla()
plt.clf()

fig = plt.figure(figsize=(16, 8))
plt.plot(time, ref_idf[0]*np.ones_like(time), label="RMSD 1ake/4ake NMP domains", color='cyan', ls='-')
plt.plot(time, ref_idf[1]*np.ones_like(time), label="RMSD 1ake/4ake LID domains", color='magenta', ls='-')

plt.plot(time, rmsd_trj_lid[:,0], label="LID", color='g')
plt.plot(time, rmsd_trj_nmp[:,0], label="NMP", color='r')
plt.title(r'Evolution of $\Delta$RMSD from 1ake LID/NMP -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'$\Delta$RMSD ($\AA$)')
plt.legend()
if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_1ake-lid-nmp_TS"+postfix+".pdf")
plt.cla()
plt.clf()

fig = plt.figure(figsize=(16, 8))
plt.plot(time, ref_idf[0]*np.ones_like(time), label="4ake NMP separation from 1ake", color='cyan', ls='-')
plt.plot(time, ref_idf[1]*np.ones_like(time), label="4ake LID separation from 1ake", color='magenta', ls='-')

plt.plot(time, rmsd_trj_lid[:,1], label="LID", color='g')
plt.plot(time, rmsd_trj_nmp[:,1], label="NMP", color='r')
plt.title(r'Evolution of $\Delta$RMSD from 4ake LID/NMP -- '+adk+' trj'+postfix)
plt.xlabel('Time (ns)')
plt.ylabel(r'$\Delta$RMSD ($\AA$)')
plt.legend()
if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_4ake-lid-nmp_TS"+postfix+".pdf")
plt.cla()
plt.clf()


#-----------------------
# RMSD w.r.t 1ake vs RMSD w.r.t 4ake plot
#-----------------------
print "\t2D CV - 1ake v 4ake"

# RMSD 1ake/4ake LID
#-----------------------
trj_rmsd_lid = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_idf-lid.txt")
#-- create figure ------
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
ax.scatter(ref_idf[1], 0, marker="x", s=80, label="RMSD 4ake LID from 1ake", c='magenta')
ax.scatter(0, ref_idf[1], marker="*", s=100, label="RMSD 1ake LID from 4ake", c='cyan')

if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(trj_rmsd_lid[:,0], trj_rmsd_lid[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(trj_rmsd_lid[:,0], trj_rmsd_lid[:,1], c='blue', s=ss, \
	  edgecolor='none')

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_ylabel(r'$\Delta$RMSD from 4ake LID domain ($\AA$)', fontsize=20)
ax.set_xlabel(r'$\Delta$RMSD from 1ake LID domain ($\AA$)', fontsize=20)
ax.set_title(r'$\Delta$RMSD LID - 1ake vs 4ake -- '+adk+' trj'+postfix)
ax.legend(loc='upper left', prop=fontP)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_1ake-4ake-lid_2D"+postfix+".pdf")
plt.cla()
plt.clf()

# RMSD 1ake/4ake NMP
#-----------------------
trj_rmsd_nmp = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_idf-nmp.txt")
#-- create figure ------
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
ax.scatter(ref_idf[0], 0, marker="x", s=80, label="RMSD 4ake NMP from 1ake", c='magenta')
ax.scatter(0, ref_idf[0], marker="*", s=100, label="RMSD 1ake NMP from 4ake", c='cyan')

if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(trj_rmsd_nmp[:,0], trj_rmsd_nmp[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(trj_rmsd_nmp[:,0], trj_rmsd_nmp[:,1], c='blue', s=ss, \
	  edgecolor='none')

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_ylabel(r'$\Delta$RMSD from 4ake NMP domain ($\AA$)', fontsize=20)
ax.set_xlabel(r'$\Delta$RMSD from 1ake NMP domain ($\AA$)', fontsize=20)
ax.set_title(r'$\Delta$RMSD NMP - 1ake vs 4ake -- '+adk+' trj'+postfix)
ax.legend(loc='upper left', prop=fontP)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_1ake-4ake-nmp_2D"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# RMSD LID vs RMSD NMP
#-----------------------
print "\t2D CV - LID v NMP"

trj_rmsd_lid = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_idf-lid.txt")
trj_rmsd_nmp = np.loadtxt(adk+"/"+adk+"_nw_1-20_all_idf-nmp.txt")

# RMSD 1ake LID/NMP
#-----------------------
#-- create figure ------
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
ax.scatter(ref_idf[0], ref_idf[1], marker="x", s=80, label="RMSD 4ake lids from 1ake", c='magenta')
ax.scatter(0, 0, marker="*", s=100, label="RMSD 1ake lids from 1ake", c='cyan')

if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(trj_rmsd_lid[:,0], trj_rmsd_nmp[:,0], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(trj_rmsd_lid[:,0], trj_rmsd_nmp[:,0], c='blue', s=ss, \
	  edgecolor='none')

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_ylabel(r'$\Delta$RMSD from 1ake NMP domain ($\AA$)', fontsize=20)
ax.set_xlabel(r'$\Delta$RMSD from 1ake LID domain ($\AA$)', fontsize=20)
ax.set_title(r'$\Delta$RMSD LID vs $\Delta$RMSD NMP 1ake -- '+adk+' trj'+postfix)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_1ake-lid-nmp_2D"+postfix+".pdf")
plt.cla()
plt.clf()

# RMSD 4ake LID/NMP
#-----------------------
#-- create figure ------
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
ax.scatter(ref_idf[0], ref_idf[1], marker="*", s=100, label="RMSD 1ake lids from 4ake", c='cyan')
ax.scatter(0, 0, marker="x", s=80, label="RMSD 4ake lids from 4ake", c='magenta')

if multicolor:
	colors = cm.RdYlGn_r( 1.0*time/np.max(time) )
	ax.scatter(trj_rmsd_lid[:,1], trj_rmsd_nmp[:,1], c=colors, s=ss, \
	    edgecolor='black', linewidth=0.2)
else:
	ax.scatter(trj_rmsd_lid[:,1], trj_rmsd_nmp[:,1], c='blue', s=ss, \
	  edgecolor='none')

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=16)
ax.grid(True)

#-- labels/legend ------
ax.set_ylabel(r'$\Delta$RMSD from 4ake NMP domain ($\AA$)', fontsize=20)
ax.set_xlabel(r'$\Delta$RMSD from 4ake LID domain ($\AA$)', fontsize=20)
ax.set_title(r'$\Delta$RMSD LID vs $\Delta$RMSD NMP 4ake -- '+adk+' trj'+postfix)

if showplots:
    plt.show()
if saveplots:
    fig.savefig(adk+"/plots/"+adk+"_nw_1-20_all_idf_4ake-lid-nmp_2D"+postfix+".pdf")
plt.cla()
plt.clf()


#################################################
# Combination plots (2D spaces)
#################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('medium')

#----------------------------
# Options
#----------------------------

showplots = False
saveplots = True

#-----------------------
# Choose frame subset
#-----------------------

subset = True
skip = 2

# default values
postfix = "" # default - no postfix
ss = 10
a1 = 0.3
a2 = 0.1

if subset:
    start = 0
    end = 500
    postfix = "_f" + str(start) + "-" + str(end)
    ss = 10

if skip:
    postfix = "_s" + str(skip)
    a1 = a2 = 0.8
else:
    skip = 1


#-----------------------
# RMSD w.r.t 1ake/4ake for 1ake and 4ake trajectories
#-----------------------
print "Plotting delta RMSDs for 1ake and 4ake trajectories"

trj_rmsd_1ake = np.loadtxt("1ake/1ake_nw_1-20_all_rmsd.txt")
trj_rmsd_4ake = np.loadtxt("4ake/4ake_nw_1-20_all_rmsd.txt")
if subset or skip:
    trj_rmsd_1ake = trj_rmsd_1ake[start:end:skip,:]
    trj_rmsd_4ake = trj_rmsd_4ake[start:end:skip,:]
time = 240*np.arange(0, len(trj_rmsd_1ake))/1000.

# Read in reference stuff
ref_rmsd = np.loadtxt("rmsd_1ake-4ake.txt")

#-- create figure ------
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
ax.scatter(ref_rmsd, 0, marker="*", s=100, label="4ake", c='cyan')
ax.scatter(0, ref_rmsd, marker="x", s=80, label="1ake", c='magenta')

p1 = ax.scatter(trj_rmsd_1ake[:,0], trj_rmsd_1ake[:,1], c='red', s=ss, \
    edgecolor='none', linewidth=0.2, label='1ake', alpha=a1)
p2 = ax.scatter(trj_rmsd_4ake[:,0], trj_rmsd_4ake[:,1], c='green', s=ss, \
    edgecolor='none', linewidth=0.2, label='4ake', alpha=a2/3)

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=20)
ax.grid(True)
leg = ax.legend(loc='center right', bbox_to_anchor=(0.95, 0.5), prop=fontP)

#-- labels/legend ------
ax.set_ylabel(r'$\Delta$RMSD from 4ake ($\AA$)', fontsize=24)
ax.set_xlabel(r'$\Delta$RMSD from 1ake ($\AA$)', fontsize=24)
ax.set_title(r'$\Delta$RMSD from End Structures (1ake vs 4ake)')

if showplots:
    plt.show()
if saveplots:
    fig.savefig("1ake4ake_nw_1-20_all_rmsd_2D"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# Angles for 1ake and 4ake trajectories
#-----------------------
print "Plotting angles for 1ake and 4ake trajectories"

angles_1ake_full = np.loadtxt("1ake/1ake_nw_1-20_all_angles.txt")
angles_4ake_full = np.loadtxt("4ake/4ake_nw_1-20_all_angles.txt")
angles_1ake = angles_1ake_full[start:end:skip,:]
angles_4ake = angles_4ake_full[start:end:skip,:]
time = 240*np.arange(0, len(angles_4ake))/1000.

ref_angles = np.loadtxt("ref_angles.txt")

#-- create figure ------
fig = plt.figure(figsize=(8, 12))
ax = fig.add_subplot(111)

#-- plot data ----------
p1 = ax.scatter(angles_1ake[:,0], angles_1ake[:,1], c='red', s=ss, \
    edgecolor='none', linewidth=0.2, label='1ake trajectory', alpha=a1)
p2 = ax.scatter(angles_4ake[:,0], angles_4ake[:,1], c='green', s=ss, \
    edgecolor='none', linewidth=0.2, label='4ake trajectory', alpha=a2/2)

ax.scatter(ref_angles[2], ref_angles[3], marker="*", s=100, label="4ake", c='cyan')
ax.scatter(ref_angles[0], ref_angles[1], marker="x", s=80, label="1ake", c='magenta')

#-- plot settings ------
xmin = 35
xmax = 85
ymin = 95
ymax = 170
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='both', which='major', labelsize=20)
ax.grid(True)
leg = ax.legend(loc='lower right', bbox_to_anchor=(0.95, 0.5), prop=fontP)

#-- labels/legend ------
ax.set_xlabel(r'$\theta_\mathrm{NMP}$', fontsize=24)
ax.set_ylabel(r'$\theta_\mathrm{LID}$', fontsize=24)
ax.set_title('NMP-CORE Angle vs LID-CORE Angle')

if showplots:
    plt.show()
if saveplots:
    fig.savefig("1ake4ake_nw_1-20_all_angles_2D"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# DMCD for 1ake and 4ake trajectories
#-----------------------
print "Plotting DCMD for 1ake and 4ake trajectories"

dmcd_1ake_full = np.loadtxt("1ake/1ake_nw_1-20_all_dmcd.txt")
dmcd_4ake_full = np.loadtxt("4ake/4ake_nw_1-20_all_dmcd.txt")
dmcd_1ake = dmcd_1ake_full[start:end:skip,:]
dmcd_4ake = dmcd_4ake_full[start:end:skip,:]
time = 240*np.arange(0, len(dmcd_4ake))/1000.

ref_dmcd = np.loadtxt("dmcd.txt")

#-- create figure ------
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
ax.set_aspect(0.4)

#-- plot data ----------
p1 = ax.scatter(dmcd_1ake[:,0], dmcd_1ake[:,1], c='red', s=ss, \
    edgecolor='none', linewidth=0.2, label='1ake', alpha=a1)
p2 = ax.scatter(dmcd_4ake[:,0], dmcd_4ake[:,1], c='green', s=ss, \
    edgecolor='none', linewidth=0.2, label='4ake', alpha=a2/2)

ax.scatter(ref_dmcd[2], ref_dmcd[3], marker="*", s=100, label="4ake", c='cyan')
ax.scatter(ref_dmcd[0], ref_dmcd[1], marker="x", s=80, label="1ake", c='magenta')

#-- plot settings ------
ax.tick_params(axis='both', which='major', labelsize=20)
ax.grid(True)
leg = ax.legend(loc='center right', bbox_to_anchor=(0.95, 0.5), prop=fontP)

#-- labels/legend ------
ax.set_xlabel(r'NMP-CORE MCD ($\AA$)', fontsize=24)
ax.set_ylabel(r'LID-CORE MCD ($\AA$)', fontsize=24)
ax.set_title('NMP-CORE vs LID-CORE Mass-Center Distances')

if showplots:
    plt.show()
if saveplots:
    fig.savefig("1ake4ake_nw_1-20_all_dmcd_2D"+postfix+".pdf")
plt.cla()
plt.clf()

#-----------------------
# Native contacts for 1ake and 4ake trajectories
#-----------------------
print "Plotting Native Contacts for 1ake and 4ake trajectories"

sc = 1.0
nc_1ake = np.loadtxt("1ake_nw_1-20_all_q1q2.dat")*sc
nc_4ake = np.loadtxt("4ake_nw_1-20_all_q1q2.dat")*sc
if subset or skip:
    nc_1ake = nc_1ake[start:end:skip,:]
    nc_4ake = nc_4ake[start:end:skip,:]
time = 240*np.arange(0, len(nc_4ake))/1000.

#-- create figure ------
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
ax.set_aspect(1)

#-- plot data ----------
p1 = ax.scatter(nc_1ake[:,0], nc_1ake[:,1], c='red', s=ss, \
    edgecolor='none', linewidth=0.2, label='1ake', alpha=a1)
p2 = ax.scatter(nc_4ake[:,0], nc_4ake[:,1], c='green', s=ss, \
    edgecolor='none', linewidth=0.2, label='4ake', alpha=a2/2)
# p1 = ax.hexbin(nc_1ake[:,0], nc_1ake[:,1], gridsize=10000, \
#     alpha=0.2)
# p2 = ax.hexbin(nc_4ake[:,0], nc_4ake[:,1], gridsize=1000)

#-- plot settings ------
xmin = ymin = 0.85*sc
xmax = ymax = 1.0*sc
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='both', which='major', labelsize=20)
ax.grid(True)
leg = ax.legend(loc='center right', bbox_to_anchor=(0.95, 0.5), prop=fontP)

#-- labels/legend ------
ax.set_ylabel('Fraction of 4ake Native Contacts', fontsize=24)
ax.set_xlabel('Fraction of 1ake Native Contacts', fontsize=24)
ax.set_title('2D Native Contacts Space')

if showplots:
    plt.show()
if saveplots:
    fig.savefig("1ake4ake_nw_1-20_all_q1q2_2D"+postfix+".pdf")
plt.cla()
plt.clf()


#################################################
# Secondary-structure analysis
#################################################
