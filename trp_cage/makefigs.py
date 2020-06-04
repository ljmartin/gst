import matplotlib.pyplot as plt
plt.style.use('ggplot')

import mdtraj as md
import numpy as np
reference_pdb = md.load_pdb('./nmr_struc_1.pdb')


standard_traj = md.load_dcd('./trajectories/trp_standardmd_traj.dcd', top='./trp-min.pdb')
#gst_equil_traj = md.load_dcd('./trp_gst_equilibration.dcd', top='./trp-min.pdb')
gst_traj =  md.load_dcd('./trajectories/trp_gst_traj.dcd', top='./trp-min.pdb')


rmsd_smd = md.rmsd(standard_traj, reference_pdb)
#rmsd = md.rmsd(gst_equil_traj, reference_pdb)
rmsd_gst = md.rmsd(gst_traj, reference_pdb)


fig, ax = plt.subplots()
fig.set_figheight(5)
fig.set_figwidth(10)

time = np.linspace(0, len(rmsd_smd)*50000/500000,len(rmsd_smd) )
ax.plot(time, rmsd_gst, label='Serial tempering')
ax.plot(time, rmsd_smd, label='Standard MD')
ax.set_ylabel('RMSD')
ax.set_xlabel('Time (ns)')
ax.legend()


fig.savefig('../figures/trpcage.tif')
fig.savefig('../figures/trpcage.svg')
fig.savefig('../figures/trpcage.png')


standard_traj[np.argmin(rmsd_smd)].save_pdb('./standard_minimum.pdb')
gst_traj[np.argmin(rmsd_gst)].save_pdb('./gst_minimum.pdb')
