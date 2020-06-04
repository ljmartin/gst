import matplotlib.pyplot as plt
from seaborn import kdeplot
plt.style.use('ggplot')

import mdtraj as md
import pymc3 as pm
import theano.tensor as tt

import numpy as np

##the dihedral idices:
indices = np.array([[4, 6, 8, 14],[6, 8, 14, 16]])


##load standard MD run:
traj_standardmd = md.load_dcd('./trajectories/diala_standardmd_traj.dcd', top='./alanine-dipeptide-implicit.pdb')
dihedrals_standardmd = md.compute_dihedrals(traj_standardmd, indices, periodic=True)
dihedrals_standardmd[:,0] = np.where(dihedrals_standardmd[:,0]<2, dihedrals_standardmd[:,0], dihedrals_standardmd[:,1]-np.pi)

##load GST run:
traj_gst = md.load_dcd('./trajectories/diala_gst_traj.dcd', top='./alanine-dipeptide-implicit.pdb')
dihedrals_gst = md.compute_dihedrals(traj_gst, indices, periodic=True)
dihedrals_gst[:,0] = np.where(dihedrals_gst[:,0]<2, dihedrals_gst[:,0], dihedrals_gst[:,1]-np.pi)
dihedrals_gst = dihedrals_gst[::2]
print('Lengths: SMD:', len(dihedrals_standardmd), len(dihedrals_gst))

fig, ax = plt.subplots(2,1, sharey=True)
fig.set_figwidth(15)
fig.set_figheight(5)

ax[0].plot(dihedrals_standardmd[:,0], label='Standard MD')
ax[0].legend()
ax[1].plot(dihedrals_gst[:,0], c='C1', label='Serial tempering')
ax[1].legend()
for a in ax:
    a.axhline(np.pi, c='k', linestyle='--')
    a.axhline(-np.pi, c='k', linestyle='--')


fig.text(0.00, 0.5, 'Slow DoF (φ)', va='center', rotation='vertical')
fig.savefig('slow_dof.tif')
fig.savefig('slow_dof.svg')

####Now fit a simple 2-state markov state model

class markov_model(pm.distributions.Discrete):
    def __init__(self, p, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.p = k = tt.as_tensor_variable(p)
        self.mode = tt.as_tensor_variable(0.5)

    def logp(self, x):
        p = self.p
        x_im1 = x[:-1]
        x_i = x[1:]
        likelihood = pm.Bernoulli.dist(p[x_im1]).logp(x_i)

        return tt.sum(likelihood)


#binarize the slow dof (phi) into two states using a priori knowledge of the states:
smd= dihedrals_standardmd[:,0] < -0.25
gst = dihedrals_gst[:,0] < -0.25

#1-arr.astype(int) means 0=ground state, 1=excited state
with pm.Model() as mo:
    a = pm.Beta('a', alpha=1, beta=1)
    b = pm.Beta('b', alpha=1, beta=1)
    bah = markov_model('bah', p = [a,b], observed=1-smd.astype(int))    
    trace_smd = pm.sample(1000)

with pm.Model() as mo:
    a = pm.Beta('a', alpha=1, beta=1)
    b = pm.Beta('b', alpha=1, beta=1)    
    bah = markov_model('bah', p = [a,b], observed=1-gst.astype(int))
    trace_gst = pm.sample(1000)



fig = plt.figure(constrained_layout=True)

gs = fig.add_gridspec(2,20)
ax0 = fig.add_subplot(gs[0, :11])
ax1 = fig.add_subplot(gs[1, :11])
ax2 = fig.add_subplot(gs[:2, 11:])

fig.set_figwidth(15)
fig.set_figheight(5)

time = np.arange(10000)/10
ax0.plot(time, dihedrals_standardmd[:,0], label='Standard MD')
ax0.legend()
ax1.plot(time, dihedrals_gst[:,0], c='C1', label='Serial tempering')
ax1.legend()

for a in [ax0, ax1]:
    a.axhline(np.pi, c='k', linestyle='--')
    a.axhline(-np.pi, c='k', linestyle='--')
    a.set_ylabel('φ')
    a.set_xlabel('Time (ns)')

#doesnt work with contsrained layout:
#fig.text(0.04, 0.5, 'Slow DoF (φ)', va='center', rotation='vertical')



for trace, col, label in zip([trace_smd, trace_gst], ['C0', 'C1'], ['Standard MD', 'Serial tempering']):
    samples = trace['a']
    kdeplot(samples, color=col, label=label)
    hpd = pm.stats.hpd(samples)
    #plt.plot( hpd,[-0.5,-0.5],linewidth=2, c=col)
    m = samples.mean()
    print(hpd)
    ax2.errorbar(m, -0.5, xerr = np.array([m-hpd[0], hpd[1]-m])[:,None],mfc='white',mew=2,
                           fmt='o',c=col, linewidth=4, markersize=15, capsize=3)

ax2.set_ylabel('Probability density function')
ax2.set_yticks([])
ax2.set_xlabel('Transition probability')
ax2.set_title('Density')
ax2.legend()


fig.savefig('all.tif')
fig.savefig('all.svg')
fig.savefig('all.png')
