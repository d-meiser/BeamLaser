import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
have_plt = True
fig_width_pt = 8.0*72/2.54
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*golden_mean
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': fig_size}
plt.rcParams.update(params)

internal_state = []
phase_space = []

i_max = 6

for i in range(i_max):
    internal_state.append(
            np.concatenate([np.loadtxt('run' + str(i) + '/internal_state_' +
                    str(j) + '_2000.txt') for j in range(4)]))
    phase_space.append(
            np.concatenate([np.transpose(np.loadtxt('run' + str(i) + '/ptcls_' +
                    str(j) + '_2000.txt')) for j in range(4)]))

for i in range(i_max):
    polarization = np.abs(internal_state[i][:,0] * internal_state[i][:,2] +
            internal_state[i][:,1] * internal_state[i][:,3])**2 + np.abs(
                internal_state[i][:,1] * internal_state[i][:,2] +
                internal_state[i][:,0] * internal_state[i][:,3])**2 
    print np.mean(polarization)
    plt.plot(phase_space[i][:,2], polarization)
plt.savefig('polarizations.pdf')
