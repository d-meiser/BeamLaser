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

i_min = 0
i_max = 6

for i in range(i_min, i_max):
    internal_state.append(
            np.loadtxt('run' + str(i) + '/internal_state_' +
                    str(0) + '_20000.txt'))
    phase_space.append(
            np.transpose(np.loadtxt('run' + str(i) + '/ptcls_' +
                    str(0) + '_20000.txt')))

for i in reversed(range(i_min, i_max)):
    polarization = np.abs(internal_state[i-i_min][:,0] *
            internal_state[i-i_min][:,2] +
            internal_state[i-i_min][:,1] * internal_state[i-i_min][:,3])**2 + np.abs(
                internal_state[i-i_min][:,1] * internal_state[i-i_min][:,2] +
                internal_state[i-i_min][:,0] * internal_state[i-i_min][:,3])**2 
    print np.mean(polarization)
    plt.plot(1.0e6*phase_space[i-i_min][::100,2], polarization[::100],
            '.',ms=1)
plt.xlabel(r'$z/\mu{\rm m}$')
plt.ylabel(r'$P$')
plt.xlim([-100, 100])
plt.ylim([0, 0.3])
plt.gcf().subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('polarizations.pdf')
