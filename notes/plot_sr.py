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

data = np.genfromtxt('sr.txt')

def inversion(internal_state):
    return ((internal_state[:,2]**2 + internal_state[:,3]**2) -
            (internal_state[:,0]**2 + internal_state[:,1]**2))

internal_states=[np.genfromtxt('internal_state_0_' + str(i) + '.txt') for i in range(0, 1000, 10)]

invs = np.array([np.mean(inversion(s)) for s in internal_states])

dt = 9.954614e-09
plt.plot(1.0e6 * data[:,0] * dt, (data[:,3]**2 + data[:,2]**2) / 5.0e5)
plt.plot(1.0e6 * data[:,0] * dt, invs)
plt.xlabel(r'$t/\mu{\rm s}$')
plt.ylabel(r'$|a|^2$, $w$')
plt.gcf().subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('sr.pdf')

