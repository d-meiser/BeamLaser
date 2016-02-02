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


data = []
for i in range(6):
    data.append(np.loadtxt('run' + str(i) + '/run' + str(i) + '.txt'))
data = np.array(data)
for i in range(6):
    plt.semilogy(1.0e-4 * data[i, :, 0], data[i, :,2]**2 + data[i, :,3]**2)
plt.xlim([0, 1.0e0])
plt.xlabel(r'$t/{\rm ms}$')
plt.ylabel(r'$|a|^2$')
plt.gcf().subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('field_intensities.pdf')
