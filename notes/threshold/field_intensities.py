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

i_min = 0
i_max = 6

data = []
for i in range(i_min, i_max):
    data.append(np.loadtxt('run' + str(i) + '/run' + str(i) + '.txt'))
data = np.array(data)
for i in range(i_min, i_max):
    plt.semilogy(1.0e-3 * data[i - i_min, :, 0],
            data[i - i_min, :,2]**2 + data[i - i_min, :,3]**2)
plt.xlim([0, 1.0e1])
plt.xlabel(r'$t/{\rm ms}$')
plt.ylabel(r'$|a|^2$')
plt.gcf().subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('field_intensities.pdf')

plt.clf()
plt.plot(1.0e-3 * data[-1, :, 0], data[-1, :,2]**2 + data[-1, :,3]**2)
plt.xlim([0, 1.0e1])
plt.xlabel(r'$t/{\rm ms}$')
plt.ylabel(r'$|a|^2$')
plt.gcf().subplots_adjust(left=0.2, bottom=0.2)
plt.savefig('field_intensity_above_threshold.pdf')
