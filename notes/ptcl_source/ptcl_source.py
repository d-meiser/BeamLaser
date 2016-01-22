import subprocess
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

beam_laser_tla_path = '../../build/examples/'
beam_laser_tla_exec = beam_laser_tla_path + 'BeamLaserTLA'

nbar = 1000
vbar = 1.0e3
dt = 1.0e-9
num_steps = 1000

output = subprocess.Popen([beam_laser_tla_exec,
                           "--nbar", str(nbar),
                           "--maxNumPtcl", str(1.5 * nbar),
                           "--vbar", str(vbar),
                           "--dt", "1.0e-9",
                           "--numSteps", str(num_steps),
                           "--dumpPeriod", "10",
                           "--dumpField",
                           "--dumpPhaseSpace"
], stdout=subprocess.PIPE).communicate()[0]

data = output.splitlines()
data = [d.split() for d in data]
data = [[float(f) for f in l] for l in data]
data = np.array(data)
print data[:,:2]
plt.plot(dt * data[:, 0]  / 1.0e-6, data[:, 1])
plt.xlabel(r'$t / \mu s$')
plt.ylabel(r'$N$')
plt.xlim([0, dt * data[-1, 0] / 1.0e-6])
plt.ylim([0, 1.2 * nbar])
plt.gcf().subplots_adjust(left = 0.18, bottom = 0.21, top=0.95, right=0.95)
plt.savefig('numPtclsVsT.pdf')

phase_space = np.loadtxt('ptcls_0_900.txt')

plt.clf()
plt.plot(phase_space[0, :] / 1.0e-6, phase_space[1,:] / 1.0e-6, '.')
plt.xlabel(r'$x/\mu {\rm m}$')
plt.ylabel(r'$y/\mu {\rm m}$')
plt.xlim([-1.2e2, 1.2e2])
plt.ylim([-1.2e2, 1.2e2])
plt.gcf().subplots_adjust(left = 0.2, bottom = 0.21, top=0.95, right=0.95)
plt.savefig('density.pdf')

plt.clf()
plt.plot(phase_space[3, :], phase_space[5,:], '.')
plt.xlabel(r'$v_x / ({\rm m} / {\rm s})$')
plt.ylabel(r'$v_z / ({\rm m} / {\rm s})$')
plt.gcf().subplots_adjust(left = 0.22, bottom = 0.21, top=0.95, right=0.95)
plt.savefig('velocity_distribution.pdf')
