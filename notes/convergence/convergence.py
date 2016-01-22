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
num_steps = 401

phase_space = []
for i in range(8):
    print dt / (2**i)
    output = subprocess.Popen([beam_laser_tla_exec,
                               "--nbar", str(nbar),
                               "--maxNumPtcl", str(1.5 * nbar),
                               "--vbar", str(vbar),
                               "--dt", str(dt / (2**i)),
                               "--numSteps", str(num_steps * 2**i),
                               "--dumpPeriod", str(400 * 2**i),
                               "--dumpInternalState"
    ], stdout=subprocess.PIPE).communicate()[0]
    phase_space.append(np.loadtxt('internal_state_0_' + str(400 * 2**i) + '.txt'))

phase_space = np.array(phase_space)

print phase_space

