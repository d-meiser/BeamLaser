import subprocess
import numpy as np
import os

top_dir = os.getcwd()

beam_laser_tla_path = top_dir + '/../../build-mpi/examples/'
beam_laser_tla_exec = beam_laser_tla_path + 'BeamLaserTLA'
mpi_exec = top_dir + '/../../mpich/bin/mpirun'

nbar_base = 10
vbar = 1.0e-1
deltaV = 1.0e-4
dt = 1.0e-7
dump_period = 1.0e4
num_steps = 50001
kappa = 1.0e8
dipole_matrix_elemen = 1.0e-31
ptcl_weight = 1.0e0
waist = 3.0e-5

for i in range(6):
    run_dir = top_dir + '/run' + str(i)
    try:
        os.mkdir(run_dir)
    except OSError:
        print "Failed to create directory -- assuming it already exists."

    os.chdir(run_dir)

    nbar = nbar_base * 10**i
    run_cmd = [
            mpi_exec,
            '-np', str(4),
            beam_laser_tla_exec,
            "--nbar", str(nbar),
            "--maxNumPtcl", str(100 * nbar),
            "--vbar", str(vbar),
            "--deltaV", str(deltaV),
            "--kappa", str(kappa),
            "--waist", str(waist),
            "--dt", str(dt),
            "--numSteps", str(num_steps),
            "--dumpPeriod", str(dump_period),
            "--dipoleMatrixElement", str(dipole_matrix_elemen),
            "--ptclWeight", str(ptcl_weight),
            "--dumpPhaseSpace",
            "--dumpField",
            "--dumpInternalState"
            ]
    log_file = open('run' + str(i) + '.txt', 'w')
    log_file.write("# Executing: \n")
    log_file.write("# " + " ".join(run_cmd))
    output = subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()[0]
    log_file.write(output)
    log_file.close()

    os.chdir(top_dir)

