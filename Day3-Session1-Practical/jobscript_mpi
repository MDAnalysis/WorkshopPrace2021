#!/bin/bash

#SBATCH -J MPI_job
#SBATCH -n 4
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node 2
#SBATCH -N 2
#SBATCH -t 00:10:00
#SBATCH -p normal
#SBATCH --reservation=parallel_mda

## SBATCH flag explanation (read more under man sbatch)
# -J MPI_job            The job name.
# -n 4                  How many tasks (programs) are you planning on running. In this case, only one.
# --cpus-per-task 1     How many cpus will be assigned to the task.
# --ntasks-per-node 2   How many MPI ranks will run on each node.
# -N 2                  Over how many different nodes (computers) do we want to run our program. Must be only 1 when using multiprocessing, as it cannot communicate over nodes.
# -t 00:10:00           How long to reserve the computing resources for. Should the job exceed this requested time, the cluster may kill it.
# -p and --reservation  The cluster resources to use. Check your cluster's documentation to know which partition/reservation you should use.

# This file loads the workshop-related modules and Python environment (which uses anaconda).
# Consult your cluster's documentation to know which modules you should load in your case.
source /project/jhlsrf005/JHL_hooks/env

# Any text output will be collected in slurm .out log file
mpirun -np 4 ./mpi_example.py
