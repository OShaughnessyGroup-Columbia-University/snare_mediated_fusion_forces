#!/bin/sh
#
# Simple "fusion" submit script for Slurm.
#
#SBATCH --account=oshaughnessy # The account name for the job.
#SBATCH --job-name=rtmd_eq # The job name.
#SBATCH --array=5 # The job ID.This labels the run number.
# #SBATCH -c 4 # The number of cpu cores to use.
#SBATCH --time=12:00:00 # The time the job will take to run.
#SBATCH --mem-per-cpu=10gb # The memory the job will use per cpu core.
#SBATCH --gres=gpu:1
#SBATCH --constraint=a100
 
module load singularity
 
#Command to execute Python program:
singularity exec --nv /burg/oshaughnessy/users/icb2114/hoomd_old.sif python3 ves_ves_nozip_snare_wlc_ini_rad_eq_movavg.py $SLURM_ARRAY_TASK_ID 0 0 0 6 1 4 

# argv 1 is the run number
# argv 2, 3, 4 are always 0 for SNARE simulations
# argv 5 controls the number of rods: n_rod_grp = [12,9,7,6,5,4,3]
# argv 6 controls the number of unzipped LD residues: N_unzip_grp = [7,10,14,18]
# argv 7 controls the LD persistence length: lp_grp = [0.3,0.35,0.4,0.45,0.5,0.55,0.6]
# NOT USED HERE argv 8 controls the initial ring radius: r_eci_arr = [9.5,8.5,7.5,6.5,5.5,4.5,3.5]

#End of script
