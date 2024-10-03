#!/bin/bash
#SBATCH -JSlurmshardinpowderdayrun                  # Job name
#SBATCH --account=gts-jw254                 # charge account
#SBATCH -N2 --ntasks-per-node=4                 # Number of nodes and cores per node required
#SBATCH --mem=120G                        # Memory per core
#SBATCH -t60:00:00                               # Duration of the job (Ex: 1 hour)
#SBATCH -qinferno                               # QOS Name
#SBATCH --array=1-9
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL              # Mail preferences
#SBATCH --mail-user=shardin31@gatech.edu        # E-mail address for notifications

cd /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/

module load anaconda3/2022.05.0.1
conda activate pd_environment_3

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	./run_all_parameters_10000.sh
elif [ $SLURM_ARRAY_TASK_ID -eq 2 ]; then
	./run_all_parameters_100000.sh
elif [ $SLURM_ARRAY_TASK_ID -eq 3 ]; then
	./run_all_parameters_1000000.sh
elif [ $SLURM_ARRAY_TASK_ID -eq 4 ]; then
	./run_all_parameters_10000000.sh

fi
