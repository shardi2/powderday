#!/bin/bash
#SBATCH -JSlurmshardinpowderdayrun                  # Job name
#SBATCH --account=gts-jw254-coda20                 # charge account
#SBATCH -N1 --ntasks-per-node=8                 # Number of nodes and cores per node required
#SBATCH --mem=700G                        # Memory per core
#SBATCH -p cpu-large 
#SBATCH -t200:00:00                               # Duration of the job (Ex: 1 hour)
#SBATCH -qinferno                               # QOS Name
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL              # Mail preferences
#SBATCH --mail-user=shardin31@gatech.edu        # E-mail address for notifications

cd /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday

module load anaconda3/2022.05.0.1
conda activate pd_environment_3
module load gcc
module load openmpi

python pd_front_end.py /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_renaissance/parameter_files_no_tables parameters_master_0_1000000 parameters_master_neb parameters_model_0_1000000
