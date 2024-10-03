#!/bin/bash
#SBATCH -JSlurmshardinpowderdayrun                  # Job name
#SBATCH --account=gts-jw254                 # charge account
#SBATCH -N2 --ntasks-per-node=4                 # Number of nodes and cores per node required
#SBATCH --mem=120G                        # Memory per core
#SBATCH -t60:00:00                               # Duration of the job (Ex: 1 hour)
#SBATCH -qinferno                               # QOS Name
#SBATCH -oReport-%j.out                         # Combined output and error messages file
#SBATCH --mail-type=BEGIN,END,FAIL              # Mail preferences
#SBATCH --mail-user=shardin31@gatech.edu        # E-mail address for notifications

cd /storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/convenience

module load anaconda3/2022.05.0.1
conda activate pd_environment_3

./run_image.sh
