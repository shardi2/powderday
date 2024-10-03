#!/bin/bash
#SBATCH -J test
#SBATCH --account=gts-jw254                     
#SBATCH -N1 --ntasks-per-node=2                
##SBATCH --mem-per-cpu=1G
#SBATCH -t0:01:000                              
#SBATCH -qinferno                               
#SBATCH -oReport-%j.out                         
##SBATCH --mail-type=BEGIN,END,FAIL              
##SBATCH --mail-user=cdbs3@gatech.edu        

echo START!
echo DONE!
