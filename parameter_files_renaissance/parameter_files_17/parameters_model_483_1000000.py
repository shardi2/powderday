#snapshot parameters
import ytree
import numpy as np

snapshot_num = 41
output_num = 661
snapnum_str = '{:04d}'.format(snapshot_num)
output_str = '{:04d}'.format(output_num)

hydro_dir = '/storage/home/hcoda1/0/jw254/data/RS-RP/RD0041/'
#'/storage/home/hcoda1/0/jw254/data/SG64-2020/DD0125/'
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/'
#'/Volumes/pegasus/gadgetruns/SIGS/G2/'

#snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'
snapshot_name = 'RedshiftOutput'+snapnum_str

#where the files should go
PD_output_dir = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_renaissance_1000000/halo_661/'
#'/Volumes/pegasus/pd_runs/test/'
Auto_TF_file = 'snap'+snapnum_str+'.logical'
Auto_dustdens_file = 'snap'+snapnum_str+'.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/example.'+output_str+'.rtin'
outputfile = PD_output_dir+'/example.'+output_str+'.rtout'


#===============================================
#GRID POSITIONS
#===============================================
a = ytree.load("/storage/home/hcoda1/0/jw254/data/RS-RP/rockstar_halos/trees/tree_0_0_0.dat")

x_cent = float(a["position"][661][0])
y_cent = float(a["position"][661][1])
z_cent = float(a["position"][661][2])

TCMB = 2.73
