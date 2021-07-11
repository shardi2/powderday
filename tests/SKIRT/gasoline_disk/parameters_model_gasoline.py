#snapshot parameters

snapnum_str = '00300'
hydro_dir = '/blue/narayanan/desika.narayanan/yt_datasets/TipsyGalaxy/'

galaxy_num = 0


galaxy_num_str = str(galaxy_num)

snapshot_name = 'galaxy.00300'

#where the files should go
PD_output_dir = '/home/desika.narayanan/pd_git/tests/SKIRT/gasoline_disk/'
Auto_TF_file = 'gasoline_test.logical'
Auto_dustdens_file = 'gasoline_test.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/pd_skirt_comparison.gasoline.rtin'
outputfile = PD_output_dir+'/pd_skirt_comparison.gasoline.rtout'

#===============================================
#GRID POSITIONS
#===============================================
x_cent=0.
y_cent=0.
z_cent=0.

#===============================================
#CMB
#===============================================
TCMB = 2.73
