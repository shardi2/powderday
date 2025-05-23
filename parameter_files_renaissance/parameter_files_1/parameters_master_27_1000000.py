# coding=utf-8
#===============================================
#HOME INFORMATION
#===============================================
pd_source_dir = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/'
#'/ufrc/narayanan/desika.narayanan/pd_git/'
import ytree
a = ytree.load("/storage/home/hcoda1/0/jw254/data/RS-RP/rockstar_halos/trees/tree_0_0_0.dat")

#===============================================
#RESOLUTION KEYWORDS
#===============================================
oref = 0 # over refine factor - should typically be set to 0
n_ref = 32 # when n_particles > n_ref, octree refines further
#zoom_box_len = 2*float(a["virial_radius"][0])
zoom_box_len = float(a["virial_radius"][28].to("kpc")*a["scale"][28])# kpc; so the box will be +/- zoom_box_len from the center
bbox_lim = 1.e5 # kpc - this is the initial bounding box of the grid (+/- bbox_lim)
                # This *must* encompass all of the particles in the
                # simulation. 

#===============================================
#PARALLELIZATION
#===============================================

n_processes = 8 # number of pool processes to run for stellar SED generation
n_MPI_processes = 8
		     # number of MPI tasks to run. for TORQUE this is
                     # best set as the same as n_processes, while for SLURM this may not be the case.

#===============================================
#RT INFORMATION
#===============================================
n_photons_initial = 1000000
n_photons_imaging = 1000000
n_photons_raytracing_sources = 1000000
n_photons_raytracing_dust = 1000000
n_photons_DIG = 1000000

FORCE_RANDOM_SEED = False
seed = -12345 # has to be an int, and negative.

#===============================================
#DUST INFORMATION 
#===============================================
#<<<<<<< HEAD
dustdir = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/hyperion-dust/dust_files/'
#'/home/desika.narayanan/hyperion-dust-0.1.0/dust_files/' #location of your dust files
#=======
#----------------
#DUST EXTINCTION
#----------------

#dustdir = '/home/desika.narayanan/hyperion-dust-0.1.0/dust_files/' #location of your dust files
#>>>>>>> 84e263f43cf8c5b736982df50d6317463866cdf2
dustfile = 'd03_3.1_6.0_A.hdf5'
PAH = True

dust_grid_type = 'dtm' # needs to be in ['dtm','rr','manual','li_bestfit','li_ml']
dusttometals_ratio = 0.4
enforce_energy_range = False # False is the default;  ensures energy conservation

SUBLIMATION = False # do we automatically kill dust grains above the
                    # sublimation temperature; right now is set to fast
                    # mode 
SUBLIMATION_TEMPERATURE = 1600. #K -- meaningliess if SUBLIMATION == False


#---------------------------------------------------------------
#Experimental Dust -- Note, these features are not fully vetted
#---------------------------------------------------------------

otf_extinction = False #flag for on the fly extinction.  If set, then we
                                    #ignore the dustdir/dustfile extinction information above. if
                                    #false, all otf_extinction* quantities are meaningless
otf_extinction_log_min_size = -4 #micron; must match what is set in the hydro simulation
otf_extinction_log_max_size = 0 #micron; must match what is set in the hydro simulation

percentile_LPAH_to_include = 0.99 #we run into memory problems if we
                                 #include every single PAH as a
                                 #source. however, you can usually
                                 #include 90-99% of the luminosity
                                 #without a problem and just ditch the
                                 #cells that don't contribute much.
                                 #lower this value if you run into
                                 #memory issues (though knowing of
                                 #course that you're missing PAH
                                 #flux...)

draine21_pah_model = True
draine21_pah_grid_write = True #this will write the PAH spectrum out
                               #to the npz file for every cell in the
                               #grid.  note, this causes the grid*.npz
                               #files to be ~10 GB or more.

draine_data_dir = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/dataverse_files/'

separate_into_dust_species = True #if set to true, we're assuming that
                                  #there are graphites, silicates, and
                                  #an aromatic fraction designated in
                                  #some way.  these are handled in
                                  #powderday/*_tributaries.

SKIP_LOGU_CALC = False #default is False -- if this is set, then we
                      #assume logU for the PAH calculations are all 0
                      #across the grid.  This can be set if there are
                      #memory troubles in calculating logU for massive
                      #grids (which can happen quickly).  Note -- logU
                      #only tends to impact things at logU>>4, which
                      #is hard to achieve in most galaxy evolution
                      #simulations (owing to resolution).


dust_density = 2.4 #units are g/cm**3.  this *only* is used in the OTF
                   #pah calculations, so if you do not have an active dust model, this will not be useful.
                               
#===============================================
#STELLAR SEDS INFO
#===============================================
FORCE_BINNED = True               # If True, force all star particles to be binned for calculating SED. 
                                  # If False, all star particles below max_age_direct (next parameter) are added 
                                  # directly without binning for calculating SED
max_age_direct  = 1.e-2           # Age (in Gyr) below which stars will be directly added without binning (works only if FORCE_BINNED is False)
                                  # To force all the star particles to be added without binning set this to an age greater than the maximum 
                                  # stellar age of the galaxy (say 16 Gyr for example)



imf_type = 2 # FSPS imf types; 0 = salpeter, 1 = chabrier; 2 = kroupa; 3 and 4 (vandokkum/dave) not currently supported
imf1 = 1.3 # Logarithmic slope of the IMF over the range 0.08 < M < 0.5. Only used if imf_type=2. (Default: 1.3)
imf2 = 2.3 # Logarithmic slope of the IMF over the range 0.5 < M < 1.0. Only used if imf_type=2. (Default: 2.3)
imf3 = 2.3 # Logarithmic slope of the IMF over the range 1.0 < M < 120. Only used if imf_type=2. (Default: 2.3)

pagb = 1 # weight given to post agb stars# 1 is the default

add_agb_dust_model = False    # add circumstellar AGB dust model (100%); Villaume, Conroy & Jonson 2015

#<<<<<<< HEAD
#===============================================
#NEBULAR EMISSION INFO
#===============================================
#add_neb_emission = True    			    # add nebular line emission (under active development)

#use_cloudy_tables = True    			    # If True, CLOUDY look up tables (dev. by Nell Byler) will be used to calculate 
                            			    # nebular emission. If False, CLOUDY models are generated individually 
                            			    # for each young star particle (under active development). 
                            			    # Note: The lookup tables work only for stars particles below 10 Myr.  (Default: True)

#use_cmdf = False                            # If True, star particles that have mass greater than cmdf_mas_mass (defined below) are broken down using a 
                                            # cluster mass distribution function (cmdf) of the form dN/dM goes as M^(beta). This works irrespecitve of whether
                                            # nebular emission is turned on or not.  The cmdf is set by the following parameters defined below: 
                                            # cmdf_min_mass, cmdf_max_mass, cmdf_bins and cmdf_beta.

#cmdf_min_mass = 3.5                         # Minimum mass of the star clusters in units of log(Msun). Note: Results might be inconsistent if
                                            # set lower than 3.5. (See Chandar et al.2014 for more info) (Default = 3.5)

#cmdf_max_mass = 5.0                         # Maximum mass of the star clusters in units of log(Msun). (Default = 5.0). Note: Only star particles that
                                            # have a mass greater than this parameter are broken down. 

#cmdf_bins = 6                               # The number of bins used for calulating the cluster mass distribution function (Default = 6.0)

#cmdf_beta = -2.0                            # Beta (power law exponent) for calculating CMDF (dN/dM goes as M^(beta)) 

#use_age_distribution = False                # Setting this to True, divides the star particles with ages between age_dist_min and age_dist_max (next parameters) into 
                                            # an ensemble of particles all of whom have the same properties except their age which is picked from a power law age 
                                            # distribution of the form dN/dt is proportional to t^-0.65 (Imp: This can only be used if use_cmdf is also set to True). 
                                            # Note: The function has a bunch of tunable parameters that can be changed though we feel that their default values
                                            # should be good enough for most cases. The function is located in cloudy_tools.py file under powderday/nebular_emission.

#age_dist_min = 3e-3                         # Star particle above this age are sub-divided into an age distribution if use_age_distribution is set to True
                                            # (Units: Gyr, Default = 3.e-3)

#age_dist_max = 1e-2                         # Star particles below this age are sub-divided into an age distribution if use_age_distribution is set to True
                                            # (Units: Gyr, Default = 1.e-2)

#***********************
# COMMON PARAMETERS
#***********************
# NOTE: These parmeters take either three or four values as an input. 
# They correspond to the value of the pararmeter for young_stars, PAGB stars, AGN and DIG respectively.

#FORCE_gas_logu = [False, False, False] 	    # If set, then we force the ionization parameter (gas_logu) to be 
                            			    # gas_logu (next parameter) else, it is taken to be variable and dependent on ionizing 
                            			    # radiation from star particles. (Default: [False, False, False])

#gas_logu = [-2.0, -2.0, -2.0]        		# Gas ionization parameter. This is only relevant 
                            			    # if add_neb_emission is set to True and FORCE_gas_logu is set to True (Default: [-2.0, -2.0, -2.0])

#gas_logu_init = [0.0, 0.0, 0.0]        	    # Force the ionization parameter to increase/decrease by this value (Scale: log). 
                            		    	# Useful if you want to run tests (Default: [0.0, 0.0, 0.0])

#FORCE_gas_logz = [False, False, False]      # If set, then we force the metallicity (gas_logz) to be gas_logz (next parameter)
                            	            # else, it is taken to be the star particles metallicity. (Default: [False, False, False])

#gas_logz = [0.0, 0.0, 0.0]  			    # Metallicity of the HII region in units of log(Z/Z_sun)
                            			    # only relevant if add_neb_emission = True and FORCE_gas_logz = True (Default: [0.0, 0.0, 0.0])

#FORCE_logq = [False, False, False]      	# If set, then we force the number of ionizing photons to be source_logq (next parameter)
                                            # else, it is taken to be variable and dependent on ionizing radiation of the source. (Default: [False, False, False])

#source_logq = [1.e47, 1.e47,1.e47]          # The number of ionizing photons emitted by the source in units of s^-1. Only relevant if add_neb_emission = True, 
    										# use_cloudy_tables = True and  FORCE_gas_logq = True (Default: [1.e47,1.e47,1.e47])  
                                            
#FORCE_inner_radius = [False, False, True]   # If set, then we force the inner radius of the cloud to be inner_radius (next parameter). 
											# IMP Note: This works only for young stars and Post-AGB stars. 
    										# For AGN we keep the inner radius fixed at whatever is set by inner_radius (next parameter) 
    										# irrespective of what this parameter is set to. (Default: [False,False,True])

#inner_radius = [1.e19, 1.e19, 2.777e+20]   	# This sets the inner radius of the cloud in cm. This is used only when add_neb_emission = True,
                            		    	# use_cloudy_tables = False and FORCE_inner_radius = True (Default: [1.e19, 1.e19, 2.777e+20], Units = cm)

#FORCE_N_O_Pilyugin = [False, False, False, False]  # If set to True, Nitrogen abundances are set according to the N/O vs O/H relation from Pilyugin et al. 2012
                                                   # If FORCE_N_O ratio (next parameter) is set to True then this parameter is ignored.(Default: [False,False,False, False])

#FORCE_N_O_ratio = [False, False, False, False]     # If set, then we force the Nitrogen abundance such that the log of N/O ratio is N_O_ratio (next parameter). 
                            			           # This can be used as a template fix adundance ratio of other elements (Default:  [False, False, False])

#N_O_ratio = [-0.85, -0.85, -0.85, -0.85]           # This sets the log of N/O ratio. This is used only when add_neb_emission = True,
                            			           # use_cloudy_tables = False, FORCE_N/O ratio = True and neb_abund = "direct" (Default: = [-0.85, -0.85, -0.85])

#neb_abund = ["dopita", "dopita", "dopita", "dopita"]  # This sets the HII region elemental abundances for generating CLOUDY models. 
                            			              # Available abundaces are.
                            			              #    dopita:    Abundances from Dopita (2001) with old solar abundances = 0.019 and ISM grains.
                            			              #    newdopita: Abundances from Dopita (2013). Solar Abundances from Grevasse 2010 - z= 0.013
                            			              #               includes smooth polynomial for N/O, C/O relationship functional form for He(z),
                            			              #               new depletion and factors in ISM grains.
                            			              #    gutkin:    Abundabces from Gutkin (2016) and PARSEC metallicity (Bressan+2012) based on Grevesse+Sauvel (1998) 
                            			              #               and Caffau+2011 
                            			              #    direct:    Abundances are taken directly from the simulation if possible. Defaults to using "dopita" if there is 
                            			              #               an error. (Note: Works only for AGNs and star particles that are added directly without binning.
                            			              #               Make sure to set FORCE_BINNED to False)
                            			              # This is used only when add_neb_emission = True and use_cloudy_tables = False. (Default: ["dopita", "dopita", "dopita"])
#***************************
# YOUNG STARS (HII regions)
#***************************

#add_young_stars = True     			        # If set, the young stars are included when calculating nebular emission (Default: True)


#HII_Rinner_per_Rs = 0.01        		    # Rinner for cloudy calculations is set to this value times the Stromgen Radius. 
                            			    # For example, if set to 0.01 Rinner is taken to be 1 % of Stromgren Radius. 
                            		    	# If FORCE_inner_radius (next parameter) is set to True then this is overridden
                            			    # and the value set by the inner_radius is used. This parameter is used 
                            			    # only when add_neb_emission = True and use_cloudy_tables = False (Default: 0.01)
   
#HII_nh = 1.e2               			    # Gas hydrogen density for calcualting nebular emission in units if cm^-3. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e2)

#HII_min_age = 1.e-3                         # Sets the minimum age limit for calculating nebular emission in units of Gyr. 
                                            # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e-3)

#HII_max_age = 1.e-2         			    # Sets the maximum age limit for calculating nebular emission in units of Gyr. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e-2)

#HII_escape_fraction = 0.0   			    # Fraction of H-ionizaing photons that escape the HII region. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 0.0)

#HII_alpha_enhacement = False                # If set, then the metallicity of star particles is set to [Fe/H] rather than the total metals. 
#=======
#alpha_enhacement = False                    # If set, then the metallicity of star particles is set to [Fe/H] rather than the total metals. 
#>>>>>>> 84e263f43cf8c5b736982df50d6317463866cdf2
                                            # Since FSPS does not support non solar abundance ratios, this parameter can be used to mimic the 
                                            # hardening of the radiation field due to alpha-enhancement. (Default: False)


#===============================================
#NEBULAR EMISSION
#===============================================
add_neb_emission = True                    # add nebular line emission (under active development)

use_cloudy_tables = False                    # If True, CLOUDY look up tables (dev. by Nell Byler) will be used to calculate 
                                            # nebular emission. Note that these tables are only applicable for star particles with ages below 10 Myr.
                                            # If False, CLOUDY models are generated individually for each star particle. This mode provides greater flexibility 
                                            # and allows for customization by tuning various free parameters. To take advantage of this flexibility, you need to 
                                            # provide a separate parameter master file (see parameters_master_neb.py example file provided with the code) when running 
                                            # Powderday. If you don't provide this file, default values will be used instead.

#===============================================
#BIRTH CLOUD INFORMATION
#===============================================

CF_on = False               # if set to true, then we enable the Charlot & Fall birthcloud models 

birth_cloud_clearing_age = 0.01 # Gyr - stars with age <
                                # birth_cloud_clearing_age have
                                # charlot&fall birthclouds meaningless
                                # of CF_on  == False


#===============================================
# Idealized Galaxy SED Parameters
#===============================================
Z_init = 0 # force a metallicity increase in the newstar particles.
           # This is useful for idealized galaxies.  The units for this
           # are absolute (so enter 0.02 for solar).  Setting to 0
           # means you use the stellar metallicities as they come in
           # the simulation (more likely appropriate for cosmological
           # runs)

           #NOTE - this is not exclusively used for idealized
           #simulations (i.e. one could use this for a cosmological
           #simulation), but the typical use case is for idealized simulations.

disk_stars_age = 8      # Gyr ;meaningless if this is a cosmological simulation
bulge_stars_age = 8     # Gyr ; meaningless if this is a cosmological simulation
disk_stars_metals = 19  # in fsps metallicity units
bulge_stars_metals = 19 # in fsps metallicity units



#===============================================
# Stellar Ages and Metallicities
#===============================================

# bins for binning the stellar ages and metallicities for SED
# assignments in cases of many (where many ==
# >N_METALLICITY_BINS*N_STELLAR_AGE_BINS) stars; this is necessary for
# reduction of memory load; see manual for details.

N_STELLAR_AGE_BINS = 100

#===============================================
#BLACK HOLES
#===============================================

BH_SED = False
BH_eta = 0.1 #bhluminosity = BH_eta * mdot * c**2.
BH_model = "Nenkova"
BH_modelfile = "/home/desika.narayanan/pd_git/agn_models/clumpy_models_201410_tvavg.hdf5"
# The Nenkova BH_modelfile can be downloaded here:
# https://www.clumpy.org/downloads/clumpy_models_201410_tvavg.hdf5
BH_var = True # Include time variations on BH luminosity (default Hickox+ 2014)

nenkova_params = [5,30,0,1.5,30,40] #Nenkova+ (2008) model parameters

#===============================================
#IMAGES AND SED PARAMETERS
#===============================================

NTHETA = 10
NPHI = 10
SED = True

SED_MONOCHROMATIC = False
FIX_SED_MONOCHROMATIC_WAVELENGTHS = False # if set, then we only use
                                          # the wavelengths in the
                                          # range between min_lam and
                                          # max_lam
SED_MONOCHROMATIC_min_lam = 0.0 # micron
SED_MONOCHROMATIC_max_lam = 1 #1000   # micron





IMAGING = True
filterdir = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/mean_throughputs_downsized/'
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/'
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/mean_throughputs_downsized/'
filterfiles = [#'one.filter', 'two.filter'
#    'arbitrary.filter',
#    'ACS_F475W.filter',
#    'ACS_F606W.filter',
#    'ACS_F814W.filter',
#    'B_subaru.filter',
     'F070W_downsized.txt',
#     'F300M_mean_system_throughput.txt',
     'F090W_downsized.txt',
#     'F322W2_mean_system_throughput.txt',
     'F115W_downsized.txt',
#     'F323N_mean_system_throughput.txt',
#     'F140M_mean_system_throughput.txt',
#     'F335M_mean_system_throughput.txt',
#     'F150W2_mean_system_throughput.txt',
     'F356W_downsized.txt',
     'F150W_downsized.txt',
#     'F360M_mean_system_throughput.txt',
#     'F162M_mean_system_throughput.txt',
#     'F405N_mean_system_throughput.txt',
#     'F164N_mean_system_throughput.txt',
#     'F410M_mean_system_throughput.txt',
#     'F182M_mean_system_throughput.txt',
#     'F430M_mean_system_throughput.txt',
#     'F187N_mean_system_throughput.txt',
     'F444W_downsized.txt',
     'F200W_downsized.txt',
#     'F460M_mean_system_throughput.txt',
#     'F210M_mean_system_throughput.txt',
#     'F466N_mean_system_throughput.txt',
#     'F212N_mean_system_throughput.txt',
#     'F470N_mean_system_throughput.txt',
#     'F250M_mean_system_throughput.txt',
#     'F480M_mean_system_throughput.txt',
     'F277W_downsized.txt'
#     'WLP4_mean_system_throughput.txt',
#     'F110W_downsized.txt'
]

# Insert additional filter files as above. In bash, the following command 
# formats the filenames for easy copying/pasting.
# $ shopt -s globstar; printf "#    '%s'\n" *.filter


npix_x = 128
npix_y = 128

#experimental and under development - not advised for use
IMAGING_TRANSMISSION_FILTER = False
filter_list = ['filters/irac_ch1.filter']
TRANSMISSION_FILTER_REDSHIFT = 0.001

#===============================================
#GRID INFORMATION  
#===============================================

MANUAL_ORIENTATION = False
THETA = 0
PHI = 0

#===============================================
#OTHER INFORMATION
#===============================================

PAH_frac = {'usg': 0.0586, 'vsg': 0.1351, 'big': 0.8063} # values will be normalized to 1

#===============================================
#DEBUGGING -THE PERFORMANCE OF THE CODE USING THESE PARAMETERS IS NOT GUARANTEED
#===============================================
SOURCES_RANDOM_POSITIONS = False
SOURCES_IN_CENTER = False
STELLAR_SED_WRITE = True
SKIP_RT = False # skip radiative transfer (i.e. just read in the grids and maybe write some diagnostics)
SUPER_SIMPLE_SED = False # just generate 1 oct of 100 pc on a side,
                         # centered on [0,0,0].  sources are added at
                         # random positions.
SKIP_GRID_READIN = False

CONSTANT_DUST_GRID = False # if set, then we don't create a dust grid by
                          # smoothing, but rather just make it the same
                          # size as the octree with a constant value of
                          # 4e-20
                          
N_MASS_BINS = 1 # this is really just a place holder that exists in
                # some loops to be able to insert some code downstream
                # for spatially varying IMFs.  right now for speed best
                # to set to 1 as it doesn't actually do anything.

FORCE_STELLAR_AGES = False
FORCE_STELLAR_AGES_VALUE = 0.05# Gyr

FORCE_STELLAR_METALLICITIES = False
FORCE_STELLAR_METALLICITIES_VALUE = 0.013 # absolute values (so 0.013 ~ solar) 

SKIRT_DATA_DUMP = True # if set, various data files useful for running SKIRT are saved.

REMOVE_INPUT_SEDS = False # If set, the hyperion input SED file is deleted after RT finishes

NEB_DEBUG = False # Dumps parameters related to nebular line emission in a file for debugging.
                  # The file includes the ionization parameter, number of ionizing photons,
                  # metallicity, inner radius, stellar mass and age for each particle.
                  # Naming convention: nebular_properties_galaxy*.txt where * is the galaxy number
OTF_EXTINCTION_MRN_FORCE = False
