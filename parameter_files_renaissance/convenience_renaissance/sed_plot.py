import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants
import h5py
import os

#========================================================
#MODIFIABLE HEADER (make this a function later with argv)
z = 11.181356517874
#ds = yt.load('/storage/home/hcoda1/0/jw254/data/SG64-2020/DD0125/output_0125') 
#run = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000000/halo_0/example.0000.rtout.sed'
run  = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_renaissance_1000000/halo_0/example.0000.rtout.sed'
#'/home/desika.narayanan/pd/examples/gadget/mw_zoom/example.135.rtout.sed'
#f = h5py.File("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/parameter_files_renaissance/final_dataset_halo_catalog.hdf5", "a")
filter_directory = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/mean_throughputs_downsized'
#========================================================







fig = plt.figure()
ax = fig.add_subplot(1,1,1)


m = ModelOutput(run)
wav,flux = m.get_sed(inclination=0,aperture=-1)

wav  = np.asarray(wav)*u.micron #wav is in micron
wav *= (1.+z)

flux = np.asarray(flux)*u.erg/u.s
dl = Planck13.luminosity_distance(z)
dl = dl.to(u.cm)
    
flux /= (4.*3.14*dl**2.)
    
nu = constants.c.cgs/(wav.to(u.cm))
nu = nu.to(u.Hz)

flux /= nu
flux = flux.to(u.mJy)

combined_array = np.column_stack((wav.to_value(), flux.to_value()))

#np.savetxt('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000000/pd_sed/sed_0000.txt', combined_array, fmt='%.15f')

#for i in range(flux.shape[0]):
#    ax.loglog(wav,flux[i,:])

ax.loglog(wav, flux)

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel('Flux (mJy)')
#ax.set_ylim([1,1e8])
#ax.set_xlim(0.05,15000)
ax.grid()

fig.savefig('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_renaissance_1000000/halo_0/sed_0000.png')

#if 'halo_0/photons_1e4/SED_wav' in f:
#    del f['halo_0/photons_1e4/SED_wav']
#    f.create_dataset('halo_0/photons_1e4/SED_wav', data=wav)
#else:
#    f.create_dataset('halo_0/photons_1e4/SED_wav', data=wav)

#if 'halo_0/photons_1e4/SED_flux' in f:
#    del f['halo_0/photons_1e4/SED_flux']
#    f.create_dataset('halo_0/photons_1e4/SED_flux', data=flux)
#else:
#    f.create_dataset('halo_0/photons_1e4/SED_flux', data=flux)

#f.close()
