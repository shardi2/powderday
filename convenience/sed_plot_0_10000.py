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
run = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/example.0000.rtout.sed'
#'/home/desika.narayanan/pd/examples/gadget/mw_zoom/example.135.rtout.sed'
#f = h5py.File("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/final_dataset_halo_catalog.hdf5", "a")
filter_directory = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/mean_throughputs_downsized'
#========================================================






fig = plt.figure()
ax = fig.add_subplot(1,1,1)


m = ModelOutput(run)
wav,flux = m.get_sed(inclination='all',aperture=-1)

wav  = np.asarray(wav)*u.micron #wav is in micron
wav *= (1.+z)

print(flux)

flux = np.asarray(flux)*u.erg/u.s
dl = Planck13.luminosity_distance(z)
dl = dl.to(u.cm)
    
flux /= (4.*3.14*dl**2.)

print(flux)

nu = constants.c.cgs/(wav.to(u.cm))
nu = nu.to(u.Hz)

flux /= nu
flux = flux.to(u.mJy)

print(flux)

width = .01    # width of top-hat filter in microns
lam_UV = .15*(1+z)  # UV wavelength in microns

min_lam_UV = np.float64(lam_UV - 0.5 * width)
max_lam_UV = np.float64(lam_UV + 0.5 * width)
min_lam_UV = min_lam_UV * u.micron
max_lam_UV = max_lam_UV * u.micron

sorted_wav = np.sort(wav)

min_idx_UV = np.searchsorted(sorted_wav, min_lam_UV, side="left")
max_idx_UV = np.searchsorted(sorted_wav, max_lam_UV, side="left")

#lam_clip = sorted_wav[min_idx_UV : max_idx_UV]
sed_clip = (flux[1][min_idx_UV : max_idx_UV].copy()) 

#print(lam_clip)
print(sed_clip)

#c = 2.99792458e14 * u.micron #in microns
#freq = c/lam_clip

#zero_point = 3.63078e6 * np.ones(len(sed_clip)) * (freq)**(-1) * u.mJy

zero_point = 3.63078e6 * np.ones(len(sed_clip)) * u.mJy

absolute_magnitude_UV = -2.5 * np.log(sed_clip/zero_point)

print(absolute_magnitude_UV)

for i in range(flux.shape[0]):
    ax.loglog(wav,flux[i,:])

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel('Flux (mJy)')
#ax.set_ylim([1,1e8])
#ax.set_xlim(0.05,15000)
ax.grid()

fig.savefig('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/pd_sed/sed_0000.png')

#f.create_dataset('halo_0/10000/SED_wav/', data=wav)

#f.create_dataset('halo_0/10000/SED_flux/', data=flux)

#f.create_dataset('halo_0/10000/UV_absolute_mag/', data=absolute_magnitude_UV)

#micron  = []
#throughput = []
#for filters in os.listdir(filter_directory):
#    f = open(filter_directory + filters, "r")
#    for line in f:
#        row = line.split()
#        micron.append(row[:-1])
#        throughput.append(row[-1])
#    min_lam = max(wav.min(), micron.min())
#    max_lam = min(wav.max(), micron.max())
#    min_idx = np.searchsorted(wavs, min_lam, side="left")
#    max_idx = np.searchsorted(wavs, max_lam, side="left")
#    lam_clip = wav[min_idx : max_idx]
#    throughput_interp = np.interp(lam_clip, micron, throughput)
#    flux *= throughput_interp
#    zero_point *= throughput_interp
#    micron.clear()
#    throughput.clear()
    
#apparent_magnitude = -2.5 * np.log10(flux / zero_point)

#f.create_dataset('halo_0/10000/apparent_magnitude/', data=apparent_magnitude)

#f.close()
