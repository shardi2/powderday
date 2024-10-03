import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants
from hyperion.util.constants import pc 
import h5py
import os

#========================================================
#MODIFIABLE HEADER (make this a function later with argv)
z = 11.181356517874
run = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/example.0000.rtout.sed'
#'/home/desika.narayanan/pd/examples/gadget/mw_zoom/example.135.rtout.sed'
f = h5py.File("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/final_dataset_halo_catalog.hdf5", "a")
#========================================================






fig = plt.figure()
ax = fig.add_subplot(1,1,1)


m = ModelOutput(run)
wav,flux = m.get_sed(inclination='all',aperture=-1, distance = 10 * pc, units = 'mJy')

flux *= u.mJy

#print(flux)
wav  = np.asarray(wav)*u.micron #wav is in micron
wav *= (1.+z)

width = .01 * (1.+z)    # width of top-hat filter in microns
lam_UV = .15 * (1.+z) # UV wavelength in microns

min_lam_UV = np.float64(lam_UV - 0.5 * width)
max_lam_UV = np.float64(lam_UV + 0.5 * width)
min_lam_UV = min_lam_UV * u.micron
max_lam_UV = max_lam_UV * u.micron

flipped_wav = wav[::-1]

min_idx_UV = np.searchsorted(flipped_wav, min_lam_UV, side="left")
max_idx_UV = np.searchsorted(flipped_wav, max_lam_UV, side="left")

ab_mag_UV = np.empty(flux.shape[0])

for i in range(flux.shape[0]):
    flipped_flux = flux[i][::-1]

    sed_clip = np.trapz(flipped_flux[min_idx_UV : max_idx_UV].copy())

    #zero_point = 3.63078e6 * np.trapz(np.ones(len(sed_clip))) * u.mJy

    zero_point = 3.63078e6 * u.mJy
    
    absolute_magnitude_UV = -2.5 * np.log10(sed_clip/zero_point)

    ab_mag_UV[i] = absolute_magnitude_UV

    print(absolute_magnitude_UV)


if 'halo_0/photons_1e4/absolute_magnitude_UV' in f:
    del f['halo_0/photons_1e4/absolute_magnitude_UV']
    f.create_dataset('halo_0/photons_1e4/absolute_magnitude_UV', data=ab_mag_UV)
else:
    f.create_dataset('halo_0/photons_1e4/absolute_magnitude_UV', data=ab_mag_UV)

f.close()
