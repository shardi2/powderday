import matplotlib
matplotlib.use('Agg')
import pdb
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
run = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/halo_0/example.0000.rtout.sed'
#'/home/desika.narayanan/pd/examples/gadget/mw_zoom/example.135.rtout.sed'
f = h5py.File("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/final_dataset_halo_catalog.hdf5", "a")
filter_directory = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/mean_throughputs_downsized/'
#========================================================







fig = plt.figure()
ax = fig.add_subplot(1,1,1)


m = ModelOutput(run)
wav,flux = m.get_sed(inclination='all',aperture=-1)

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


micron  = []
throughput = []
x = 0

for i in range(flux.shape[0]):
    for filters in os.listdir(filter_directory):
        f = open(filter_directory + filters, "r")
        for line in f:
            row = line.split()
            micron.append(row[0])
            throughput.append(row[1])
        micron = np.asarray(micron)
        throughput = np.asarray(throughput)
        micron = micron.astype(np.float64) * (1+z)
        throughput = throughput.astype(np.float64)
        micron *= u.micron
        min_lam = max(wav.min(), micron.min())
        max_lam = min(wav.max(), micron.max())
        flipped_wav = wav[::-1]
        min_idx = np.searchsorted(flipped_wav, min_lam, side="left")
        max_idx = np.searchsorted(flipped_wav, max_lam, side="left")
        lam_clip = flipped_wav[min_idx : max_idx]
        throughput_interp = np.interp(lam_clip, micron, throughput)
        flipped_flux = flux[i][::-1]
        sed_clip = flipped_flux[min_idx : max_idx]
        print("sed clip", sed_clip)
        sed_final = np.trapz(sed_clip * throughput_interp)
        print("sed final", sed_final)
        zero_point = 3.63078e6 * np.trapz(throughput_interp) * u.mJy
        apparent_magnitude = -2.5 * np.log10(sed_final / zero_point)
        print("app mag", apparent_magnitude)
        del micron
        del throughput
        micron = []
        throughput = []
        x += 1


#if 'halo_0/photons_1e4/apparent_mag' in f:
#    del f['halo_0/photons_1e4/apparent_mag']
#    f.create_dataset('halo_0/photons_1e4/apparent_mag', data=wav)
#else:
#    f.create_dataset('halo_0/photons_1e4/apparent_mag', data=wav)

f.close() 
