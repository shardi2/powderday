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

import sys
np.set_printoptions(threshold=sys.maxsize)
from matplotlib.collections import LineCollection

#========================================================
#MODIFIABLE HEADER (make this a function later with argv)
z = 11.181356517874
#ds = yt.load('/storage/home/hcoda1/0/jw254/data/SG64-2020/DD0125/output_0125') 
#run_100 = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_100000/halo_0/example.0000.rtout.sed'
run_1mil = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_3000000/halo_1/example.0001.rtout.sed'
run_100 = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_300000/halo_1/example.0001.rtout.sed'
#'/home/desika.narayanan/pd/examples/gadget/mw_zoom/example.135.rtout.sed'
filter_directory = '/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/filters/mean_throughputs_downsized'
#========================================================







fig = plt.figure()
ax = fig.add_subplot(1,1,1)


m_1mil = ModelOutput(run_1mil)
wav,flux_1mil = m_1mil.get_sed(inclination=0,aperture=-1)

wav  = np.asarray(wav)*u.micron #wav is in micron
wav *= (1.+z)

flux_1mil = np.asarray(flux_1mil)*u.erg/u.s
dl = Planck13.luminosity_distance(z)
dl = dl.to(u.cm)
    
flux_1mil /= (4.*3.14*dl**2.)
    
nu = constants.c.cgs/(wav.to(u.cm))
nu = nu.to(u.Hz)

flux_1mil /= nu
flux_1mil = flux_1mil.to(u.mJy)


m_100 = ModelOutput(run_100)
wav_100,flux_100 = m_100.get_sed(inclination=0,aperture=-1)

flux_100 = np.asarray(flux_100)*u.erg/u.s

flux_100 /= (4.*3.14*dl**2.)

flux_100 /= nu
flux_100 = flux_100.to(u.mJy)


flux_diff = flux_1mil - flux_100

frac_diff = flux_diff/flux_1mil

filtered_diff = frac_diff[~np.isnan(frac_diff)]

with open('halo_1_differences.txt', 'a') as file:
    np.savetxt(file, abs(filtered_diff.value))

#print(frac_diff)

#print(wav)

print("norms over whole redshifted wav range")

L_inf = np.max(abs(filtered_diff))
L_1 = abs(filtered_diff).sum()/len(filtered_diff)
L_2 = (((abs(filtered_diff)**2).sum())/len(filtered_diff))**(1/2)

print("inf norm = ", L_inf)
print("L_1 = ", L_1)
print("L_2 = ", L_2)


print("\nnorms for rest-frame UV region")

wav_index = np.array(wav/(1.+z))
UV_end_index = (np.abs(wav_index - 0.01)).argmin()
UV_start_index = (np.abs(wav_index - 0.4)).argmin()

L_inf_UV = np.max(abs(filtered_diff[UV_start_index:UV_end_index + 1]))
L_1_UV = abs(filtered_diff[UV_start_index:UV_end_index + 1]).sum()/len(filtered_diff[UV_start_index:UV_end_index + 1])
L_2_UV = (((abs(filtered_diff[UV_start_index:UV_end_index + 1])**2).sum())/len(filtered_diff[UV_start_index:UV_end_index + 1]))**(1/2)

print("inf norm = ", L_inf_UV)
print("L_1 = ", L_1_UV)
print("L_2 = ", L_2_UV)

print("\nnorms for rest-frame visibile region")

vis_end_index = (np.abs(wav_index - 0.4)).argmin()
vis_start_index = (np.abs(wav_index - 0.7)).argmin()

L_inf_vis = np.max(abs(filtered_diff[vis_start_index:vis_end_index + 1]))
L_1_vis = abs(filtered_diff[vis_start_index:vis_end_index + 1]).sum()/len(filtered_diff[vis_start_index:vis_end_index + 1])
L_2_vis = (((abs(filtered_diff[vis_start_index:vis_end_index + 1])**2).sum())/len(filtered_diff[vis_start_index:vis_end_index + 1]))**(1/2)

print("inf norm = ", L_inf_vis)
print("L_1 = ", L_1_vis)
print("L_2 = ", L_2_vis)

#for i in range(frac_diff.shape[0]):
#    flux_diff = flux_1mil[i,:] - flux_100[i,:]
#    ax.loglog(wav,frac_diff[i,:])

points = np.array([wav, abs(frac_diff)]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
colors = ['red' if value < 0 else 'blue' for value in frac_diff]

lc = LineCollection(segments, colors=colors, linewidth=2)
ax.add_collection(lc)
ax.set_yscale('log')
ax.set_xscale('log')

#ax.loglog(wav,frac_diff)

#for i in range(frac_diff.size):
#    ax.loglog(wav[i],frac_diff[i])
#    print("wave", np.log(wav[i]/u.micron))
#    print("frac_diff", np.log(frac_diff[i]))

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel('Fractional difference')
#ax.set_ylim([1,1e8])
#ax.set_xlim(0.05,15000)
ax.grid()

fig.savefig('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_3000000/pd_sed/sed_compare_frac_diff_halo_1.png')


