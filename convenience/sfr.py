import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.units as u
import h5py
import ytree
import yt
import unyt
yt.enable_plugins()


# ------------------------
# modifiable header
# ------------------------

a = ytree.load("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/arbor/arbor.h5")
#m = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/halo_0/example.0000.rtout.image')
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/pd_image/example.0000.rtout.image')
wav = 0.55816454  # micron
f = h5py.File("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/final_dataset_halo_catalog.hdf5", "a")
#f = h5py.File("final_dataset_halo_catalog.hdf5", "a")
# ------------------------

#initialize variables
stellar_mass_10 = 0
stellar_mass_100 = 0
total_stellar_mass = 0

#total stellar mass of all halos 
for i in range(0, a.size, 1):
    stellar_mass_10 += a["stellar_mass_10Myr"][i]
    stellar_mass_100 += a["stellar_mass_100Myr"][i]
    total_stellar_mass += a["stellar_mass"][i]

SFR_10 = stellar_mass_10 / (10 * unyt.Myr)

sSFR_10 = SFR_10 / (total_stellar_mass)

SFR_100 = stellar_mass_100 / (100 * unyt.Myr)

sSFR_100 = SFR_100 / (total_stellar_mass)

print(SFR_10)

print(SFR_100)

if 'halo_0/SFR_10Myr' in f:
    del f['halo_0/SFR_10Myr']
    f['halo_0'].attrs['SFR_10Myr'] = SFR_10
else:
    f['halo_0'].attrs['SFR_10Myr'] = SFR_10

if 'halo_0/sSFRi_10Myr' in f:
    del f['halo_0/sSFR_10Myr']
    f['halo_0'].attrs['sSFR_10Myr'] = sSFR_10
else:
    f['halo_0'].attrs['sSFR_10Myr'] = sSFR_10

if 'halo_0/SFR_100Myr' in f:
    del f['halo_0/SFR_100Myr']
    f['halo_0'].attrs['SFR_100Myr'] = SFR_100
else:
    f['halo_0'].attrs['SFR_100Myr'] = SFR_100

if 'halo_0/sSFRi_100Myr' in f:
    del f['halo_0/sSFR_100Myr']
    f['halo_0'].attrs['sSFR_100Myr'] = sSFR_100
else:
    f['halo_0'].attrs['sSFR_100Myr'] = sSFR_100

f.close()

