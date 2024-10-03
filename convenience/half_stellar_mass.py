import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.units as u
import h5py
import ytree
import yt
yt.enable_plugins()


# ------------------------
# modifiable header
# ------------------------

a = ytree.load("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/arbor/arbor.h5")
m = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/halo_0/example.0000.rtout.image')
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/pd_image/example.0000.rtout.image')
wav = 0.55816454  # micron
f = h5py.File("/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/Powderday/powderday/final_dataset_halo_catalog.hdf5", "a")
#f = h5py.File("final_dataset_halo_catalog.hdf5", "a")
# ------------------------

#getting star filter
ds = yt.load('/storage/home/hcoda1/0/jw254/data/SG64-2020/DD0125/output_0125')

def p2(pfilter, data):
    return (data[('nbody', 'particle_type')] == 7)
yt.add_particle_filter(function=p2, name='p2', requires=['particle_mass', 'particle_type'])

ds.add_particle_filter("p2")
ds.field_list

#getting sphere around halo
radius = ds.quan(a['virial_radius'][0], 'kpccm')
radius_kpc = radius.to('kpc')
center = a['position'][0]
ad = ds.sphere(center, radius)

#positions of stars in sphere
star_radii = ad["p2", "particle_radius"]
star_radii_kpc = star_radii.to('kpc')

#mass of stars in spheres
star_masses = ad["p2", "particle_mass"]

#getting indexes of sorted positions of stars
ir = np.argsort(star_radii_kpc)

#sorting positions of stars
r = star_radii_kpc[ir]

#sorting masses of stars
m = star_masses[ir]

cumulative_sum_m = np.cumsum(m)

#find half mass
half_mass = cumulative_sum_m[-1] * 0.5

#find radius of half mass
half_mass_rad = np.interp(half_mass, cumulative_sum_m, r)
print('half stellar mass radius', half_mass_rad) 

# Open figure and create axes
fig = plt.figure()
ax = fig.add_subplot(111)

#marking half stellar mass radius
ax.axvline(x=half_mass_rad, color='r', linestyle='--', label='Half Stellar Mass Radius')

#plot cumulative sum
ax.plot(r, cumulative_sum_m, label='Cumulative Sum')

# Adding labels and legend
ax.set_xlabel('Radius')
ax.set_ylabel('Cumulative Mass')
ax.legend()

fig.savefig('half_mass_check.png')

if 'halo_0/photons_1e4/half_stellar_mass_rad' in f:
    del f['halo_0/photons_1e4/half_stellar_mass_rad']
    f.attrs['halo_0/photons_1e4/half_stellar_mass_rad'] = half_mass_rad
else:
    f.attrs['halo_0/photons_1e4/half_stellar_mass_rad'] = half_mass_rad

f.close()

