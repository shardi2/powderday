import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.units as u
import h5py
import ytree

# ------------------------
# modifiable header
# ------------------------

a = ytree.load("/storage/home/hcoda1/0/jw254/data/SG64-2020/rockstar_halos-jhw/trees/tree_0_0_0.dat")
m = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000_20/example.0000.rtout.image')
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/pd_image/example.0000.rtout.image')
wav = 0.55816454  # micron
#f = h5py.File("final_dataset_halo_catalog.hdf5", "a")
# ------------------------


#radial profile = averaged brightness of all of the pixels at a certain radius
def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    #print("tbin", tbin)
    #nr = np.bincount(r.ravel())
    #print("nr", nr)
    #radialprofile = tbin / nr
    radialprofile = tbin
    #print(radialprofile)
    unique_r = np.unique(r)
    return radialprofile, unique_r

# Get the image from the ModelOutput object
image = m.get_image(units='ergs/s')

# Open figure and create axes
fig = plt.figure()
ax = fig.add_subplot(111)

# Find the closest wavelength
iwav = np.argmin(np.abs(wav - image.wav))

# Calculate the image width in kpc
w = image.x_max * u.cm
w = w.to(u.kpc)

#find radial profile
virial_radius = a["virial_radius"][0].to("kpc")*a["scale"][0]
nx, ny = image.val[0, :, :, iwav].shape
center = (nx / 2, ny /2)
rad_profile, r = radial_profile(image.val[0, :, :, iwav], center)

print("rad profile:", rad_profile)
print("radius:", r)
#gets sum of all x values (brightness) up to each x value; adds up shell brightnesses up to each radius
cumulative_sum = np.cumsum(rad_profile)
print("cumulative sum:", cumulative_sum)

half_light = cumulative_sum[-1] * 0.5

#find half light radius
#half_light_rad = np.interp(half_light, cumulative_sum, np.arange(len(cumulative_sum)))
half_light_rad = np.interp(half_light, cumulative_sum, r)
print("half light", half_light_rad)

#convert to kpc
kpc_half_light_rad = half_light_rad * ((virial_radius * 0.20) / len(rad_profile))

# Marking Half-Light Radius
ax.axvline(x=half_light_rad, color='r', linestyle='--', label='Half-Light Radius')

#plot half light radius
ax.plot(r, rad_profile, label='Radial Profile')

#plot cumulative sum
ax.plot(r, cumulative_sum, label='Cumulative Sum')

# Adding labels and legend
ax.set_xlabel('Radius')
ax.set_ylabel('Brightness')
ax.legend()

fig.savefig('half_light_check_20.png')

print("half light kpc:", kpc_half_light_rad)
