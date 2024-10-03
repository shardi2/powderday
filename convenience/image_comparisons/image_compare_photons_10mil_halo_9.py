import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from hyperion.model import ModelOutput
import astropy.units as u
import h5py

import sys
np.set_printoptions(threshold=sys.maxsize)
from matplotlib.colors import ListedColormap

# ------------------------
# modifiable header
# ------------------------

m_1mil = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_1000000/halo_63/example.0063.rtout.image')
m_10mil = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000000_1angle/halo_63/example.0063.rtout.image')
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/pd_image/example.0000.rtout.image')
wav = 0.55816454  # micron
#f = h5py.File("final_dataset_halo_catalog.hdf5", "a")
# ------------------------

# Get the image from the ModelOutput object
image_1mil = m_1mil.get_image(units='ergs/s')
image_10mil = m_10mil.get_image(units='ergs/s')

# Open figure and create axes
fig = plt.figure()
ax = fig.add_subplot(111)

# Find the closest wavelength
iwav = np.argmin(np.abs(wav - image_1mil.wav))

# Calculate the image width in kpc
w = image_1mil.x_max * u.cm
w = w.to(u.kpc)

image_diff = (image_10mil.val[0, :, :, iwav]) - (image_1mil.val[0, :, :, iwav])

frac_diff = image_diff / image_10mil.val[0, :, :, iwav].mean()

#pos_mask = frac_diff > 0
#neg_mask = frac_diff < 0

#pos_frac_diff = np.where(pos_mask, frac_diff, np.nan)
#neg_frac_diff = np.where(neg_mask, frac_diff, np.nan)

#print(frac_diff)

minv = np.min(image_diff)

maxv = np.max(image_diff)

filtered_diff = frac_diff[~np.isinf(frac_diff)]

filtered_diff_2 = filtered_diff[~np.isnan(filtered_diff)]

with open('halo_63_differences.txt', 'a') as file:
    np.savetxt(file, abs(filtered_diff_2))

#print(filtered_diff_2)

L_inf = np.max(abs(filtered_diff_2))
L_1 = abs(filtered_diff_2).sum()/len(filtered_diff_2)
L_2 = (((abs(filtered_diff_2)**2).sum())/len(filtered_diff_2))**(1/2)

print("inf norm = ", L_inf)
print("L_1 = ", L_1)
print("L_2 = ", L_2)

print("difference in total luminosities:", (image_10mil.val[0, :, :, iwav]).sum() - (image_1mil.val[0, :, :, iwav]).sum())

#my_cmap = plt.cm.PuOr(np.arange(plt.cm.PuOr.N))
#my_cmap[:,3] = 0.5
#my_cmap = ListedColormap(my_cmap)

# plot the beast
#cax = ax.imshow(final_image_log, cmap=plt.cm.viridis, vmin=minv, vmax=maxv,
#                origin='lower', extent=[-w.value, w.value, -w.value, w.value])
cax = ax.imshow(frac_diff, cmap=plt.cm.PuOr, norm=mpl.colors.SymLogNorm(linthresh=0.01, vmin=-10, vmax=10), origin='lower', extent=[-w.value, w.value, -w.value, w.value])
#, norm=colors.LogNorm())
#cax = ax.imshow(pos_frac_diff, cmap=plt.cm.PuOr, norm=mpl.colors.SymLogNorm(linthresh=0.001, vmin=0.1, vmax=10), origin='lower', extent=[-w.value, w.value, -w.value, w.value])
#cax = ax.imshow(abs(neg_frac_diff), cmap=my_cmap, norm=mpl.colors.SymLogNorm(linthresh=0.001, vmin=0.1, vmax=10), origin='lower', extent=[-w.value, w.value, -w.value, w.value])

# Finalize the plot
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel('x (kpc)')
ax.set_xlabel('y (kpc)')

plt.colorbar(cax, label='Difference Between Runs')

fig.savefig('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000000_1angle/pd_image/pd_image_compare_diff_halo_63.png', bbox_inches='tight', dpi=150)

