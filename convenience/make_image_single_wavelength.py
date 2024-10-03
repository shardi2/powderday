import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.units as u
import h5py
import sys
import numpy
numpy.set_printoptions(threshold=sys.maxsize)

# ------------------------
# modifiable header
# ------------------------

m = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000000/halo_0/example.0000.rtout.image')
#'/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000/pd_image/example.0000.rtout.image')
wav = 0.55816454  # micron
#f = h5py.File("final_dataset_halo_catalog.hdf5", "a")
# ------------------------

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

# plot the beast
#cax = ax.imshow(final_image_log, cmap=plt.cm.viridis, vmin=minv, vmax=maxv,
#                origin='lower', extent=[-w.value, w.value, -w.value, w.value])
cax = ax.imshow(np.log10(image.val[0, :, :, iwav]), cmap=plt.cm.viridis,
                origin='lower', extent=[-w.value, w.value, -w.value, w.value])
#, norm=colors.LogNorm())

#print(image.val[0, :, :, iwav])

print(np.log10(image.val[0, :, :, iwav]))

# Finalize the plot
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel('x (kpc)')
ax.set_xlabel('y (kpc)')

plt.colorbar(cax, label='log10 Luminosity (ergs/s)')

fig.savefig('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_10000000/pd_image/pd_image_100000.png', bbox_inches='tight', dpi=150)
