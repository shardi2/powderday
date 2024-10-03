import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.units as u
import h5py

# ------------------------
# modifiable header
# ------------------------

m = ModelOutput('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_100000/example.0000.rtout.image')
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


#maxv = image.max()
#print(maxv)
#minv = maxv - 15
final_image = image.val[0, :, :, iwav]
final_image[final_image == 0] = 1e-300
final_image_log = np.log10(final_image)
maxv = final_image_log.max()
#minv = np.min(final_image_log[final_image_log>1])
minv = maxv - 15
# plot the beast
cax = ax.imshow(final_image_log, cmap=plt.cm.viridis, vmin=minv, vmax=maxv,
                origin='lower', extent=[-w.value, w.value, -w.value, w.value])


# Finalize the plot
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel('x (kpc)')
ax.set_xlabel('y (kpc)')

plt.colorbar(cax, label='log Luminosity (ergs/s)')

fig.savefig('pd_image_100000.png', bbox_inches='tight', dpi=150)


f.create_dataset('halo_3/1000000/image', data=final_image_log)