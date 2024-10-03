import matplotlib.pyplot as plt
import numpy as np
import h5py
import math

f = h5py.File('/storage/home/hcoda1/7/shardin31/p-jw254-0/Research/summer2023/pd_test/run_1000/convolved.0125.hdf5', 'r')


convolved_image = f['image_data'][0]
filter_name = f['filter_names'][0].astype(str)

fig = plt.figure()
ax = fig.add_subplot(111)
#fig, axes = plt.subplots(len(convolved_image), 1, figsize=(6, 6 * len(convolved_image)))

w = f['image_data'].attrs['width']
w_unit = f['image_data'].attrs['width_unit'].astype(str)

final_convolved = []
x = 0
for i in convolved_image:
    #i = list(filter(lambda a: a != 0, i))
    #i = [float(math.log(x)) for x in i]
    #print(type(i))
    i[i == 0] = 1e-300
    i = np.array(i, dtype=float)
    #print(i)
    #print(len(i))
    #x += 1
    #print(x)
    #cax = ax.imshow(i.reshape(1,-1), cmap=plt.cm.viridis, origin='lower', extent=[-w, w, -w, w])
    final_convolved.append(i)		
		
#print(final_convolved)

#final_convolved = np.array(final_convolved_log)
#print(convolved_image[0])

#final_convolved_log = [np.log(item) for item in final_convolved]
final_convolved = np.array(final_convolved, dtype=float)

#print(convolved_image.size)
#final_convolved_masked = np.ma.array(final_convolved_log, dtype=float)

#print(type(final_convolved))

cax = ax.imshow(np.log(final_convolved), cmap=plt.cm.viridis,
                origin='lower', extent=[-w, w, -w, w])
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel('x ({})'.format(w_unit))
ax.set_ylabel('y ({})'.format(w_unit))
plt.colorbar(cax, label='log Luminosity (ergs/s)', format='%.0e')
plt.title("Convolved image: {}".format(filter_name))
plt.tight_layout()
plt.show()
fig.savefig("new_image.png")
#plt.show()
