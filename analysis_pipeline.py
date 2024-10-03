import yt
import os, sys
import time
from yt.config import ytcfg
import ytree
import numpy as np
yt.enable_parallelism()

nproc = ytcfg.get("yt", "internals", "topcomm_parallel_size")

zcat = dict(z=[], fn=[])
with open('/storage/home/hcoda1/0/jw254/data/SG64-2020/redshifts.dat', 'r') as fp:
    for line in fp:
        zcat['z'].append(float(line.split("=")[-1]))
        zcat['fn'].append(line.split(":")[0])
    zcat['z'] = np.array(zcat['z'])

def get_filename_from_redshift(z):
    i = np.argmin(np.abs(z - zcat['z']))
    return zcat['fn'][i]

if os.path.exists("arbor/arbor.h5"):
    a = ytree.load("arbor/arbor.h5")
else:
    a = ytree.load("/storage/home/hcoda1/0/jw254/data/SG64-2020/rockstar_halos-jhw/trees/tree_0_0_0.dat")

if "total_stellar_mass" not in a.field_list:
    a.add_analysis_field("total_stellar_mass", default=-1, units="Msun")

#max_array_size = 205
#if "star_masses" not in a.field_list:
#    a.add_analysis_field("star_masses", default=np.full(max_array_size, -1), units="g", dtype = np.ndarray)

#if "star_positions" not in a.field_list:
#    a.add_analysis_field("star_positions", default=-1, units="unitary", dtype = np.ndarray)

#if "star_ages" not in a.field_list:
#    a.add_analysis_field("star_ages", default=-1, units="yr", dtype = np.ndarray)

#if "star_radii" not in a.field_list:
#    a.add_analysis_field("star_radii", default=-1, units="cm", dtype = np.ndarray)

def p2(pfilter, data):
    return (data[('nbody', 'particle_type')] == 7)

yt.add_particle_filter(function=p2, name='p2', requires=['particle_mass', 'particle_type'])

def calculate_stellar_mass(node):
    sphere = node.sphere
    node["total_stellar_mass"] = sphere.quantities.total_quantity(("p2", "particle_mass"))

def save_star_masses(node):
    print("MAX_SIZE:", max_array_size)
    sphere = node.sphere
    #node_array = np.array(sphere["p2","particle_mass"])
    #size_difference = max_array_size - len(node_array)
    #size_difference = max_array_size - len(sphere["p2","particle_mass"])
    #print("SIZE DIFFERENCE", size_difference)
    print("ARRAY LENGTH", len(sphere["p2","particle_mass"]))
    print("ARRAY", sphere["p2","particle_mass"])
    #print("NODE", node["star_masses"])
    nstars = sphere["p2","particle_mass"].size
    node["star_masses"] = np.zeros(max_array_size)
    print("NODE", node["star_masses"].size)
    #node["star_masses"][:nstars] = sphere["p2","particle_mass"]
    #node["star_masses"] = np.concatenate((np.array(sphere["p2","particle_mass"]),np.zeros(size_difference)), dtype=float) 

def save_star_positions(node):
    sphere = node.sphere
    node["star_positions"] = np.array(sphere["p2","particle_position"])

def save_star_ages(node):
    sphere = node.sphere
    node["star_ages"] = np.array(sphere["p2","age"])

def save_star_radii(node):
    sphere = node.sphere
    node["star_radii"] = np.array(sphere["p2","particle_radius"])

def max_stellar_mass_size(node, maximum_array_size):
    max_size = maximum_array_size
    #for tree in ytree.parallel_trees(trees, dynamic=dynamic):
    #    for node in tree["forest"]:
    get_yt_dataset(node)
    get_yt_sphere(node)
    sphere = node.sphere
    if max_size < np.array(sphere["p2","particle_mass"]).size:
        max_size = np.array(sphere["p2","particle_mass"]).size
    return max_size

def stellar_mass_recipe(pipeline):
    pipeline.add_operation(get_yt_dataset)
    pipeline.add_operation(get_yt_sphere)
    pipeline.add_operation(calculate_stellar_mass)
#    pipeline.add_operation(max_stellar_mass_size) 
    pipeline.add_operation(save_star_masses)
#    pipeline.add_operation(save_star_positions)
#    pipeline.add_operation(save_star_ages)
#    pipeline.add_operation(save_star_radii)
    #pipeline.add_operation(delete_stellar_mass_field)
    pipeline.add_operation(delete_attributes, ["ds", "sphere"])

def delete_attributes(node, attributes):
    for attr in attributes:
        if hasattr(node, attr):
            delattr(node, attr)

def delete_dataset(node):
    if hasattr(node, "ds"):
        node.ds.index.clear_all_data()
        del node.ds.index.grid_dimensions
        del node.ds.index.grid_left_edge
        del node.ds.index.grid_right_edge
        del node.ds.index.grid_levels
        del node.ds.index.grid_particle_count
        del node.ds.index.grids

def delete_stellar_mass_field(node):
    if "stellar_mass" in node.field_list:
        #node["total_stellar_mass"] = node["stellar_mass"]
        node.remove_field("stellar_mass")

def get_yt_dataset(node):
    # assume you have something like this
    filename = get_filename_from_redshift(node["redshift"])
    # attach it to the node for later use
    node.ds = yt.load("/storage/home/hcoda1/0/jw254/data/SG64-2020/"+str(filename))
    #yt.add_particle_filter(function=p2, name='p2', requires=['particle_mass', 'particle_type'])
    node.ds.add_particle_filter("p2")


def get_yt_sphere(node):
    # this works if get_yt_dataset has been called first
    ds = node.ds

    center = node["position"].to("unitary").v
    radius = node["virial_radius"].to("unitary").to_value()
    node.sphere = ds.sphere((center), (radius, "unitary"))

ap = ytree.AnalysisPipeline()
#ap.add_operation(max_stellar_mass_size, always_do=True)
ap.add_recipe(stellar_mass_recipe)
# later, in the pipeline
ap.add_operation(delete_dataset, always_do=True)

#global max_array_size
max_array_size = 205
trees = list(a[:])
dynamic = nproc > 2

#for tree in ytree.parallel_trees(trees, dynamic=dynamic):
#    for node in tree["forest"]:
#        max_array_size = max_stellar_mass_size(node, maximum_array_size = max_array_size)

if "star_masses" not in a.field_list:
    a.add_analysis_field("star_masses", default=np.zeros(max_array_size), units="g", dtype = np.ndarray)

for tree in ytree.parallel_trees(trees, dynamic=dynamic):
#for node in ytree.parallel_trees(trees):   
    for node in tree["forest"]:
#        max_array_size = max_stellar_mass_size(node, maximum_array_size = max_array_size)
        ap.process_target(node)
        break
    break

#if yt.is_root():
#    a.save_arbor(trees=trees)
