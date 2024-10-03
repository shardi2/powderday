import yt
import os, sys
import time
from yt.config import ytcfg
import ytree
import numpy as np
import glob

yt.enable_parallelism()
yt.set_log_level(30)

nproc = ytcfg.get("yt", "internals", "topcomm_parallel_size")

zcat = dict(z=[], fn=[])
with open('/storage/home/hcoda1/0/jw254/data/SG64-2020/redshifts.dat', 'r') as fp:
    for line in fp:
        zcat['z'].append(float(line.split("=")[-1]))
        zcat['fn'].append(line.split(":")[0])
    zcat['z'] = np.array(zcat['z'])

def get_output_list(search_pattern):
    files = glob.glob(search_pattern)
    files.sort()
    allds = {}
    for f in files:
        allds[f] = yt.load(f)
    return allds

def get_filename_from_redshift(z):
    i = np.argmin(np.abs(z - zcat['z']))
    return zcat['fn'][i]

if os.path.exists("arbor/arbor.h5"):
    a = ytree.load("arbor/arbor.h5")
else:
    a = ytree.load("/storage/home/hcoda1/0/jw254/data/SG64-2020/rockstar_halos-jhw/trees/tree_0_0_0.dat")

if "stellar_mass" not in a.field_list:
    a.add_analysis_field("stellar_mass", default=-1, units="Msun")

if "gas_mass" not in a.field_list:
    a.add_analysis_field("gas_mass", default=-1, units="Msun")

if "total_pop3_mass" not in a.field_list:
    a.add_analysis_field("total_pop3_mass", default=-1, units="Msun")

if "total_cold_gas_mass" not in a.field_list:
    a.add_analysis_field("total_cold_gas_mass", default=-1, units="Msun")

if "total_HI_mass" not in a.field_list:
    a.add_analysis_field("total_HI_mass", default=-1, units="Msun")

if "total_HII_mass" not in a.field_list:
    a.add_analysis_field("total_HII_mass", default=-1, units="Msun")

if "stellar_mass_10Myr" not in a.field_list:
    a.add_analysis_field("stellar_mass_10Myr", default=-1, units="Msun")

if "stellar_mass_100Myr" not in a.field_list:
    a.add_analysis_field("stellar_mass_100Myr", default=-1, units="Msun")

if "total_stellar_angular_momentum" not in a.field_list:
    a.add_analysis_field("total_stellar_angular_momentum", default=-1, units='cm**2*g/s')

def p2(pfilter, data):
    return (data[('nbody', 'particle_type')] == 7)

yt.add_particle_filter(function=p2, name='p2', requires=['particle_mass', 'particle_type'])

def p3(pfilter, data):
    return ((data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') < 1e-10)) \
        |  ((data[('nbody', 'particle_type')] == 1) & (data[('nbody', 'creation_time')] > 0) & \
            (data[('nbody', 'particle_mass')].in_units('Msun') > 1)) \
        |  ((data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') > 1e-3))

yt.add_particle_filter(function=p3, name='p3', requires=['particle_mass', 'particle_type', 'creation_time'])

def p3_sn(pfilter, data):
    return (data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') < 1e-10)

yt.add_particle_filter(function=p3_sn, name='p3_sn', requires=['particle_mass', 'particle_type'])

def p3_bh(pfilter, data):
    return (data[('nbody', 'particle_type')] == 1) & (data[('nbody', 'creation_time')] > 0) & \
        (data[('nbody', 'particle_mass')].in_units('Msun') > 1e-3)

yt.add_particle_filter(function=p3_bh, name='p3_bh', requires=['creation_time', 'particle_mass', 'particle_type'])

def p3_living(pfilter, data):
    return (data[('nbody', 'particle_type')] == 5) & (data[('nbody', 'particle_mass')].in_units('Msun') > 1e-3)

yt.add_particle_filter(function=p3_living, name='p3_living', requires=['particle_mass', 'particle_type'])

def calculate_stellar_mass(node):
    sphere = node.sphere
    node["stellar_mass"] = sphere.quantities.total_quantity(("p2", "particle_mass"))

def calculate_gas_mass(node):
    sphere = node.sphere
    node["gas_mass"] = sphere.quantities.total_quantity(("gas", "mass"))

def calculate_total_pop3_mass(node):
    sphere = node.sphere
    p3_liv = sphere.quantities.total_quantity(("p3_living", "particle_mass"))
    p3_bh = sphere.quantities.total_quantity(("p3_bh", "particle_mass"))
    p3_sn = sphere.quantities.total_quantity(("p3_sn", "particle_mass"))
    node["total_pop3_mass"] = p3_liv + p3_bh + (p3_sn * (1e20))

def calculate_total_cold_gas_mass(node):
    ds = node.ds
    sphere = node.sphere
    cr = ds.cut_region(sphere, ["obj[('gas', 'temperature')] < 1e3"])
    node["total_cold_gas_mass"] = cr.quantities.total_quantity(('gas', 'mass'))
#    node["total_cold_gas_mass"] = 0
#    for i, temp in enumerate(sphere["gas", "temperature"]):
#        if temp < 1000:
#            node["total_cold_gas_mass"] += sphere["gas", "mass"][i]

def calculate_total_HI_mass(node):
    sphere = node.sphere
    node['total_HI_mass'] = sphere.quantities.total_quantity(("gas", "H_p0_mass"))
    
def calculate_total_HII_mass(node):
    sphere = node.sphere
    node['total_HII_mass'] = sphere.quantities.total_quantity(("gas", "H_p1_mass"))

def calculate_stellar_mass_10Myr(node):
    ds = node.ds
    sphere = node.sphere
    cr = ds.cut_region(sphere, ["obj[('p2', 'age')].to('Myr') < 10"])
    node["stellar_mass_10Myr"] = cr.quantities.total_quantity(('p2', 'particle_mass'))

def calculate_stellar_mass_100Myr(node):
    ds = node.ds
    sphere = node.sphere
    cr = ds.cut_region(sphere, ["obj[('p2', 'age')].to('Myr') < 100"])
    node["stellar_mass_100Myr"] = cr.quantities.total_quantity(('p2', 'particle_mass'))

def calculate_angular_momentum(node):
    sphere = node.sphere
    node["total_stellar_angular_momentum"] = sphere.quantities.total_quantity(("p2", "particle_angular_momentum")) 

def stellar_mass_recipe(pipeline, outputs):
    #pipeline.add_operation(p2)
    pipeline.add_operation(get_yt_dataset, outputs)
    pipeline.add_operation(get_yt_sphere)
    pipeline.add_operation(calculate_stellar_mass)
    pipeline.add_operation(calculate_gas_mass)
    pipeline.add_operation(calculate_total_pop3_mass)
    pipeline.add_operation(calculate_total_cold_gas_mass)
    pipeline.add_operation(calculate_total_HI_mass)
    pipeline.add_operation(calculate_total_HII_mass)
    pipeline.add_operation(calculate_stellar_mass_10Myr)
    pipeline.add_operation(calculate_stellar_mass_100Myr)
    pipeline.add_operation(calculate_angular_momentum)
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

def get_yt_dataset(node, outputs):
    # assume you have something like this
    filename = get_filename_from_redshift(node["redshift"])
    # attach it to the node for later use
    node.ds = yt.load("/storage/home/hcoda1/0/jw254/data/SG64-2020/"+str(filename))
    #yt.add_particle_filter(function=p2, name='p2', requires=['particle_mass', 'particle_type'])
    if "p2" not in node.ds.particle_types:
        node.ds.add_particle_filter("p2")
    if "p3" not in node.ds.particle_types:
        node.ds.add_particle_filter("p3")
    if "p3_sn" not in node.ds.particle_types:
        node.ds.add_particle_filter("p3_sn")
    if "p3_bh" not in node.ds.particle_types:
        node.ds.add_particle_filter("p3_bh")
    if "p3_living" not in node.ds.particle_types:
        node.ds.add_particle_filter("p3_living")


def get_yt_sphere(node):
    # this works if get_yt_dataset has been called first
    ds = node.ds

    center = node["position"].to("unitary").v
    radius = node["virial_radius"].to("unitary").to_value()
    node.sphere = ds.sphere((center), (radius, "unitary"))

allds = get_output_list("DD*/output_????")
ap = ytree.AnalysisPipeline()
ap.add_recipe(stellar_mass_recipe, allds)
# later, in the pipeline
ap.add_operation(delete_dataset, always_do=True)

trees = list(a[:])
dynamic = nproc > 2
for tree in ytree.parallel_trees(trees, dynamic=dynamic):
#for node in ytree.parallel_trees(trees):   
    for node in tree["forest"]:
        ap.process_target(node)

#if yt.is_root():
#    a.save_arbor(trees=trees)
