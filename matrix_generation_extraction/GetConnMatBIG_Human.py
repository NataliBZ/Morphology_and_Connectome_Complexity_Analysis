# Natali Barros Zulaica
# 12.03.2024

#This python script loads the connectome for the human STRUCTURAL microcircuit. 
#This means that synapses are placed at ALL touch locations.

#It is based on the following notebook 
# https://bbpgitlab.epfl.ch/conn/papers/sscx_and_em/-/blob/main/notebooks/microns/Schneider-Mizell/100%20um%20microns%20struc.ipynb?ref_type=heads

#We will load connectomes of a square superfice of 960um x 960um
#so we could obtain 20000 neurons (nodes) approx.


import numpy, pandas, conntility
import numpy as np
import bluepysnap as snap
import voxcell
import json
from matplotlib import pyplot as plt
import h5py
import tqdm
import LoadConnMat as lcm # package with useful functions to load connectivity matrices for this case.

## 1. load the circuit and the positions and orientations of cells in the respective flatmaps

fn_h = "/gpfs/bbp.cscs.ch/project/proj159/circuit_building/270923/run/sonata/circuit_config.json"
h_atlas_dir = "/gpfs/bbp.cscs.ch/project/proj159/circuit_building/270923/atlas/"

fm_h = voxcell.VoxelData.load_nrrd(h_atlas_dir + "flatmap.nrrd")
orient_h = voxcell.VoxelData.load_nrrd(h_atlas_dir + "orientation.nrrd")
circ_h = snap.Circuit(fn_h)

syn = 'EXC'
layers_h = ['L2','L3a','L3b']
range_x_h = [700, 1660] # square space similar to human circuit radius for 20000 nodes
range_y_h = [700, 1660]
pixel_sz_h = 16.0

nrn_h = lcm.load_conn_node(circ_h, layers_h, syn, range_x_h, range_y_h, fm_h, orient_h, pixel_sz_h)

## 2. load the connectivity

con_fn_h = "/gpfs/bbp.cscs.ch/project/proj159/circuit_building/270923/run/sonata/networks/edges/structural/hncx_neurons__hncx_neurons_chemical_synapse/edges.h5"
h5_h = h5py.File(con_fn_h, "r")
ranges_h = h5_h["edges/hncx_neurons__hncx_neurons_chemical_synapse/indices/target_to_source/node_id_to_ranges"][nrn_h["node_ids"], :]
edge_ids_h = h5_h["edges/hncx_neurons__hncx_neurons_chemical_synapse/indices/target_to_source/range_to_edge_id"][ranges_h[:, 0]]

E_h = lcm.load_conn_edges(h5_h, edge_ids_h, 'human')

## 3. save flat locations into the dataframe E_rat

xyz = ["afferent_center_x", "afferent_center_y", "afferent_center_z"]

E_h = lcm.save_flat_loc(E_h, xyz, fm_h, orient_h, pixel_size=pixel_sz_h)

## 4. Build connectivity object, compress it and save it

M_h = lcm.build_conn_mat(E_h, nrn_h)

# Safe
fname_h_raw = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Human_960um_squared_struc_conmat.h5"
M_h.to_h5(fname_h_raw)

# Compress and safe
fname_h_raw = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Human_960um_squared_struc_conmat_compressed.h5"
M_h_comp = M_h.compress()
M_h_comp.to_h5(fname_h_raw)

# As many of the synapses will be out of the selected range we need to filter them 
# The ConnectivityMatrix object expects the indices of pre- and post-synaptic neurons to be 
# in the range of 0, 1, 2, ..., N, where N is the size of the population. 
# However, in this case we were giving it the gids of the neurons which can be all over the 
# place (0, 5, 20, 25, 100, ....), so we need to re-index

# filter, re-index, compress for faster loading and safe

fname_h_filt = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Human_960um_squared_struc_conmat_filtered.h5"
M_h = lcm.filter_and_reindex(M_h)
M_h.to_h5(fname_h_filt)

fname_h_comp = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Human_960um_squared_struc_conmat_filtered_compressed.h5"
M_h_comp = M_h.compress()
M_h_comp.to_h5(fname_h_comp)