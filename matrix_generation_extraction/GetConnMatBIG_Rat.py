# Natali Barros Zulaica
# 12.03.2024

#This python script loads the connectome for the rat STRUCTURAL microcircuit. 
#This means that synapses are placed at ALL touch locations.

#It is based on the following notebook 
# https://bbpgitlab.epfl.ch/conn/papers/sscx_and_em/-/blob/main/notebooks/microns/Schneider-Mizell/100%20um%20microns%20struc.ipynb?ref_type=heads

#We will load connectomes of a square superfice of 623um x 623um
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

fn_rat = "/gpfs/bbp.cscs.ch/project/proj83/jira-tickets/NSETM-1948-extract-hex-O1/data/S1_data/circuit_config.json"
rat_atlas_dir = "/gpfs/bbp.cscs.ch/project/proj83/jira-tickets/NSETM-1948-extract-hex-O1/data/O1_data_physiology/atlas/"

fm_rat = voxcell.VoxelData.load_nrrd(rat_atlas_dir + "flatmap.nrrd")
orient_rat = voxcell.VoxelData.load_nrrd(rat_atlas_dir + "orientation.nrrd")
circ_rat = snap.Circuit(fn_rat)

syn = 'EXC'
layers_rat = ['2','3']
range_x_rat = [1000,1623] #square to match same num nodes (cells) as in human (aprox 20000)
range_y_rat = [1000,1623]
pixel_sz_rat = 34.0
node_pop_rat = 'S1nonbarrel_neurons'

nrn_rat = lcm.load_conn_node(circ_rat, layers_rat, syn, range_x_rat, range_y_rat, fm_rat, orient_rat, pixel_sz_rat, node_population_name=node_pop_rat)

## 2. load the connectivity

con_fn_rat = "/gpfs/bbp.cscs.ch/project/proj83/circuits/Bio_M/20200805/connectome/structural/edges.sonata"
h5_rat = h5py.File(con_fn_rat, "r")
ranges_rat = h5_rat["edges/default/indices/target_to_source/node_id_to_ranges"][nrn_rat["node_ids"], :]
edge_ids_rat = h5_rat["edges/default/indices/target_to_source/range_to_edge_id"][ranges_rat[:, 0]]

E_rat = lcm.load_conn_edges(h5_rat, edge_ids_rat, 'rat')

## 3. save flat locations into the dataframe E_rat

xyz = ["afferent_center_x", "afferent_center_y", "afferent_center_z"]

E_rat = lcm.save_flat_loc(E_rat, xyz, fm_rat, orient_rat, pixel_size=pixel_sz_rat)

## 4. Build connectivity object, compress it and save it

M_rat = lcm.build_conn_mat(E_rat, nrn_rat)

# Safe
fname_rat_raw = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Rat_623um_squared_struc_conmat.h5"
M_rat.to_h5(fname_rat_raw)

# Compress and safe
fname_rat_raw = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Rat_623um_squared_struc_conmat_compressed.h5"
M_rat_comp = M_rat.compress()
M_rat_comp.to_h5(fname_rat_raw)

# As many of the synapses will be out of the selected range we need to filter them 
# The ConnectivityMatrix object expects the indices of pre- and post-synaptic neurons to be 
# in the range of 0, 1, 2, ..., N, where N is the size of the population. 
# However, in this case we were giving it the gids of the neurons which can be all over the 
# place (0, 5, 20, 25, 100, ....), so we need to re-index

# filter, re-index, compress for faster loading and safe

fname_r_filt = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix/Rat_623um_squared_struc_conmat_filtered.h5"
M_rat = lcm.filter_and_reindex(M_rat)
M_rat.to_h5(fname_r_filt)

fname_r_comp = "/gpfs/bbp.cscs.ch/project/proj159/home/barros/conn_matrix Rat_623um_squared_struc_conmat_filtered_compressed.h5"
M_r_comp = M_rat.compress()
M_r_comp.to_h5(fname_r_comp)
