# Author: Daniela Egas Santander
# Last updated: 02.2024

# Computing size of the intersections

from pathlib import Path
import pickle
import itertools
import scipy.sparse as sp
import numpy as np 
from conntility import ConnectivityMatrix
import sys
sys.path.append("/Users/danielaegas/repos/Human_connectivity_paper/matrix_generation_extraction")
from prunning import *

# Paths to data
root=Path("/Users/danielaegas/")/"OneDrive - Open Brain Institute"/"Human_Connectivity_Paper"/"conn_matrix"
fnames= {"structural":{"rat":
                       "Rat_623um_squared_struc_conmat_filtered_compressed.h5",
                       "human":
                       "Human_NEWstruct_conmat_filtered_compressed.h5"},
         "functional":{"rat":
                       "Rat_623um_squared_funct_conmat_filtered_compressed.h5",
                       "human":
                       "Human_NEWfunct_conmat_filtered_compressed.h5"}}


def int_functional(animal, mat_types):
    print("Computing union with functional connectome")
    # Load functional matrix 
    fname_func=fnames["functional"][animal]
    A_functional=ConnectivityMatrix.from_h5((root/fname_func)).matrix
    intersection_functional={prunning_name:{} for prunning_name in mat_types.keys()}
    for prunning_name, prunning in mat_types.items():
        path_in=f"{root}/mat_intersection/{prunning.__name__}_multiple_{animal}.pkl"
        with open(path_in, "rb") as f:
            mats = pickle.load(f)
        for seed in mats.keys():
            intersection_functional[prunning_name][seed]=(mats[seed].astype("bool")+A_functional.astype("bool")).nnz
        print(f"Done with {prunning_name}")
    path_out=f"{root}/mat_intersection/union_with_funct_{animal}.pkl"
    with open(path_out, "wb") as f:
        pickle.dump(intersection_functional, f)
    print(f"Done with unions with functional for {animal}")

   
def int_controls(animal, mat_types, reduced=False):
    prunning_types=list(mat_types.keys())
    print("Computing union across controls")
    # Intersection size between controls 
    intersections={}
    for i, j in itertools.combinations_with_replacement(prunning_types, 2):
        intersections[i,j]={}
        path_i=f"{root}/mat_intersection/{mat_types[i].__name__}_multiple_{animal}.pkl"
        with open(path_i, "rb") as f:
            mats_i = pickle.load(f)
        seeds=mats_i.keys()
        if i==j:
            l=0
            for s in seeds:
                intersections[(i,j)][(0,s)]=(mats_i[l].astype("bool")+mats_i[s].astype("bool")).nnz
        if reduced: 
            if i!=j:
                path_j=f"{root}/mat_intersection/{mat_types[j].__name__}_multiple_{animal}.pkl"
                with open(path_j, "rb") as f:
                    mats_j = pickle.load(f)
                l=0
                for s in seeds:
                    intersections[(i,j)][(l,s)]=(mats_i[l].astype("bool")+mats_j[s].astype("bool")).nnz
        else: 
            if i!=j:
                path_j=f"{root}/mat_intersection/{mat_types[j].__name__}_multiple_{animal}.pkl"
                with open(path_i, "rb") as f:
                    mats_j = pickle.load(f)
                for l,s in itertools.combinations_with_replacement(seeds, 2):
                    intersections[(i,j)][(l,s)]=(mats_i[l].astype("bool")+mats_j[s].astype("bool")).nnz
        print(f"Done with {i}, {j}")
    path_out=f"{root}/mat_intersection/union_across_controls_{animal}.pkl"
    with open(path_out, "wb") as f:
        pickle.dump(intersections, f)
    print(f"Done with unions across controls {animal}")
   




if __name__ == "__main__":
    mat_types={"T_to_T":prune_touches_to_target_touches,
               "E_to_T":prune_edges_to_target_touches,
               "T_to_E":prune_touches_to_target_edges,
               "E_to_E":prune_edges_to_target_edges
              }
    for animal in ["rat", "human"]:
        print(f"Analyzing {animal}")
        int_functional(animal, mat_types)
        int_controls(animal, mat_types, reduced=True)























        