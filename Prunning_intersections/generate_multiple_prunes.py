# Author: Daniela Egas Santander
# Last updated: 02.2024

# Generating multiple prunned controls

from pathlib import Path
import pickle
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


def generate_prunned_controls(animal, mat_types, seeds):

    # Load functional matrices 
    fname_func=fnames["functional"][animal]
    A_functional=ConnectivityMatrix.from_h5((root/fname_func)).matrix.tocsr()
    
    # Load structural data 
    fname_struct=fnames["structural"][animal]
    A_structural=ConnectivityMatrix.from_h5(root/fname_struct).matrix.tocsr()
    
    # For testing we use a small submatrix, comment out the next two lines to generate the full matrices
    #A_structural=A_structural[np.ix_(range(100), range(100))] #just for testing  
    #A_functional=A_functional[np.ix_(range(100), range(100))] #just for testing
    
    #Target values from the functional matrices
    target_edges=A_functional.nnz
    target_touches=A_functional.sum()
    
    print(f"Number of edges in structural {A_structural.nnz}")
    print(f"Target number of edges {target_edges}")
    
    print(f"\nNumber of touches in structural {A_structural.sum()}")
    print(f"Target number of touches {target_touches}")

    
    # Prunning to match touches 
    for prunning_name, prunning in mat_types.items():
        mats={}
        target = target_touches if prunning_name[-4:]=="to_T" else target_edges
        print(prunning_name, target)
        for seed in seeds:
            mats[seed]=prunning(A_structural, target, seed=seed)
        path_out=f"{root}/mat_intersection/{prunning.__name__}_multiple_{animal}.pkl"
        with open(path_out, "wb") as f:
            pickle.dump(mats, f)
        print(f"Done with {prunning.__name__} for {animal}")

if __name__ == "__main__":
    seeds = np.arange(30) 

    mat_types={"T_to_T":prune_touches_to_target_touches,
               "E_to_T":prune_edges_to_target_touches,
               "T_to_E":prune_touches_to_target_edges,
               "E_to_E":prune_edges_to_target_edges
              }
    for animal in ["rat", "human"]:
        print(f"Running analysis for {animal}")
        generate_prunned_controls(animal, mat_types, seeds)