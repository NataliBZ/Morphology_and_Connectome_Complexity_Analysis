from pathlib import Path
import pickle
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from conntility import ConnectivityMatrix

import sys
sys.path.append("/Users/danielaegas/repos/Human_connectivity_paper/matrix_generation_extraction")
from prunning import *

def edge_survival_probability_for_touch_count(m_before, m_after):
    """
    Returns for each touch count the probability that the edge survives pruning.
    Args:
      m_before: sparse matrix of touches per pair before pruning
      m_after: sparse matrix of touches per pair after pruning
    """
    m_before = m_before.tocoo(); m_after = m_after.tocoo()
    df1 = pd.DataFrame(
        {"row": m_before.row, "col": m_before.col, "count": m_before.data}
    ).set_index(["row", "col"])
    df2 = pd.DataFrame(
        {"row": m_after.row, "col": m_after.col, "active": m_after.data > 0}
    ).set_index(["row", "col"])
    df = pd.concat([df1, df2], axis=1).fillna(False)
    return df.groupby("count")["active"].mean()

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

# Prunning types 
mat_types={"T_to_T":prune_touches_to_target_touches,
           "E_to_T":prune_edges_to_target_touches,
           "T_to_E":prune_touches_to_target_edges,
           "E_to_E":prune_edges_to_target_edges
          }




animals=["rat", "human"]
survival_prob={}
for animal in animals:
    # Load functional and structural matrices
    fname_func=fnames["functional"][animal]
    A_functional=ConnectivityMatrix.from_h5((root/fname_func)).matrix
    fname_func=fnames["structural"][animal]
    A_structural=ConnectivityMatrix.from_h5((root/fname_func)).matrix
    
    #Load control matrices
    mats={}
    for prunning_name, prunning in mat_types.items():
        path_in=f"{root}/mat_intersection/{prunning.__name__}_multiple_{animal}.pkl"
        with open(path_in, "rb") as f:
            mats[prunning_name] = pickle.load(f)
    
    print(f"Analyzing {animal}")
    survival_prob[animal]={}
    survival_prob[animal]["functional"]=edge_survival_probability_for_touch_count(A_structural,A_functional)
    print(f"Done with functional") 
    for prunning_name in mats.keys():
        survival_prob[animal][prunning_name]={}
        for seed in mats["T_to_T"].keys(): 
            survival_prob[animal][prunning_name][seed]=edge_survival_probability_for_touch_count(A_structural,mats[prunning_name][seed])
        print(f"Done with {prunning_name}")   

path_out=f"{root}/mat_intersection/survival_probabilites.pkl"
with open(path_out, "wb") as f:
    pickle.dump(survival_prob, f)
print(f"Done")