"""Compute simplex counts for functional matrices and controls
Last updated: 05.2025
Author: Daniela Egas
"""

from pathlib import Path
import pickle
import numpy as np 
import pandas as pd
from conntility import ConnectivityMatrix
from connalysis.network.topology import simplex_counts
from connalysis.randomization import ER_model, configuration_model, run_DD2, bishuffled_model
import time
# Paths to data
animals=["rat", "human"]
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
mat_types={"prune_touches_to_target_edges":"T_to_E",
           "prune_edges_to_target_edges":"E_to_E"}

# Setting up controls 
# Loading optimized DD2 parameters - precomputed 

DD2_params={"human":{"a":0.33377772104625897, "b":0.008021800520398505}, 
            "rat":{"a":0.13442783930481592, "b":0.00838843470489593}}
for animal in DD2_params.keys():
    connectome =ConnectivityMatrix.from_h5((root/fnames["functional"][animal]))
    DD2_params[animal]["n"]=connectome._shape[0]
    DD2_params[animal]["xyz"]=connectome.vertices[["x","y","z"]].to_numpy()

# Control dict of controls to run simplex counts on
no_seeds=30
controls_config={"ER":{"function":ER_model, 
                       "seeds":[(i,i) for i in np.arange(no_seeds)],
                        },
                "CM":{"function":configuration_model,
                      "seeds":np.arange(no_seeds)
                        },
                "bishuffuled":{"function":bishuffled_model,
                               "seeds":np.arange(no_seeds)}, 
                "DD2":{"function":run_DD2,
                        "seeds":[(i,i) for i in np.arange(no_seeds)]
                        },
                "E_to_E":{"function":"precomputed",
                          "seeds":np.arange(no_seeds)
                          },
                "T_to_E":{"function":"precomputed",
                          "seeds":np.arange(no_seeds)
                          },
                        
}
# Compute simplex counts on functional matrices and all controls 
for animal in animals:
        start_global=time.time() 
        print(f"Running {animal}")
        # Load data
        # Functional matrix
        fname_func=fnames["functional"][animal]
        connectome=ConnectivityMatrix.from_h5((root/fname_func))
        A_functional=connectome.matrix.astype(bool).astype(int)
        #Load pregenerated control matrices
        mats={}
        for path_name,prunning__name in mat_types.items():
                print(f"Loading {prunning__name} for {animal}")
                path_in=f"{root}/mat_intersection/{path_name}_multiple_{animal}.pkl"
                with open(path_in, "rb") as f:
                        mats[prunning__name] = pickle.load(f)
        # Computing simplex counts
        sc = {}
        # Compute original simplex counts 
        start=time.time()
        sc["functional"]=simplex_counts(A_functional)
        print(f"Functional simplex counts took {time.time()-start} seconds")

        # Compute simplex counts for each control
        for model_name, model_config in controls_config.items():
                start=time.time()
                print(f"Running {model_name}")
                sc[model_name]={}
                func = model_config["function"]
                for seed in model_config["seeds"]:
                        if model_name=="DD2":
                                adj_control = run_DD2(DD2_params[animal]["n"],DD2_params[animal]["a"],DD2_params[animal]["b"], 
                                                DD2_params[animal]["xyz"], seed=seed)
                        elif model_name in ["E_to_E", "T_to_E"]:
                                try:
                                        adj_control = mats[model_name][seed]
                                except KeyError as e:
                                        print(f"KeyError: {e} - Seed {seed} not found for model {model_name}")
                                        continue
                                
                        else:
                                adj_control = func(A_functional, seed=seed)
                        sc[model_name][seed] = simplex_counts(adj_control)
                print(f"{model_name} simplex counts took {(time.time()-start)/60:.2f} minutes, for {len(model_config['seeds'])} seeds")

        # Forma and save results 
        sc_formatted={}
        for key in sc.keys():
                if key =="functional":
                        sc_formatted[key]=sc[key]
                else:
                        sc_formatted[key]=pd.DataFrame.from_dict(sc[key], orient="index")

                # Save the dictionary to a pickle file
                path_out=f"{root}/simplex_counts/{animal}_full_sc.pkl"
                with open(path_out, 'wb') as f:
                        pickle.dump(sc_formatted, f)
        print(f"Done with {animal} in {(time.time()-start_global)/60:.2f} minutes")
        print("========================================")