# Author: Daniela Egas Santander
# Last updated: 10.2024

# Functions to prune connectomes (i.e., weighted matrices) in different ways


import scipy.sparse as sp
import numpy as np 

def random_sub_coo_matrix(n, data, row, col,n_target, seed):
    """Utility function that randomly samples data, row and column to the target size n_target
    """
    np.random.seed(seed) # To generate reproducible results
    N=data.size
    assert n_target<=N, "The size of the subsample must be smaller or equal than the total size"
    sample_idx=np.sort(np.random.choice(np.arange(N), size=n_target, replace=False))
    A_prunned=sp.coo_matrix((data[sample_idx], (row[sample_idx], col[sample_idx])), shape=(n, n))
    A_prunned.sum_duplicates()
    return A_prunned.astype(int)

def prune_touches_to_target_touches(A, n_touches, seed): 
    """Randomly prunes the touches of connectome A until it has n_touches in total.  
    Parameters
    ----------
    A : sparse matrix 
    n_touches: target number of touches i.e., target sum of the entries 
    seed: seed for randomization 
    
    Returns
    -------
    sparse matrix
        With entries natural numbers smaller or equal than those in A and whose sum is n_touches.
    """
    # Get data from A
    A=A.tocoo()
    N=A.shape[0]
    row=A.row
    col=A.col
    data=A.data
    # Transform to one with multiple entries according to the number of touches
    new_row, new_col = [], []
    for i, val in enumerate(data):
        new_row.extend([row[i]] * val)
        new_col.extend([col[i]] * val)
    # Create a new data array where all entries are 1
    new_data = np.ones(len(new_row))
    new_row=np.array(new_row)
    new_col=np.array(new_col)
    
    # Returned subsampled matrix
    return random_sub_coo_matrix(N, new_data, new_row, new_col,n_touches, seed)
    
def prune_edges_to_target_edges(A, n_edges, seed): 
    """Randomly prunes the edges of connectome A until it has n_edges in total.  
    Parameters
    ----------
    A : sparse matrix 
    n_edges: target number of edges 
    seed: seed for randomization 
    
    Returns
    -------
    sparse matrix
        Submatrix of A with n_edges
    """
    # Get data from A
    A=A.tocoo()
    N=A.shape[0]
    row=A.row
    col=A.col
    data=A.data
    
    # Returned subsampled matrix
    return random_sub_coo_matrix(N,data, row, col,n_edges, seed)

def prune_touches_to_target_edges(A, n_edges, seed, tol=0.01):
    """Randomly prunes touches of connectome A until it has (approximately) n_edges in total.  
    Parameters
    ----------
    A : sparse matrix 
    n_edges: target number of edges 
    seed: seed for randomization
    tol: error tolerance for the number of edges (fraction of the total)
    
    Returns
    -------
    sparse matrix
        Submatrix of A with approximately n_edges
    """
    #Iteratively prune touches of a to approximate functional density 
    A_touches=A.copy()
    counter=0
    while A_touches.nnz/n_edges>1+tol:
        prunning_guess=int(A_touches.sum()*n_edges/A_touches.nnz)
        A_touches=prune_touches_to_target_touches(A_touches,prunning_guess, seed)
        if counter%5==0: 
            print(f"Ran {counter}-loop.  Number of edges in loop {A_touches.nnz}")
        counter+=1
    print(f"Ran {counter} loops.  The result has {A_touches.nnz*100/n_edges:.2f}% of the target edges")
    return A_touches


def prune_edges_to_target_touches(A, n_touches, seed, tol=0.01):
    """Randomly prunes the edges of connectome A until it has n_touches in total.  
    Parameters
    ----------
    A : sparse matrix 
    n_touches: target number of touches i.e., target sum of the entries 
    seed: seed for randomization
    tol: error tolerance for the number of touches (fraction of the total)
    
    Returns
    -------
    sparse matrix
        With entries natural numbers smaller or equal than those in A and whose sum is approximately n_touches.
    """
    A_edges=A.copy()
    counter=0
    while A_edges.sum()/n_touches>1+tol:
        mean_weight=A_edges.sum()/A_edges.nnz
        guess_n_edges=np.ceil(n_touches/mean_weight).astype(int)
        A_edges=prune_edges_to_target_edges(A, guess_n_edges, seed=seed)
        if counter%5==0: 
            print(f"Ran {counter}-loop.  Number of touches in loop {A_edges.sum()}")
        counter+=1
    print(f"Ran {counter} loops.  The result has {A_edges.sum()*100/n_touches:.2f}% of the target edges")
    return A_edges
