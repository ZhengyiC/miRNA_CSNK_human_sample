#this file contains some calculation functions frequently used in scRNA data analysis

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt



def normalization( dat_ct):
    """this function normalize the data so that each cell has the same 
    number of total counts as the median value of the total counts among all cells.
    The data will also be log-like transformed
    Count values will also be transformed to z-scores for each gene"""
    sc.pp.normalize_total(dat_ct) 
    dat_ct.X = np.arcsinh(dat_ct.X).copy()
    sc.pp.scale(dat_ct)
    
    return dat_ct




def save_h5ad_compresseed(dat_ls, name_ls, compression_opt = None):
    """This function save a list of anndata objs with names/paths specified in name_ls
    Files will be saved in compressed form, and the compression option is 'gzip for all files unless specified 
    The function will return the number of files saved successfully'"""
    
    if (not compression_opt):
        compression_opt = ['gzip' for i in range(len(dat_ls))]
    saved_num = 0
    for j in range(len(dat_ls)):
        try:
            dat_ls[j].write(name_ls[j], compression = compression_opt[j])
        except BaseException as err:
            print(err)
            continue
        else:
            saved_num +=1
        
    return saved_num

