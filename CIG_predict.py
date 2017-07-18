import numpy as np
import scipy.optimize as spo
from wilcoxontest import logP_wilcoxon
from selector import CIG_selecter

def optimize_allocs(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_dfs, criteria):
    bounds = ((-500000, 500000), (1000,1000000))
    # bounds =((5000, -5000), (1000,2000,3000,4000), (10,20,30,40))
    options = {'disp': True, 'eps': 1e4, 'ftol': 1e-06}
    # options = {'disp': False}

    # initial guess
    allocs = np.asarray([500000, 10000])
    cur_cutoff = 10

    allocs = spo.minimize(CIG_distance, allocs, args=(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_dfs, criteria, cur_cutoff),
                          method='SLSQP', bounds=bounds, options=options).x

    print allocs, 'this is the final'
    # allocs = allocs/np.sum(allocs)
    return allocs

def CIG_distance(allocs, CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_dfs, criteria, cur_cutoff):
    cur_up_stream_distance, cur_window_size = allocs

    cur_CIG_results_df, cur_non_CIG_results_df = CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                                              cur_up_stream_distance, cur_window_size,
                                                              all_dfs, cur_cutoff)

    # print cur_CIG_results_df.shape, cur_non_CIG_results_df.shape
    if cur_CIG_results_df[criteria].mean() < cur_non_CIG_results_df[criteria].mean():
        cur_logP = logP_wilcoxon(cur_CIG_results_df[criteria],
                                 cur_non_CIG_results_df[criteria])
    else:
        cur_logP = logP_wilcoxon(cur_non_CIG_results_df[criteria],
                                 cur_CIG_results_df[criteria])

    print cur_up_stream_distance, cur_window_size, cur_cutoff, cur_logP
    return cur_logP