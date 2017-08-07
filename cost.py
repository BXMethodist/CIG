from selector import CIG_selecter, CIG_selecter_all
from stat_test import logP_wilcoxon

def wilcoxon_cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                      cur_up_stream_distance, cur_window_size,
                      all_dfs, cur_cutoff, criteria):
    cur_CIG_results_df, cur_non_CIG_results_df = CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                                              cur_up_stream_distance, cur_window_size,
                                                              all_dfs, cur_cutoff, criteria)

    # print cur_CIG_results_df.shape, cur_non_CIG_results_df.shape
    if cur_CIG_results_df[criteria].mean() < cur_non_CIG_results_df[criteria].mean():
        cur_logP = logP_wilcoxon(cur_CIG_results_df[criteria],
                                 cur_non_CIG_results_df[criteria])
    else:
        cur_logP = logP_wilcoxon(cur_non_CIG_results_df[criteria],
                                 cur_CIG_results_df[criteria])
    return cur_logP

def fisher_cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                      cur_up_stream_distance, cur_window_size,
                      all_dfs, cur_cutoff, criteria):
    all_gene_results_df = CIG_selecter_all(all_gene_GTF, cur_up_stream_distance, cur_window_size,
                                                              all_dfs, cur_cutoff, criteria)

    return cur_logP