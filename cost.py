from selector import CIG_selecter, CIG_selecter_all
from stat_test import logP_wilcoxon, logP_fisher

def wilcoxon_cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                          cur_up_stream_distance, cur_down_stream_distance,
                          all_dfs, cur_cutoff, criteria, marker,
                          TSS_pos, TTS_pos, wigs):
    cur_CIG_results_df, cur_non_CIG_results_df = CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                                              cur_up_stream_distance, cur_down_stream_distance,
                                                              all_dfs, cur_cutoff, criteria,
                                                              TSS_pos, TTS_pos,
                                                              wigs)

    # print 'wilcoxon'
    print cur_CIG_results_df[criteria].mean(),  cur_non_CIG_results_df[criteria].mean(), 'average'

    if marker == 'h3k27me3' and criteria != 'kurtosis' and criteria != 'skewness':
        cur_logP = logP_wilcoxon(cur_CIG_results_df[criteria],
                                 cur_non_CIG_results_df[criteria])
    elif criteria == 'kurtosis' or criteria == 'skewness':
        if cur_CIG_results_df[criteria].mean() < cur_non_CIG_results_df[criteria].mean():
            cur_logP = logP_wilcoxon(cur_CIG_results_df[criteria],
                                     cur_non_CIG_results_df[criteria])
        else:
            cur_logP = logP_wilcoxon(cur_non_CIG_results_df[criteria],
                                     cur_CIG_results_df[criteria])
    else:
        cur_logP = logP_wilcoxon(cur_non_CIG_results_df[criteria],
                                 cur_CIG_results_df[criteria])
    return cur_logP

def fisher_cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                         cur_up_stream_distance, cur_down_stream_distance,
                         all_dfs, cur_cutoff, criteria, marker,
                         TSS_pos, TTS_pos, wigs):
    all_gene_results_df = CIG_selecter_all(CIG_gene_df, all_gene_GTF, cur_up_stream_distance, cur_down_stream_distance,
                                           all_dfs, cur_cutoff, criteria,
                                           TSS_pos, TTS_pos, wigs)
    # print 'fisher'
    cur_logP = logP_fisher(CIG_gene_df, all_gene_results_df, criteria, top_enrich=500)

    return cur_logP