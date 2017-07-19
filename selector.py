import pandas as pd, numpy as np
import random
from ranges import get_range_absolute

def get_stats(gene_df, df):
    """
    This function will select the target gene's peaks stats from the dataframe
    :param gene_df: gene:range dictionary, get the best result for each gene from transcripts range
    :param df: dataframe contain peaks from certain cutoff
    :return:
    """
    results = []

    for gene in gene_df['gene'].unique():
        cur_ranges = gene_df[gene_df['gene'] == gene]['range']
        # if gene == 'GFI1':
        #     print gene, cur_ranges
        best_total_width = 0
        best_height = 0
        best_single_width = 0
        best_total_signal = 0
        best_single_signal = 0
        best_skewness = 0
        best_kurtosis = 0
        for r in cur_ranges:
            chr_name, start, end = r
            # print chr_name, start, end
            ## Here is the problem, danpos selector will consider the entire overlapped peaks
            ## The other approach is using self designed peak calling, to make sure each parameter will return different value

            # old slow solution
            cur_df = df[(df['chr']==chr_name) &
                        ((df['start'].between(start, end)) |
                         (df['end'].between(start, end)))]
            cur_total_width = cur_df['width_above_cutoff'].sum()
            cur_height = cur_df['height'].max()
            cur_single_width = cur_df['width_above_cutoff'].max()
            cur_total_signal = cur_df['total_signal'].sum()
            cur_single_signal = cur_df['total_signal'].max()

            # This is for kurtosis and skewness
            if cur_df.shape[0] > 0 and 'skewness' in cur_df.columns and 'kurtosis' in cur_df.columns:
                cur_skewness = cur_df.ix[cur_df['total_signal'].argmax(),'skewness']
                cur_kurtosis = cur_df.ix[cur_df['total_signal'].argmax(),'kurtosis']
            else:
                cur_skewness = 0
                cur_kurtosis = 0

            # new indexing solution
            # cur_df = df.get_peaks(chr_name, start, end)
            # cur_total_width = np.sum([p.width for p in cur_df])
            # cur_height = np.max([p.height for p in cur_df])
            # cur_single_width = np.max([p.width for p in cur_df])
            # cur_total_signal = np.sum([p.total_signal for p in cur_df])
            # cur_single_signal = np.max([p.total_signal for p in cur_df])

            # print cur_total_width, cur_total_signal, cur_height

            if cur_total_width > best_total_width:
                best_total_width = cur_total_width
            if cur_height > best_height:
                best_height = cur_height
            if cur_single_width > best_single_width:
                best_single_width = cur_single_width
            if cur_total_signal > best_total_signal:
                best_total_signal = cur_total_signal
            if cur_single_signal > best_single_signal:
                best_single_signal = cur_single_signal

            # this is for kurtosis and skewness
            if cur_skewness is not None and (best_skewness is None or abs(cur_skewness) > abs(best_skewness)):
                best_skewness = cur_skewness
            if cur_kurtosis is not None and (best_kurtosis is None or abs(cur_kurtosis) > abs(best_kurtosis)):
                best_kurtosis = cur_kurtosis
        results.append((gene, best_total_width, best_height, best_single_width,
                        best_total_signal, best_single_signal,
                        best_skewness, best_kurtosis))
    return results

def random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed):
    """
    random select genes for each cell type

    :param CIG_gene_table: table path
    :param exclude_list: exclude_list path
    :param all_genes: all gene GTF
    :param random_seed: random seed number to make sure each return the same random gene set
    :return:
    """
    results = []
    candidates = set()
    for gene in all_genes:
        if gene not in exclude_list and gene not in CIG_gene_df['gene']:
            candidates.add(gene)
    candidates = list(candidates)
    random.seed(random_seed)
    negative_control_genes = random.sample(candidates, CIG_gene_df.shape[0])

    for i in range(len(negative_control_genes)):
        results.append([negative_control_genes[i]]+list(CIG_gene_df.iloc[i, 1:]))
    negative_control_genes_df = pd.DataFrame(results)
    negative_control_genes_df.columns = CIG_gene_df.columns
    return negative_control_genes_df


def CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, up_stream_distance, widow_size, all_dfs, cutoff):
    """

    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    CIG_gene_list = list(CIG_gene_df['gene'].values)
    non_CIG_gene_list = list(non_CIG_gene_df['gene'].values)
    CIG_gene_ranges = get_range_absolute(CIG_gene_list, all_gene_GTF, up_stream_distance, widow_size)
    non_CIG_gene_ranges = get_range_absolute(non_CIG_gene_list, all_gene_GTF, up_stream_distance, widow_size)

    CIG_results = []
    non_CIG_results = []

    for cell_type in CIG_gene_df['cell_type'].unique():
        cur_CIG_gene_list = list(CIG_gene_df[CIG_gene_df['cell_type'] == cell_type]['gene'].values)
        cur_CIG_gene_range = CIG_gene_ranges[CIG_gene_ranges['gene'].isin(cur_CIG_gene_list)]

        cur_non_CIG_gene_list = list(non_CIG_gene_df[non_CIG_gene_df['cell_type'] == cell_type]['gene'].values)
        cur_non_CIG_gene_range = non_CIG_gene_ranges[non_CIG_gene_ranges['gene'].isin(cur_non_CIG_gene_list)]

        # print cell_type, cutoff
        # print all_dfs[cell_type].keys()
        cur_df = all_dfs[cell_type][cutoff]
        cur_CIG_result = get_stats(cur_CIG_gene_range, cur_df)
        cur_non_CIG_result = get_stats(cur_non_CIG_gene_range, cur_df)

        CIG_results += cur_CIG_result
        non_CIG_results += cur_non_CIG_result
    CIG_results_df = pd.DataFrame(CIG_results)
    # print CIG_results_df
    CIG_results_df.columns = ['gene', 'total_width', 'height', 'single_width', 'total_signal', 'single_signal', 'skewness', 'kurtosis']
    non_CIG_results_df = pd.DataFrame(non_CIG_results)
    non_CIG_results_df.columns = ['gene', 'total_width', 'height', 'single_width', 'total_signal', 'single_signal', 'skewness', 'kurtosis']
    return CIG_results_df, non_CIG_results_df

if __name__ == "__main__":
    from collections import defaultdict
    import os

    ### step 1 get random negative control set
    all_gene_GTF = pd.read_csv('/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()

    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())

    exclude_list_file = open('/home/tmhbxx3/scratch/CIG/test/merge_ten_cells_relative_genes.txt', 'r')
    exclude_list = [x.strip().upper() for x in exclude_list_file]
    exclude_list_file.close()
    exclude_list = set(exclude_list)

    CIG_gene_df = pd.read_csv('CIG_strong_trimmed.csv')
    CIG_gene_df['gene'] = CIG_gene_df['gene'].str.upper()

    # print CIG_gene_df.shape, 'CIG genes'

    non_CIG_gene_df = random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed=9000)

    # print non_CIG_gene_df.shape, 'non CIG gens'

    ### step 2 get all_dfs

    cutoff_range = [15]

    marker = 'h3k27ac'
    dfs_path = '/home/tmhbxx3/scratch/CIG/test/'+ marker + '_sk_peaks/'
    dfs = os.listdir(dfs_path)
    all_dfs = defaultdict(dict)
    for table_name in dfs:
        info = table_name.split('_')
        cell_type = info[0]
        cutoff = int(info[-1][:-4])
        # print cell_type, cutoff
        if cutoff in cutoff_range:
            all_dfs[cell_type][cutoff] = pd.read_csv(dfs_path+table_name)

    CIG_gene_list = list(CIG_gene_df['gene'].values)
    non_CIG_gene_list = list(non_CIG_gene_df['gene'].values)

    print all_dfs.keys()

    df1, df2 = CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                 25000, 25000,
                 all_dfs, 15)

    from wilcoxontest import logP_wilcoxon

    print logP_wilcoxon(df1['kurtosis'], df2['kurtosis'])



# python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py selector /home/tmhbxx3/scratch/CIG/h3k4me3_peaks/pooled/h3k4me3_CD34_h3k4me3_sample_ENCFF357RXZ_ENCFF102IJI_input_ENCFF084HEK_ENCFF088FNX_ENCFF599JOR_sample.bgsub.Fnor.peaks_10.xls --gene_file /home/tmhbxx3/scratch/CIG/test/all_gene_GTF.xls --gene_out ./test_danpos_selector.txt --genicSelector TSS:-3000:3000