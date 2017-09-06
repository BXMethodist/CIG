import pandas as pd, numpy as np
import random
from ranges import get_range_absolute
from collections import defaultdict

def df_to_index_danpos(df, bin=3000):
    # print df
    results = {}
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.split()
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def df_to_index_sk(df, bin=3000):
    results = {}
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.strip().split(',')
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             float(line[8]),
             float(line[9]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()
    return results

def get_stats(gene_df, df_path, criteria, bin=3000, df_function=df_to_index_danpos):
    """
    This function will select the target gene's peaks stats from the dataframe
    :param gene_df: gene:range dictionary, get the best result for each gene from transcripts range
    :param df: dataframe contain peaks from certain cutoff
    :return:
    """
    if criteria != 'skewness' and criteria != 'kurtosis':
        table_dict = df_function(df_path)
    else:
        df_function = df_to_index_sk
        table_dict = df_function(df_path)

    results = defaultdict(float)

    for k in range(gene_df.shape[0]):

        # print cur_ranges
        # if gene == 'GFI1':
        #     print gene, cur_ranges
        gene_name = gene_df.iloc[k, 0]

        chr_name, start, end = gene_df.iloc[k, 1], gene_df.iloc[k, 2], gene_df.iloc[k, 3]
        ## Here is the problem, danpos selector will consider the entire overlapped peaks
        ## The other approach is using self designed peak calling, to make sure each parameter will return different value
        cur_table = set()

        if end < start:
            mid = (start + end) / 2
            start = mid
            end = mid

        for i in range(int(start/bin), int(end/bin) + 1):
            if chr_name in table_dict and i in table_dict[chr_name]:
                table = table_dict[chr_name][i]
                cur_table = cur_table.union(table)

        if len(cur_table) == 0:
            continue

        # if gene_name == 'A2BP1':
        #     print cur_table

        selected_table = []
        for t in cur_table:
            if start < t[1] < end:
                selected_table.append(t)
            elif start < t[2] < end:
                selected_table.append(t)
            elif t[1] <= start and end <= t[2]:
                selected_table.append(t)

        # if gene_name == 'CCL13':
        #     print selected_table, start, end

        if len(selected_table) == 0:
            continue

        cur_df = pd.DataFrame(list(selected_table))

        if cur_df.shape[1] == 6:
            cur_df.columns = ['chr',
                          'start',
                          'end',
                          'width_above_cutoff',
                          'total_signal',
                          'height',]
        else:
            cur_df.columns = ['chr',
                              'start',
                              'end',
                              'width_above_cutoff',
                              'total_signal',
                              'height',
                              'skewness',
                              'kurtosis']

        if criteria == 'total_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.sum()
        elif criteria == 'height':
            cur_value = cur_df['height'].max()
        elif criteria == 'single_width':
            cur_col = cur_df['end'] - cur_df['start']
            cur_value = cur_col.max()
        elif criteria == 'total_signal':
            cur_value = cur_df['total_signal'].sum()
        elif criteria == 'single_signal':
            cur_value = cur_df['total_signal'].max()
        #
        # # This is for kurtosis and skewness
        elif cur_df.shape[0] > 0 and criteria == 'skewness' and 'skewness' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(),'skewness']
        elif cur_df.shape[0] > 0 and criteria == 'kurtosis' and 'kurtosis' in cur_df.columns:
            cur_value = cur_df.ix[cur_df['total_signal'].argmax(), 'kurtosis']

        if cur_value > results[gene_name] and criteria != 'skewness' and criteria != 'kurtosis':
            results[gene_name] = cur_value
        # this is for kurtosis and skewness
        elif criteria == 'kurtosis':
            if abs(cur_value) > abs(results[gene_name]):
                results[gene_name] = cur_value
        elif criteria == 'skewness':
            if abs(cur_value) > results[gene_name]:
                results[gene_name] = abs(cur_value)

    final = []
    for gene_name in gene_df['gene'].unique():
        final.append((gene_name, results[gene_name]))
    # print final, len(final)
    return final

def call_stats(gene_df, wigs, criteria, cutoff):
    results = defaultdict(float)

    for k in range(gene_df.shape[0]):
        gene_name = gene_df.iloc[k, 0]
        chr_name, start, end = gene_df.iloc[k, 1], gene_df.iloc[k, 2], gene_df.iloc[k, 3]
        if chr_name in wigs.genome:
            cur_signals = wigs.genome[chr_name].get_signals(start, end)
            # criterias = ['total_width']
            # criterias = ['single_width']
            # criterias = ['height']
            # criterias = ['total_signal']
            # criterias = ['single_signal']
            # criterias = ['skewness']
            # criterias = ['kurtosis']
            if criteria == 'total_width':
                cur_value = len(np.where(cur_signals > cutoff)[0])
                if results[gene_name] < cur_value:
                    results[gene_name] = cur_value
            elif criteria == 'total_width':
                pass
            elif criteria == 'single_width':
                pass
            elif criteria == 'height':
                pass
            elif criteria == 'total_signal':
                pass
            elif criteria == 'single_signal':
                pass

        else:
            continue

    final = []
    for gene_name in gene_df['gene'].unique():
        final.append((gene_name, results[gene_name]))
    return final


def random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed, number_genes=None):
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

    if number_genes is None:
        number_genes = len(candidates)

    negative_control_genes = random.sample(candidates, number_genes)

    for i in range(len(negative_control_genes)):
        try:
            results.append([negative_control_genes[i]]+list(CIG_gene_df.iloc[i, 1:]))
        except:
            results.append([negative_control_genes[i]] + list(CIG_gene_df.iloc[0, 1:]))
    negative_control_genes_df = pd.DataFrame(results)
    negative_control_genes_df.columns = CIG_gene_df.columns
    return negative_control_genes_df

def CIG_selecter(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, criteria,
                 TSS_pos, TTS_pos, wigs):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
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
    CIG_gene_ranges = get_range_absolute(CIG_gene_list, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos,
                                         TTS_pos)
    non_CIG_gene_ranges = get_range_absolute(non_CIG_gene_list, all_gene_GTF, up_stream_distance, down_stream_distance,
                                             TSS_pos,
                                             TTS_pos)

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
        # print cur_CIG_gene_range
        cur_CIG_result = get_stats(cur_CIG_gene_range, cur_df, criteria)
        cur_non_CIG_result = get_stats(cur_non_CIG_gene_range, cur_df, criteria)

        CIG_results += cur_CIG_result
        non_CIG_results += cur_non_CIG_result
    CIG_results_df = pd.DataFrame(CIG_results)
    # print CIG_results_df
    CIG_results_df.columns = ['gene', criteria]
    non_CIG_results_df = pd.DataFrame(non_CIG_results)
    non_CIG_results_df.columns = ['gene', criteria]
    return CIG_results_df, non_CIG_results_df

def CIG_selecter_all(CIG_gene_df, all_gene_GTF, up_stream_distance, down_stream_distance, all_dfs, cutoff, criteria,
                     TSS_pos, TTS_pos, wigs):
    """
    get the genes status and return a data frame with two columns, gene name and criteria.
    :param CIG_gene_df:
    :param non_CIG_gene_df:
    :param all_gene_GTF:
    :param up_stream_distance:
    :param widow_size:
    :param all_dfs: dictionary of dictionary, cell type and cutoff
    :param cutoff:
    :return:
    """
    all_gene_ranges = get_range_absolute(None, all_gene_GTF, up_stream_distance, down_stream_distance,
                                         TSS_pos, TTS_pos)

    all_gene_results = []

    for cell_type in CIG_gene_df['cell_type'].unique():
        cur_df = all_dfs[cell_type][cutoff]

        # overlap selecter
        cur_CIG_result = get_stats(all_gene_ranges, cur_df, criteria)

        # non overlap selecter
        # cur_CIG_result = call_stats(all_gene_ranges, wigs, criteria, cutoff)

        all_gene_results += cur_CIG_result

    all_gene_results_df = pd.DataFrame(all_gene_results)

    all_gene_results_df.columns = ['gene', criteria]

    return all_gene_results_df

if __name__ == "__main__":
    import os
    all_gene_GTF = pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
    CIG_gene_df = pd.read_csv('top500.tuson.oncogene_keji.csv')
    dfs_path = '/home/tmhbxx3/archive/h3k27me3/' + 'h3k27me3' + '_regions/pooled/'
    dfs = [x for x in os.listdir(dfs_path) if x.endswith('.xls')]
    all_dfs = defaultdict(dict)
    for table_name in dfs:
        info = table_name.split('_')


        cell_type = info[1]
        # print info
        # print info[-1]
        cutoff = float(info[-1][:-4])
        # print cell_type, cutoff
        if cutoff in [10]:
            # print 'start load table', table_name
            # all_dfs[cell_type][cutoff] = peak_table(dfs_path+table_name)

            all_dfs[cell_type][cutoff] = dfs_path + table_name

    result = CIG_selecter_all(CIG_gene_df, all_gene_GTF, -5000, 5000, all_dfs, 10, 'total_width',
                     'TSS', 'TSS', None)

    result.to_csv('bo_selector.tsv', sep='\t', index=False)




# ('chr17', 32688470, 32689280, 590, 12239.2809126, 41.2235371253)