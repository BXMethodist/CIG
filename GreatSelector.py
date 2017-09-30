import os
import pandas as pd, numpy as np
from Wig import Wig
from selector import random_control_genes
from stat_test import logP_wilcoxon
from collections import defaultdict


def danpos_Great():
    pbs = open("Great_selector.pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N Great_selector\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")

    markers =['h3k4me3', 'h3k4me1', 'h3k27ac', 'h3k27me3']

    for marker in markers:
        cur_folder = '/home/tmhbxx3/scratch/CIG/'+marker+'default_peaks/pooled/'
        names = [x for x in os.listdir(cur_folder) if x.endswith('_10.xls')]
        tables = [cur_folder + x for x in names]
        for i in range(len(names)):
            table = tables[i]
            name = names[i]
            pbs.write('python /archive/tmhbxx3/tools/danposTemp_multi_q/danpos.py selector ' + table + ' --gene_file /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls --out ' + './'+name+'out.xls' + ' --gene_out ./'+name+'gene_out.xls --GREATSelector GREAT:-5000:1000:1000000\n')

def load_danpos_table(path, wig):
    f = open(path, 'r')
    info = f.readlines()[1:]
    f.close()
    results = {}
    kur_results = []
    skew_results = []
    for line in info:
        cur_info = line.split()
        chr, start, end = cur_info[0:3]
        # print cur_info
        if chr in wig.genome.keys():
            skewness, kurtosis = wig.genome[chr].get_sk(int(start), int(end))
        else:
            skewness, kurtosis = 0, 0

        results[tuple(cur_info[0:4])] = {'width': int(cur_info[4]),
                                         'signal': float(cur_info[5]),
                                         'height': float(cur_info[6]),
                                         'kurtosis':kurtosis,
                                         'skewness':skewness}
        kur_results.append(kurtosis)
        skew_results.append(skewness)
    path_df = pd.read_csv(path, sep='\t')
    path_df['kurtosis'] = kur_results
    path_df['skewness'] = skew_results
    path_df.to_csv(path, sep='\t', index=None)
    return results

def get_genes_peakfeatures(peaks, selects):
    df = pd.read_csv(selects, sep='\t')
    total_width = []
    single_width = []
    total_signal = []
    single_signal = []
    height = []
    kurtosis = []
    skewness = []
    for i in range(df.shape[0]):
        selected = df.ix[i, 'selections']
        cur_total_width = 0
        cur_total_signal = 0
        cur_single_signal = 0
        cur_single_width = 0
        cur_height = 0
        cur_skewness = 0
        cur_kurtosis = 0
        selected = selected.replace('TSS:', '').split(',')

        for p in selected:
            p = tuple(p.split('/'))

            cur_total_width += int(peaks[p]['width'])
            cur_total_signal += float(peaks[p]['signal'])
            ss = float(peaks[p]['signal'])
            if ss > cur_single_signal:
                cur_single_signal = ss
                cur_kurtosis = peaks[p]['kurtosis']
                cur_skewness = peaks[p]['skewness']
            sw = int(peaks[p]['width'])
            if sw > cur_single_width:
                cur_single_width = sw
            h = float(peaks[p]['height'])
            if h > cur_height:
                cur_height = h


        total_width.append(cur_total_width)
        total_signal.append(cur_total_signal)
        single_signal.append(cur_single_signal)
        single_width.append(cur_single_width)
        height.append(cur_height)
        kurtosis.append(cur_kurtosis)
        skewness.append(cur_skewness)

    df['total_width'] = total_width
    df['total_signal'] = total_signal
    df['single_width'] = single_width
    df['single_signal'] = single_signal
    df['height'] = height
    df['kurtosis'] = kurtosis
    df['skewness'] = skewness
    df.to_csv(selects, index=None, sep='\t')
    return df

def gene_features(table, gtf='/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'):
    pass

def get_gene_transcript(gtf='/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls'):
    df = pd.read_csv(gtf, sep='\t')
    from collections import defaultdict
    results = defaultdict(set)
    for i in range(df.shape[0]):
        results[df.ix[i, 'hg19.kgXref.geneSymbol']].add(df.ix[i, '#hg19.knownGene.name'])
    return results

if __name__ == "__main__":
    ## Get the features of defautl Great Algo -5000:1000:1000000
    # markers =['h3k4me3', 'h3k4me1', 'h3k27ac', 'h3k27me3']
    # for marker in markers:
    #     cur_folder = '/home/tmhbxx3/scratch/CIG/' + marker + 'default_peaks/pooled/'
    #     names = [x for x in os.listdir(cur_folder) if x.endswith('_10.xls')]
    #     tables = [cur_folder + x for x in names]
    #     for i in range(len(names)):
    #         table = tables[i]
    #         name = names[i]
    #         wig_name = name[name.find('_')+1:name.find('.peaks')]
    #         wig_path = '/home/tmhbxx3/scratch/CIG/' + marker + '/' + wig_name + '.wig'
    #         cur_wig = Wig(wig_path)
    #         cur_peaks = load_danpos_table(table, cur_wig)
    #
    #         out_name = name + 'gene_out.xls'
    #         cur_result_df = get_genes_peakfeatures(cur_peaks, out_name)
    #         cur_result_df.to_csv(name[:-4]+'_features.csv', index=None, sep='\t')
    #
    # for marker in markers:
    #     os.system('mv *'+marker+'*_features.csv '+'/home/tmhbxx3/scratch/CIG/GreatAlgoDefault/'+marker +'/')

    ### Get CIG and non-CIG features of Great Algo -5000:1000:1000000
    all_gene_GTF = pd.read_csv('/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls', sep='\t')

    # all_gene_GTF = pd.read_csv('hg19.GREATgene2UCSCknownGenes.table.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()

    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())

    exclude_list_file = open('/home/tmhbxx3/scratch/CIG/test/merge_ten_cells_relative_genes.txt', 'r')
    exclude_list = [x.strip() for x in exclude_list_file]
    exclude_list_file.close()
    exclude_list = set(exclude_list)

    CIG_gene_df = pd.read_csv('CIG_strong_bo.csv')
    # CIG_gene_df = pd.read_csv('top500.tuson.oncogene_keji.csv')
    # CIG_gene_df = pd.read_csv('top500.tuson.tumorSuppressor_keji.csv')
    CIG_gene_df['gene'] = CIG_gene_df['gene'].str.upper()

    # print CIG_gene_df.shape, 'CIG genes'

    non_CIG_gene_df = random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed=9000, number_genes=CIG_gene_df.shape[0])
    # non_CIG_gene_df = random_control_genes(CIG_gene_df, [], all_genes, random_seed=9000, number_genes=None)
    # non_CIG_gene_df = pd.read_csv('top500.tuson.tumorSuppressor_keji.csv')
    non_CIG_gene_df['gene'] = non_CIG_gene_df['gene'].str.upper()

    gene_to_transcript = get_gene_transcript()

    markers = ['h3k4me3', 'h3k4me1', 'h3k27me3', 'h3k27ac']
    features = ['total_width', 'total_signal', 'height', 'skewness', 'kurtosis', 'single_signal', 'single_width', ]

    all_dfs = defaultdict(dict)
    for marker in markers:
        dfs_path = '/home/tmhbxx3/scratch/CIG/GreatAlgoDefault/'+marker +'/'
        dfs = os.listdir(dfs_path)

        for table_name in dfs:
            info = table_name.split('_')
            cell_type = info[1]
            all_dfs[marker][cell_type] = dfs_path + table_name

    cig_results = []
    non_cig_results = []


    for marker in markers:
        cur_cig_results = []
        cur_non_cig_results = []
        for cell_type in CIG_gene_df['cell_type'].unique():
            cur_CIG_genes = CIG_gene_df[CIG_gene_df['cell_type'] == cell_type]['gene'].unique()
            cur_non_CIG_genes = non_CIG_gene_df[non_CIG_gene_df['cell_type'] == cell_type]['gene'].unique()
            # print len(cur_CIG_genes), len(cur_non_CIG_genes)
            cur_table = pd.read_csv(all_dfs[marker][cell_type], sep='\t', index_col=0)

            for cig_gene in cur_CIG_genes:
                cur_transcripts = gene_to_transcript[cig_gene]
                cur_df = cur_table[cur_table.index.isin(cur_transcripts)]

                if cur_df.shape[0] == 0:
                    cur_cig_results.append([cig_gene, 0, 0, 0, 0, 0, 0, 0])
                else:
                    cur_cig_results.append([cig_gene,
                                        cur_df['total_width'].max(),
                                        cur_df['total_signal'].max(),
                                        cur_df['height'].max(),
                                        cur_df.ix[cur_df['skewness'].abs().argmax(), 'skewness'],
                                        cur_df.ix[cur_df['kurtosis'].abs().argmax(), 'kurtosis'],
                                        cur_df['single_signal'].max(),
                                        cur_df['single_width'].max()])

            for non_cig_gene in cur_non_CIG_genes:
                cur_transcripts = gene_to_transcript[non_cig_gene]
                cur_df = cur_table[cur_table.index.isin(cur_transcripts)]
                if cur_df.shape[0] == 0:
                    cur_non_cig_results.append([non_cig_gene, 0, 0, 0, 0, 0, 0, 0])
                else:
                    cur_non_cig_results.append([non_cig_gene,
                                            cur_df['total_width'].max(),
                                            cur_df['total_signal'].max(),
                                            cur_df['height'].max(),
                                            cur_df.ix[cur_df['skewness'].abs().argmax(), 'skewness'],
                                            cur_df.ix[cur_df['kurtosis'].abs().argmax(), 'kurtosis'],
                                            cur_df['single_signal'].max(),
                                            cur_df['single_width'].max()])
        cig_results.append(cur_cig_results)
        non_cig_results.append(cur_non_cig_results)

    cig_results = np.concatenate(cig_results, axis=1)
    non_cig_results = np.concatenate(non_cig_results,axis=1)


    columns = ['gene', 'total_width', 'total_signal', 'height', 'skewness', 'kurtosis', 'single_signal', 'single_width']

    columns = [marker + '_' + column for marker in markers for column in columns]
    print columns

    cig_result_df = pd.DataFrame(cig_results)
    cig_result_df.columns = columns

    non_cig_result_df = pd.DataFrame(non_cig_results)
    non_cig_result_df.columns = columns

    cig_result_df.to_csv('CIG_Great_Algo_stats.csv', index=False)
    non_cig_result_df.to_csv('non_CIG_Great_Algo_stats.csv', index=False)








# result = get_gene_transcript()
# print result['DDX11L1']


