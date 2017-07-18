"""
This module is the main module for the CIG scipy minimizer
"""
from selector import random_control_genes
from table import *
from CIG_predict import optimize_allocs

if __name__ == "__main__":
    ### step 1 get random negative control set
    all_gene_GTF = pd.read_csv('/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()

    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())

    exclude_list_file = open('/home/tmhbxx3/scratch/CIG/test/merge_ten_cells_relative_genes.txt', 'r')
    exclude_list = [x.strip() for x in exclude_list_file]
    exclude_list_file.close()
    exclude_list = set(exclude_list)

    CIG_gene_df = pd.read_csv('CIG_strong.csv')
    CIG_gene_df['gene'] = CIG_gene_df['gene'].str.upper()

    # print CIG_gene_df.shape, 'CIG genes'

    non_CIG_gene_df = random_control_genes(CIG_gene_df, exclude_list, all_genes, random_seed=9000)

    # print non_CIG_gene_df.shape, 'non CIG gens'

    ### step 2 get all_dfs

    # cutoff_range = range(1, 301)
    cutoff_range = [10]

    marker = 'h3k4me3'
    dfs_path = '/home/tmhbxx3/scratch/CIG/h3k4me3_peaks/pooled/'
    dfs = os.listdir(dfs_path)
    all_dfs = defaultdict(dict)
    for table_name in dfs:
        info = table_name.split('_')
        cell_type = info[1]
        cutoff = int(info[-1][:-4])
        # print cell_type, cutoff
        if cutoff in cutoff_range:
            print 'start load table', table_name
            # all_dfs[cell_type][cutoff] = peak_table(dfs_path+table_name)
            all_dfs[cell_type][cutoff] = pd.read_csv(dfs_path+table_name, sep='\t')

    allocs = optimize_allocs(CIG_gene_df, non_CIG_gene_df, all_gene_GTF, all_dfs, 'total_width')
    print allocs

