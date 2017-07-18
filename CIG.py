"""
This module is the main module for the CIG grid search
"""

from selector import random_control_genes
from grid_search import grid_search
from table import *
import pickle

def load_obj(name):
    f = open(name, 'rb')
    result = pickle.load(f)
    f.close()
    return result

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


    # cutoff_range = [12]
    # markers = ['h3k4me3']  #, 'h3k27ac', 'h3k27me3', 'h3k4me1']
    # markers = ['h3k4me3']
    # markers = ['h3k27ac']
    # markers = ['h3k27me3']
    markers = ['h3k4me1']
    # criterias = ['total_width', 'single_width', 'height']
    # criterias = ['total_signal', 'single_signal']
    # criterias = ['skewness', 'kurtosis']
    # criterias = ['total_width']
    # criterias = ['single_width']
    criterias = ['height']
    # criterias = ['total_signal']
    # criterias = ['single_signal']

    if 'skewness' in criterias or 'kurtosis' in criterias:
        option = True
        cutoff_range = range(5, 101, 5)
    else:
        option = False
        cutoff_range = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5]

    for marker in markers:
        # dfs_path = '/home/tmhbxx3/scratch/CIG/'+marker+'_peaks/pooled/'

        # this is for kurtosis and skewness
        if option:
            dfs_path = '/home/tmhbxx3/scratch/CIG/test/'+ marker + '_sk_peaks/'
        else:
            dfs_path = '/home/tmhbxx3/scratch/CIG/' + marker+ '_peaks/pooled/'

        dfs = os.listdir(dfs_path)
        all_dfs = defaultdict(dict)
        for table_name in dfs:
            info = table_name.split('_')
            # cell_type = info[1]
            # this is for kurtosis
            if option:
                cell_type = info[0]
            else:
                cell_type = info[1]
            cutoff = float(info[-1][:-4])
            # print cell_type, cutoff
            if cutoff in cutoff_range:
                # print 'start load table', table_name
                # all_dfs[cell_type][cutoff] = peak_table(dfs_path+table_name)
                if option:
                    all_dfs[cell_type][cutoff] = pd.read_csv(dfs_path + table_name)
                else:
                    all_dfs[cell_type][cutoff] = pd.read_csv(dfs_path+table_name, sep='\t')

        print len(all_dfs[cell_type].keys())
        # import pickle
        # with open('all_h3k4me3_peaks' + '.pkl', 'wb') as f:
        #     pickle.dump(all_dfs, f, pickle.HIGHEST_PROTOCOL)

        # print all_dfs
        ## step 3 do the grid and get the best CIG gene stats
        up_stream_distance_range = range(-250000, 250000, 1000)
        window_size_range = range(-250000, 250000, 1000)
        # window_size_range = [10000]

        for criteria in criterias:
            grid_path = grid_search(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                          up_stream_distance_range, window_size_range,
                                          all_dfs, cutoff_range, criteria, process=8)

            grid_path_df = pd.DataFrame(grid_path)

            grid_path_df.to_csv('grid_path_' + marker + '_' + criteria + '.csv')

            grid_path_results = []

            for path in grid_path:
                grid_path_results.append(path[0]+[path[1]])

            grid_path_results_df = pd.DataFrame(grid_path_results)

            grid_path_results_df.to_csv('grid_path_results'+marker+'_'+criteria+'.csv')
            # print grid_path
            # break
        # break
