"""
This module is to calculate the ranges for different combinations of parameters
Such TSS, mid of TSS and TTS. relative size for gene body
"""
import pandas as pd

def get_range_absolute(gene_list, all_gene_GTF, left_distance, right_distance, TSS_pos='TSS', TTS_pos='TSS'):
    """

    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """
    results = set()
    for gene in gene_list:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'] == gene]
        # print cur_df
        for transcript in range(cur_df.shape[0]):
            cur_chr = cur_df.iloc[transcript, 1]
            cur_strand = cur_df.iloc[transcript, 2]
            cur_start = cur_df.iloc[transcript, 3]
            cur_end = cur_df.iloc[transcript, 4]

            if TSS_pos == 'TSS' and TTS_pos == 'TSS':
                if cur_strand == '+':
                    cur_left = cur_start + left_distance if cur_start + left_distance >=0 else 0
                    cur_right = cur_start + right_distance
                elif cur_strand == '-':
                    cur_right = cur_end - left_distance
                    cur_left = cur_end - right_distance if cur_end - right_distance >= 0 else 0
            if TSS_pos == 'TSS' and TTS_pos == 'TTS':
                if cur_strand == '+':
                    cur_left = cur_start + left_distance if cur_start + left_distance >=0 else 0
                    cur_right = cur_end + right_distance
                elif cur_strand == '-':
                    cur_right = cur_end - left_distance
                    cur_left = cur_start - right_distance if cur_start - right_distance >=0 else 0

            # print cur_left, cur_right, cur_start, cur_end, cur_strand
            # if cur_left > cur_right:
            #     continue
            results.add((gene, (cur_chr, cur_left, cur_right)))
    # print results
    results = list(results)
    result_df = pd.DataFrame(results)
    result_df.columns = ['gene', 'range']
    return result_df

def get_range_relative(gene_list, all_gene_GTF, left_relative, right_relative, TSS_pos='TSS', TTS_pos='TSS'):
    """

    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """
    results = set()
    for gene in gene_list:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'] == gene]
        # print cur_df
        for transcript in range(cur_df.shape[0]):
            cur_chr = cur_df.iloc[transcript, 1]
            cur_strand = cur_df.iloc[transcript, 2]
            cur_start = cur_df.iloc[transcript, 3]
            cur_end = cur_df.iloc[transcript, 4]

            cur_size = cur_end - cur_start

            left_distance = cur_size * left_relative
            right_distance = cur_size * right_relative

            if TSS_pos == 'TSS' and TTS_pos == 'TSS':
                if cur_strand == '+':
                    cur_left = cur_start + left_distance
                    cur_right = cur_start + right_distance
                elif cur_strand == '-':
                    cur_right = cur_end - left_distance
                    cur_left = cur_end - right_distance
            if TSS_pos == 'TSS' and TTS_pos == 'TTS':
                if cur_strand == '+':
                    cur_left = cur_start + left_distance
                    cur_right = cur_end + right_distance
                elif cur_strand == '-':
                    cur_right = cur_end - left_distance
                    cur_left = cur_start - right_distance

            # print cur_chr, cur_left, cur_right, cur_start, cur_strand
            if cur_left > cur_right:
                continue
            results.add((gene, (cur_chr, cur_left, cur_right)))
    # print results
    results = list(results)
    result_df = pd.DataFrame(results)
    result_df.columns = ['gene', 'range']
    return result_df

if __name__ == "__main__":
    all_gene_GTF = pd.read_csv('/archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.xls', sep='\t')
    all_gene_GTF['hg19.kgXref.geneSymbol'] = all_gene_GTF['hg19.kgXref.geneSymbol'].str.upper()

    all_genes = set(all_gene_GTF['hg19.kgXref.geneSymbol'].unique())

    exclude_list_file = open('/home/tmhbxx3/scratch/CIG/test/merge_ten_cells_relative_genes.txt', 'r')
    exclude_list = [x.strip() for x in exclude_list_file]
    exclude_list_file.close()
    exclude_list = set(exclude_list)

    CIG_gene_df = pd.read_csv('CIG_strong.csv')
    CIG_gene_df['gene'] = CIG_gene_df['gene'].str.upper()

    CIG_gene_list = list(CIG_gene_df['gene'].values)

    print get_range_absolute(all_gene_GTF, CIG_gene_list, 75000, 75000)