"""
This module is to calculate the ranges for different combinations of parameters
Such TSS, mid of TSS and TTS. relative size for gene body
"""
import pandas as pd

def get_range_absolute(gene_list, all_gene_GTF, left_distance, right_distance, TSS_pos, TTS_pos):
    """

    :param all_gene_GTF:
    :param gene_list:
    :param left_distance: upstream is -, downstream is +
    :param right_distance: upstream is -, downstream is +
    :param left_pos:
    :param right_pos:
    :return:
    """
    if gene_list is not None:
        cur_df = all_gene_GTF[all_gene_GTF['hg19.kgXref.geneSymbol'].isin(gene_list)]
    else:
        cur_df = all_gene_GTF

    positive_df = cur_df[cur_df['hg19.knownGene.strand'] == '+'].copy()
    negative_df = cur_df[cur_df['hg19.knownGene.strand'] == '-'].copy()

    if TSS_pos == 'TSS' and TTS_pos == 'TSS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df[positive_df['left_range'] < 0] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txStart'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txEnd'] - right_distance
        negative_df[negative_df['left_range'] < 0] = 0

    elif TSS_pos == 'TSS' and TTS_pos == 'TTS':
        positive_df['left_range'] = positive_df['hg19.knownGene.txStart'] + left_distance
        positive_df[positive_df['left_range'] < 0] = 0
        positive_df['right_range'] = positive_df['hg19.knownGene.txEnd'] + right_distance

        negative_df['right_range'] = negative_df['hg19.knownGene.txEnd'] - left_distance
        negative_df['left_range'] = negative_df['hg19.knownGene.txStart'] - right_distance
        negative_df[negative_df['left_range'] < 0] = 0

    elif TSS_pos == 'MID':
        positive_df['MID'] = (positive_df['hg19.knownGene.txStart'] + positive_df['hg19.knownGene.txEnd'])/2
        negative_df['MID'] = (negative_df['hg19.knownGene.txStart'] + negative_df['hg19.knownGene.txEnd'])/2

        positive_df['left_range'] = positive_df['MID'] + left_distance
        positive_df[positive_df['left_range'] < 0] = 0
        positive_df['right_range'] = positive_df['MID'] + right_distance

        negative_df['right_range'] = negative_df['MID'] - left_distance
        negative_df['left_range'] = negative_df['MID'] - right_distance
        negative_df[negative_df['left_range'] < 0] = 0

    new_df = positive_df.append(negative_df)
    # print new_df.columns
    result_df = new_df[['hg19.kgXref.geneSymbol', 'hg19.knownGene.chrom', 'left_range', 'right_range']]
    result_df.columns = ['gene', 'chr', 'left_range', 'right_range']

    return result_df

def get_range_relative(gene_list, all_gene_GTF, left_relative, right_relative, TSS_pos, TTS_pos):
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
    from time import time
    start = time()
    print get_range_absolute(all_genes, all_gene_GTF, -3000, 3000)
    end = time()
    print end - start