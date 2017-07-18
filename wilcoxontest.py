from scipy import stats
import pandas as pd, numpy as np

def logP_wilcoxon(groupA, groupB):
    """
    return the negative log P value for two groups
    :param groupA:
    :param groupB:
    :return:
    """
    try:
        rank_diff, p = stats.mannwhitneyu(groupA, groupB, alternative='less')
        # print p
        if np.log10(p) < -80:
            return 0
        return np.log10(p)
    except:
        return 0

