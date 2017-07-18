from scipy import stats
import pandas as pd, numpy as np

def logP_wilcoxon(groupA, groupB, bags=20):
    """
    return the negative log P value for two groups
    :param groupA:
    :param groupB:
    :return:
    """
    p_values = []

    for i in range(bags):
        cur_groupA = np.random.choice(groupA, int(len(groupA)*0.8), replace=True)
        cur_groupB = np.random.choice(groupB, int(len(groupB)*0.8), replace=True)
        try:
            rank_diff, p = stats.mannwhitneyu(cur_groupA, cur_groupB, alternative='less')
            # print p
            # There is a problem in multiple processing. for some reason, the function will return a very small number
            # However, it could pass the test without problem if just run the corresponding dataframe
            if np.log10(p) < -80:
                p_values.append(0)

            p_values.append(np.log10(p))
        except:
            p_values.append(0)

    return np.mean(p_values)