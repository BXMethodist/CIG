from scipy import stats
import numpy as np

def peak_skewness(peak_signals):
    """
    :param peak_signals: has two rows, first rows is genome index, 2nd is the signals
    :return:
    """
    values = []
    for i in range(peak_signals.shape[1]):
        index = peak_signals[0, i]
        count = int(peak_signals[1, i] * 10)
        values += [index] * count

    return stats.skew(np.asarray(values))