from Wig import Wig
import os
import numpy as np

if __name__ == "__main__":
    # marker = 'h3k4me3'
    # marker = 'h3k27ac'
    marker = 'h3k4me1'
    # marker = 'h3k27me3'

    path = '/home/tmhbxx3/scratch/CIG/' + marker + '/'
    wigs = [x for x in os.listdir(path) if x.endswith('.wig')]
    # print wigs
    # for w in wigs:
    cur_wig = Wig(path+wigs[9])
    # cur_wig.CallPeak(0.25)
    for cutoff in np.arange(0.25, 5.25, 0.25):
        cur_wig.CallPeak(cutoff)
