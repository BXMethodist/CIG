from Wig import Wig
import os, sys
import numpy as np

if __name__ == "__main__":
    index = int(sys.argv[1])
    marker = sys.argv[2]

    # print index

    # marker = 'h3k4me3'
    # marker = 'h3k27ac'
    # marker = 'h3k4me1'
    # marker = 'h3k27me3'

    path = '/home/tmhbxx3/scratch/CIG/' + marker + '/'
    wigs = [x for x in os.listdir(path) if x.endswith('.wig')]

    cur_wig = Wig(path+wigs[index])
    for cutoff in range(1, 100, 1):
        cur_wig.CallPeak(cutoff)
