import pandas as pd, os
from bisect import bisect_left, bisect_right
from collections import defaultdict

class peak_table:
    def __init__(self, table_path):
        f = open(table_path, 'r')
        info = f.readlines()[1:]
        f.close()
        self.peaks = {}
        self.peaks_index = defaultdict(list)

        for line in info:
            cur_line = [x.strip() for x in line.split('\t')]
            chromosome, start, end = cur_line[0:3]
            width, total_signal, height = cur_line[4:7]
            cur_peak = peak(chromosome, int(start), int(end), int(width), float(total_signal), float(height))
            if chromosome in self.peaks:
                self.peaks[chromosome][start].add(cur_peak)
                self.peaks[chromosome][end].add(cur_peak)
            else:
                self.peaks[chromosome] = defaultdict(set)
                self.peaks[chromosome][start].add(cur_peak)
                self.peaks[chromosome][end].add(cur_peak)
        for chromosome in self.peaks.keys():
            self.peaks_index[chromosome] = sorted(self.peaks[chromosome].keys())

    def get_peaks(self, chromosome, start, end):
        indexes = self.peaks_index[chromosome]
        start_index = bisect_left(indexes, start)
        end_index = bisect_right(indexes, end)
        candidates = indexes[start_index:end_index]
        result_peaks = set()
        for c in candidates:
            for p in self.peaks[chromosome][c]:
                result_peaks.add(p)

        return result_peaks


class peak:
    def __init__(self, chromosome, start, end, width, total_signal, height):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.width = width
        self.total_signal = total_signal
        self.height =height
