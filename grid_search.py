"""
This module will implement the grid_search
"""

import numpy as np
from copy import deepcopy
from multiprocessing import Process, Queue

def next_grid(all_range, range_step, cur_step, cur_center, reduction_step, step_limit, search):
    new_step = cur_step/reduction_step
    if new_step < step_limit:
        new_step = step_limit
    center_index = all_range.index(cur_center)

    # print center_index
    if search:
        if range_step == cur_step:
            start_index = center_index - 2 if center_index - 2 >= 0 else 0
            end_index = center_index + 3 if center_index + 3 <= len(all_range) else len(all_range)
            # print start_index, end_index
            return all_range[start_index:end_index], new_step
        else:
            start_index = center_index - int(cur_step*2/range_step) if center_index - int(cur_step*2/range_step) >=0 else 0
            end_index = center_index + int(cur_step*2/range_step) + 1 if center_index + int(cur_step*2/range_step) + 1 <= len(all_range) else len(all_range)
            # print start_index, end_index

            return all_range[start_index:end_index][::int(new_step/range_step)], new_step
    else:
        start_index = center_index - 4 if center_index - 4 >= 0 else 0
        end_index = center_index + 5 if center_index + 5 <= len(all_range) else len(all_range)
        return all_range[start_index:end_index][::int(new_step / range_step)], new_step

def grid_search(CIG_gene_df, non_CIG_gene_df,
                all_gene_GTF,up_stream_distance_range, window_size_range,
                all_dfs, cutoff_range, criteria, cost_function,
                TSS_pos, TTS_pos,
                up_stream_distance_grid=50000, window_size_grid=50000, cutoff_grid=4,
                up_stream_distance_range_step=1000, window_size_range_step=1000, cutoff_range_step=1,
                up_stream_distance_step=2, window_size_step=2, cutoff_step=2,
                up_stream_distance_limit=1000, window_size_limit=1000, cutoff_limit=1,
                process=8, wigs=None):

    ## track the parameters and logP path
    path = []

    new_up_stream_distance_range = deepcopy(up_stream_distance_range)
    new_window_size_range = deepcopy(window_size_range)
    new_cutoff_range = deepcopy(cutoff_range)

    search = True

    # iter = 1
    # count = 0

    grid_up_stream_distance_range = new_up_stream_distance_range[
                                    up_stream_distance_grid / up_stream_distance_range_step / 2::(
                                    up_stream_distance_grid / up_stream_distance_range_step)]
    grid_window_size_range = new_window_size_range[
                             window_size_grid / window_size_range_step / 2::(window_size_grid / window_size_range_step)]
    grid_cutoff_range = new_cutoff_range[int(cutoff_grid / cutoff_range_step / 2)::int(cutoff_grid / cutoff_range_step)]

    # grid_up_stream_distance_range = [4000]
    # grid_window_size_range = range(1000, 500000)
    # grid_cutoff_range = [10]

    past_path = {}

    while True:
        # determine the new grid

        print grid_up_stream_distance_range
        print grid_window_size_range
        print grid_cutoff_range

        # if count >= iter:
        #     break

        combinations = np.stack(np.meshgrid(grid_up_stream_distance_range, grid_window_size_range, grid_cutoff_range), -1).reshape(-1, 3)

        # print combinations
        if TSS_pos == TTS_pos:
            new_combinations = []
            for comb in combinations:
                if comb[1] > comb[0]:
                    new_combinations.append(comb)
            combinations = np.asarray(new_combinations)
        print combinations.shape, 'current grid size is '

        # combinations = np.asarray([[-3000, 1000000, 2]])

        # print combinations
        #
        # return

        best_log_P = None
        best_comb = None

        chunks = []
        cur_index = 0
        reminder = len(combinations) % process
        chunk_size = len(combinations) / process
        for i in range(process):
            if reminder > 0:
                chunks.append(combinations[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
                cur_index += 1
                reminder -= 1
            else:
                chunks.append(combinations[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

        total_chunk_size = 0
        for chunk in chunks:
            total_chunk_size += len(chunk)
        if total_chunk_size != len(combinations):
            print 'multiple processes chunk size is not correct'
            return None

        queue = Queue()
        processes = []

        cur_past_path = deepcopy(past_path)
        for i in range(process):
            cur_chunk = chunks[i]
            p = Process(target=CIG_process,
                        args=(queue, cur_chunk, CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                              all_dfs, criteria, cur_past_path, cost_function, TSS_pos, TTS_pos,
                              wigs))
            processes.append(p)
            p.start()

        cur_path = []

        for i in range(process):
            cur_best_logP, cur_best_comb, cur_process_path = queue.get()
            cur_path += cur_process_path
            # print cur_best_logP, cur_best_logP< best_log_P, best_log_P is None or cur_best_logP < best_log_P
            if (best_log_P is None or cur_best_logP < best_log_P) and cur_best_logP is not None:
                best_log_P = cur_best_logP
                best_comb = cur_best_comb
            print best_comb, best_log_P
        print best_comb, best_log_P
        for p in processes:
            p.join()

        for trace in cur_path:
            if tuple(trace[0]) not in past_path:
                past_path[tuple(trace[0])] = trace[1]

        path += cur_path

        if not search:
            break

        # break
        if up_stream_distance_grid == up_stream_distance_limit and window_size_grid == window_size_limit and cutoff_grid == cutoff_limit:
            search = False
            # break
        # use best comb to determine new range

        up_stream_center, window_size_center, cutoff_center = best_comb
        grid_up_stream_distance_range, up_stream_distance_grid = next_grid(up_stream_distance_range,
                                                                           up_stream_distance_limit,
                                                                           up_stream_distance_grid,
                                                                           up_stream_center,
                                                                           up_stream_distance_step,
                                                                           up_stream_distance_limit,
                                                                           search)
        grid_window_size_range, window_size_grid = next_grid(window_size_range,
                                           window_size_limit,
                                           window_size_grid,
                                           window_size_center,
                                           window_size_step,
                                           window_size_limit,
                                           search)
        grid_cutoff_range, cutoff_grid = next_grid(cutoff_range,
                                      cutoff_limit,
                                      cutoff_grid,
                                      cutoff_center,
                                      cutoff_step,
                                      cutoff_limit,
                                      search)

        # break

    return path


def CIG_process(queue, combinations, CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                all_dfs, criteria, cur_past_path, cost_function, TSS_pos, TTS_pos,
                wigs):
    path = []
    best_log_P = None
    best_comb = None

    for comb in combinations:
        # print comb[0], comb[1], comb[1] <= comb[0]

        cur_up_stream_distance, cur_window_size, cur_cutoff = comb
        if (cur_up_stream_distance, cur_window_size, cur_cutoff) in cur_past_path:
            cur_logP = cur_past_path[(cur_up_stream_distance, cur_window_size, cur_cutoff)]
        else:
            cur_logP = cost_function(CIG_gene_df, non_CIG_gene_df, all_gene_GTF,
                                     cur_up_stream_distance, cur_window_size,
                                     all_dfs, cur_cutoff, criteria, TSS_pos, TTS_pos,
                                     wigs)
        print comb, cur_logP
        path.append((comb, cur_logP))
        if (best_log_P is None or best_log_P > cur_logP) and cur_logP is not None:
            best_log_P = cur_logP
            best_comb = comb
    # print best_comb, best_log_P, 'current process best value is'
    #     break
    # print 'for current process', best_log_P, best_comb
    queue.put((best_log_P, best_comb, path))
    return