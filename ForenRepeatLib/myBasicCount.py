"""

    - This is based on the quantification used in CE result
    - Zilin Ren @Author
    - 20200831  @Data
    - zilin.ren@outlook.com @Email

"""


import re

import pysam

import numpy as np
import pandas as pd

def func_copy_number_by_regex(read_seq_str, pattern, dist_cutoff=8):
    repeat_intervals = []
    repeat_lst = []
    start_previous = -1
    end_previous   = -1
    
    for inds, repeat_inds in enumerate(re.finditer(pattern, read_seq_str)):
        start, end = repeat_inds.span()

        if inds == 0:
            repeat_lst.append([read_seq_str[start:end]])
            repeat_intervals.append([start, end])

        elif start == end_previous:
            repeat_lst[-1].append(read_seq_str[start:end])
            repeat_intervals[-1][1] = end

        elif start - end_previous <= dist_cutoff and inds != 0:
            repeat_lst[-1].append(read_seq_str[end_previous:start])
            repeat_lst[-1].append(read_seq_str[start:end])
            repeat_intervals[-1][1] = end

        else:
            repeat_lst.append([read_seq_str[start:end]])
            repeat_intervals.append([start, end])
        
        start_previous = start
        end_previous   = end

    if len(repeat_intervals) == 0:
            return -1, -1

    else:
        
        copy_nums = [len(i) for i in repeat_lst]
        unit_count_lst = []
        unit_pattern_lst = []
        for unit in repeat_lst[copy_nums.index(max(copy_nums))]:
            if unit_count_lst == []:
                unit_pattern_lst.append(unit)
                unit_count_lst.append(1)
            else:
                if unit_pattern_lst[-1] == unit:
                    unit_count_lst[-1] += 1
                else:
                    unit_pattern_lst.append(unit)
                    unit_count_lst.append(1)
        ## output pattern
        seq_pattern = ''
        for ind, unit in enumerate(unit_pattern_lst):
            seq_pattern += '[%s]%s' % (unit, unit_count_lst[ind])
        
        return max(copy_nums), seq_pattern
    
def func_extract_reads_in_bam(bam_file_path, chrom, start, end, pattern):
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    iters   = samfile.fetch(chrom, start, end, until_eof=True)
    copy_number_lst = []
    for rank, line in enumerate(iters):
        if line.seq == None:
            # print("## the flag of no seq in bam file: %s" % line.flag)
            continue
        else:
            read_seq = line.seq
        
        copy_number, seq_pattern = func_copy_number_by_regex(read_seq, pattern)
        copy_number_lst.append(copy_number)
    return copy_number_lst

def func_str_genotyper(count_lst, cutoff = .5):
    
    infos = {}
    for i in count_lst:
        if i not in infos:
            infos[i] = 1
        else:
            infos[i] += 1
    
    total_cov = len(count_lst)
    count_reads_lst  = []
    count_number_lst = []
    
    for i in count_lst:
        if i not in count_number_lst:
            count_number_lst.append(i)
            count_reads_lst.append(1)
        else:
            count_reads_lst[count_number_lst.index(i)] += 1
    
    count_df = pd.DataFrame(list( zip( *[count_number_lst, count_reads_lst] ) ), columns=['copy_number', 'reads_support'])
    count_df['total_cov'] = total_cov
    count_df['cov_ratio'] = count_df['reads_support'] / count_df['total_cov']
    count_df['max_ratio'] = count_df['cov_ratio'] / max(count_df['cov_ratio'])
    count_left_df = count_df.loc[count_df['max_ratio'] >= cutoff, :]
    
    if count_left_df.shape[0] == 1:
        return count_left_df["copy_number"].tolist() * 2, 2.1, total_cov, str(infos)
    
    elif count_left_df.shape[0] == 2:
        return count_left_df["copy_number"].tolist(), 2.2, total_cov, str(infos)

    elif count_left_df.shape[0] == 3:
        temp_df = count_left_df.sort_values(by = 'max_ratio', ascending=False)
        return temp_df["copy_number"].tolist()[:2], 2.3, total_cov, str(infos)
    else:
        temp_df = count_left_df.sort_values(by = 'max_ratio', ascending=False)
        return temp_df["copy_number"].tolist()[:1] + [pd.NA], 1, total_cov, str(infos)