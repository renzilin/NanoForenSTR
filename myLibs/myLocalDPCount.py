import os
import re
import pathlib


import numpy as np
import pandas as pd

import pysam

import math

from ForenRepeatLib import myLocalDPAlign


def func_inference_on_str_locus(chrom, start, end, pattern, bam_file_path):
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    iters   = samfile.fetch(chrom, start, end, until_eof=False)
    copy_number_lst = []
    for rank, line in enumerate(iters):
        if line.seq == None:
            print("## the flag of no seq in bam file %s" % line.flag)
            continue
        else:
            read_seq = line.seq
            repeat_intervals, _ = get_interval(pattern, read_seq)
            break

            dist = [repeat_intervals[i][0] - repeat_intervals[i-1][1] for i in range(1, len(repeat_intervals))]
            dists.append(dist)
    return
    

def func_read_in_pattern_file(pattern_file_path):
    pattern_dict = {}
    with open(pattern_file_path, 'r') as file:
        for line in file:
            line_lst = line.strip().split(",")
            pattern_dict[line_lst[0]] = [line_lst[1], line_lst[2], line_lst[3], line_lst[4]]
    return pattern_dict



def func_reads_covering_str_locus(chrom, start, end, pattern, bam_file_path):
    
    reads_lst = []
    
    samfile = pysam.AlignmentFile(bam_file_path, 'rb')
    iters   = samfile.fetch(chrom, start, end, until_eof=False)
    for rank, line in enumerate(iters):
        if line.seq == None:
            print("## the flag with no seq in bam file %s" % line.flag)
            continue
        
        else:
            read_seq = line.seq
            reads_lst.append(read_seq)    
    return reads_lst


def get_interval(pattern, read_seq):

    intervals_lst = []
    repeat_lst = []

    for inds, repeat_inds in enumerate(re.finditer(pattern, read_seq)):
        start, end = repeat_inds.span()
        if inds == 0:
            repeat_lst.append([read_seq[start:end]])
            intervals_lst.append([start, end])

        elif start == end_previous:
            repeat_lst[-1].append(read_seq[start:end])
            intervals_lst[-1][1] = end

        else:
            repeat_lst.append([read_seq[start:end]])
            intervals_lst.append([start, end])

        start_previous = start
        end_previous   = end
    
    return intervals_lst, repeat_lst


def func_get_skip_intervals(intervals_lst):
    
    skip_intervals = []

    for ind, interval in enumerate(intervals_lst):
        if ind == 0:
            skip_intervals.append([interval[1]])

        elif ind + 1 == len(intervals_lst):
            skip_intervals[-1].append(interval[0])

        else:
            skip_intervals[-1].append(interval[0])
            skip_intervals.append([interval[1]])
    
    return skip_intervals



def func_score_on_unit(ref_unit, seq_unit):
    mismatch_num = 0
    gap_num = 0
    for ind, ref_base in enumerate(ref_unit):
        seq_base = seq_unit[ind]
        
        if seq_base == '-':
            gap_num += 1
        elif seq_base != ref_base:
            mismatch_num += 1
        
    return mismatch_num, gap_num



def func_if_aligned_repeat(aligned_ref, aligned_seq, pattern):
    
    if '-' in aligned_ref:
#         print(aligned_ref, '\n', aligned_seq)
        return False
    
    if len(aligned_ref) == len(pattern):
        mismatch_num, gap_num = func_score_on_unit(aligned_ref, aligned_seq)
        if gap_num <=2 and mismatch_num == 0:
            return True
        elif gap_num == 0 and mismatch_num <= 2:  # 2 or 1?
            return True
        else:
            return False
    
    ref_unit_lst = re.findall(r'.{%s}' % len(pattern), aligned_ref)
    seq_unit_lst = re.findall(r'.{%s}' % len(pattern), aligned_seq)
    
    
    for ind, ref_unit in enumerate(ref_unit_lst):
        seq_unit = seq_unit_lst[ind]
        mismatch_num, gap_num = func_score_on_unit(ref_unit, seq_unit)
        
        if mismatch_num > 1 or gap_num > 1:
            return False
            
    return True
    



def func_repeat_interval(intervals_lst, skip_intervals_lst, read_seq, pattern):
    if_connect = []
    for ind, skip_interval_coor in enumerate(skip_intervals_lst):

        skip_start = skip_interval_coor[0]
        skip_end   = skip_interval_coor[1]

        skip_seq   = read_seq[skip_start : skip_end]
        local_ref  = pattern * math.ceil( len(skip_seq) / len(pattern) )

        align_score, aligned_ref, aligned_seq = myLocalDPAlign.align(local_ref, skip_seq, 0, -5, -5)

        if func_if_aligned_repeat(aligned_ref, aligned_seq, pattern):
            if_connect.append(1)
        else:
            if_connect.append(0)
    
    
    updated_intervals_lst = []
    for ind in range(len(intervals_lst)):
        if ind == 0:
            updated_intervals_lst.append([intervals_lst[ind]])
        elif if_connect[ind - 1]:
            updated_intervals_lst[-1].append(intervals_lst[ind])
        else:
            updated_intervals_lst.append([intervals_lst[ind]])
    
    
    repeat_interval = []
    for update_intervals in updated_intervals_lst:
        if len(update_intervals) == 1:
            repeat_interval.append(update_intervals[0])
        else:
            repeat_interval.append([update_intervals[0][0], update_intervals[-1][-1]])
    
    repeat_interval_len = [i[1] - i[0] for i in repeat_interval]
    
    max_interval = repeat_interval[repeat_interval_len.index(max(repeat_interval_len))]
    max_interval_len = max_interval[1] - max_interval[0]
    
    return repeat_interval, repeat_interval_len, max_interval_len

def func_get_repeat_allele(read_lst, pattern):
    copy_number_lst = []
    
    if len(read_lst) <= 10:   ## coverage
        return copy_number_lst
    
    for read_ind, read_seq in enumerate(read_lst):
        
        intervals_lst, repeat_lst = get_interval(pattern, read_seq)
        if len(intervals_lst) > 1:
            skip_intervals_lst = func_get_skip_intervals(intervals_lst)
            repeat_interval, repeat_interval_len, max_interval_len = func_repeat_interval(intervals_lst, skip_intervals_lst, read_seq, pattern)
        
        elif len(intervals_lst) == 1:
            repeat_interval, repeat_interval_len, max_interval_len = intervals_lst, 1, (intervals_lst[0][1]-intervals_lst[0][0])
        
        else:
            repeat_interval, repeat_interval_len, max_interval_len = 0, 0, 0
        
        copy_number = max_interval_len//len(pattern)
        
        copy_number_lst.append(copy_number)
    return copy_number_lst

def func_str_genotyper(count_lst, cutoff = .5, coverage = 10):
    
    if len(count_lst) == 0:
        return [0, 0], -1, 0, -1
    
    elif len(count_lst) <= coverage and len(count_lst) > 0:
        infos = {}
        for i in count_lst:
            if i not in infos:
                infos[i] = 1
            else:
                infos[i] += 1
        
        return [0, 0], -1, 0, str(infos)


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