

"""


    - ForenRepeat ver0.0.1
    - Written by Ren Zilin; zilin.ren@outlook.com
    - For reads from the forenseq kit and nanopore platform 


"""


## 
import os
import argparse
import pathlib

import pandas as pd

from tqdm import tqdm
## local pack
# import sys
# sys.path.insert(0, "myLibs")
from myLibs import myBasicCount
from myLibs.myLocalDPCount import func_reads_covering_str_locus, func_get_repeat_allele, func_str_genotyper
from myLibs import myPatternReader
# import myplot

## parameters
# parser = argparse.ArgumentParser(description='Repeat quantification for reads from forenseq kit. Still working on')
# parser.add_argument('--BAM', help='input1: *.bam file', default = None)
# parser.add_argument('--PAT', help='input2: repeat pattern file', default = None)
# parser.add_argument('--ID', help='output file name', default = None)
# args = parser.parse_args()

# bam_file_path = args.BAM
# pattern_file_path = args.PAT
# sample = args.ID

def quick(args):
    
    bam_file_path = args.BAM
    pattern_file_path = args.PAT
    sample = args.ID
    
    ## create dir
    path1  = 'nanoforenstr_output/%s' % sample
    path2  = 'nanoforenstr_output/%s/figs' % sample
    pathlib.Path(path2).mkdir(parents=True, exist_ok=True)
    
    
    ## pattern file
    pattern_dict = myPatternReader.func_read_in_pattern_file(pattern_file_path)
    
    
    result_lst = []
    for str_name in pattern_dict:
        chrom   = pattern_dict[str_name][0]
        start   = int(pattern_dict[str_name][1])
        end     = int(pattern_dict[str_name][2])
        pattern = pattern_dict[str_name][3]
    
        
        count_lst = myBasicCount.func_extract_reads_in_bam(bam_file_path, chrom, start, end, pattern)
        genotype, valid_alleles_num, total_cov, infos = myBasicCount.func_str_genotyper(count_lst)
        
        
        ## plot hist
        # fig = myplot.func_save_plot(count_lst, str_name, path2)
        
        
        ## output allele
        if valid_alleles_num < 2:
            result_lst.append([str_name, chrom, pattern, genotype[1], genotype[0], valid_alleles_num, total_cov, infos])
        else :
            result_lst.append([str_name, chrom, pattern, sorted(genotype)[0], sorted(genotype)[1], valid_alleles_num, total_cov, infos])


    pd.DataFrame(result_lst, 
                 columns=["# STR", "CHR", "pattern", 
                          "allele1", "allele2", "valid_alleles_num", 
                          "total_cov", "infos"]).to_csv('%s/forenRep.csv' % path1, index=None, na_rep="NA")
    return
    

######################  LA Local Align
def LA(args):
    
    def fetch_flanking_seq(chrom, start, end, flanking_len, reference_genome, samtools):
        pre_flanking_coor = "%s:%s-%s" % (chrom, int(start) - flanking_len, int(start)-1)
        suf_flanking_coor = "%s:%s-%s" % (chrom, int(end) + 1, int(end) + flanking_len)

        # print("## USE samtools to extract the flanking seq;\n## the command: %s faidx %s %s" % (samtools, reference_genome, pre_flanking_coor))
        faidx_output = os.popen("%s faidx %s %s" % (samtools, reference_genome, pre_flanking_coor)).read()
        pre_flanking = faidx_output.split('\n')[1]

        # print("## USE samtools to extract the flanking seq;\n## the command: %s faidx %s %s" % (samtools, reference_genome, suf_flanking_coor))
        faidx_output = os.popen("%s faidx %s %s" % (samtools, reference_genome, suf_flanking_coor)).read()
        sup_flanking = faidx_output.split('\n')[1]
        
        return pre_flanking, sup_flanking


    bam_file_path = args.BAM
    pattern_file_path = args.PAT
    sample = args.ID
    reference_genome = args.ref 
    samtools = args.samtools
    
    ## create dir
    
    path1  = 'nanoforenstr_output/%s' % sample
    pathlib.Path(path1).mkdir(parents=True, exist_ok=True)
    # path2  = 'nanoforenstr_output/%s/figs' % sample
    # pathlib.Path(path2).mkdir(parents=True, exist_ok=True)
    
    
    ## pattern file
    pattern_dict = myPatternReader.func_read_in_pattern_file(pattern_file_path)
    
    
    result_lst = []

    pbar = tqdm(pattern_dict)
    for str_ind, str_name in enumerate(pbar):
        pbar.set_description('## The locus %s is under inference {%d/%d}' % (str_name, str_ind + 1, len(pattern_dict)))



    # for str_ind, str_name in enumerate(pattern_dict):
        # print('## The locus %s is under inference {%d/%d}' % (str_name, str_ind + 1, len(pattern_dict)), flush = True)
        chrom   = pattern_dict[str_name][0]
        start   = int(pattern_dict[str_name][1])
        end     = int(pattern_dict[str_name][2])
        pattern = pattern_dict[str_name][3]
        
        
        ## reads
        flanking_len = 20
        seq_prefix, seq_suffix = fetch_flanking_seq(chrom, start, end, flanking_len, reference_genome, samtools)
        read_lst = func_reads_covering_str_locus(chrom, start, end, bam_file_path, flanking_len, seq_prefix, seq_suffix)

        ## copy number list
        copy_number_lst = func_get_repeat_allele(read_lst, pattern)
        
        ## Get genotype
        genotype, valid_alleles_num, total_cov, infos = func_str_genotyper(copy_number_lst, .56)
        
        ## plot hist
        # fig = myplot.func_save_plot(copy_number_lst, str_name, path2)
        
        ## output allele
        result_lst.append([str_name, chrom, pattern, genotype[0], genotype[1], valid_alleles_num, total_cov, infos])


    pd.DataFrame(result_lst, 
                 columns=["# STR", "CHR", "pattern", 
                          "allele1", "allele2", "valid_alleles_num", 
                          "total_cov", "infos"]).to_csv('%s/%s.csv' % (path1, sample), index=None, na_rep="NA")
    return
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Genotyping forensic STR', 
                                     usage='''python NanoForenRepeat.py <options> [<args>]

Available options are:    
    quick    Genotyping forensic STR by quick mode
    LA       Genotyping forensic STR by local align (recommend)
    
    ''')
    subparsers = parser.add_subparsers(help='sub-command help')
    subparsers.required = True

    ## quick
    qucik_parser = subparsers.add_parser('quick', help='Genotyping forensic STR by quick mode')
    qucik_parser.add_argument('--BAM', help='input1: *.bam file', default = None)
    qucik_parser.add_argument('--PAT', help='input2: repeat pattern file', default = None)
    qucik_parser.add_argument('--ID', help='output file name', default = None)
    qucik_parser.set_defaults(func=quick)

    ## localalign
    window_parser = subparsers.add_parser('LA', help='Genotyping forensic STR by local align')
    window_parser.add_argument('--BAM', help='input1: *.bam file', default = None)
    window_parser.add_argument('--PAT', help='input2: repeat pattern file', default = None)
    window_parser.add_argument('--ID', help='output file name', default = None)
    window_parser.add_argument('--ref', help='reference genome', default = None)
    window_parser.add_argument('--samtools', help='samtools', default = None)
    window_parser.set_defaults(func=LA)
        
    ## run
    args = parser.parse_args()
    args.func(args)
