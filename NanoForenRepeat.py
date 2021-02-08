

"""


    - ForenRepeat ver0.1
    - Written by Ren Zilin; zilin.ren@outlook.com
    - For reads from the forenseq kit and nanopore platform 


"""


## 
import argparse
import pathlib

import pandas as pd

## local pack
# import sys
# sys.path.insert(0, "myLibs")
from myLibs import myBasicCount
from myLibs import myLocalDPCount
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
    path1  = 'forenRepeat/%s' % sample
    path2  = 'forenRepeat/%s/figs' % sample
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
    
    bam_file_path = args.BAM
    pattern_file_path = args.PAT
    sample = args.ID
    
    ## create dir
    pathlib.Path('forenRepeat_output_LA/').mkdir(parents=True, exist_ok=True)
    path1  = 'forenRepeat_output_LA/%s' % sample
    path2  = 'forenRepeat_output_LA/%s/figs' % sample
    pathlib.Path(path2).mkdir(parents=True, exist_ok=True)
    
    
    ## pattern file
    pattern_dict = myPatternReader.func_read_in_pattern_file(pattern_file_path)
    
    
    result_lst = []
    for str_ind, str_name in enumerate(pattern_dict):
        print('## The locus %s is under inference {%d/%d}' % (str_name, str_ind + 1, len(pattern_dict)), flush = True)
        chrom   = pattern_dict[str_name][0]
        start   = int(pattern_dict[str_name][1])
        end     = int(pattern_dict[str_name][2])
        pattern = pattern_dict[str_name][3]

        ## reads
        flanking_len = 30
        read_lst = myLocalDPCount.func_reads_covering_str_locus(chrom, start, end, bam_file_path, flanking_len)
        
        ## copy number list
        copy_number_lst = myLocalDPCount.func_get_repeat_allele(read_lst, pattern)
        
        ## Get genotype
        genotype, valid_alleles_num, total_cov, infos = myLocalDPCount.func_str_genotyper(copy_number_lst)
        
        ## plot hist
        # fig = myplot.func_save_plot(copy_number_lst, str_name, path2)
        
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
    window_parser.set_defaults(func=LA)
        
    ## run
    args = parser.parse_args()
    args.func(args)
