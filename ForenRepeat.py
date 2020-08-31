

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
from ForenRepeatLib import myBasicCount
from ForenRepeatLib import myPatternReader
from ForenRepeatLib import myplot

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
        fig = myplot.func_save_plot(count_lst, str_name, path2)
        
        
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
    parser = argparse.ArgumentParser(description='Repeat quantification for reads from forenseq kit. Still working on', 
                                     usage='''python ForenRepeat.py <command> [<args>]

Available commands are:    
    quick    Repeat quantification by quick mode
    window   Repeat quantification by window scoring
    align    Repeat quantification by alignment scoring
    
    ''')
    subparsers = parser.add_subparsers(help='sub-command help')
    subparsers.required = True

    ## quick
    qucik_parser = subparsers.add_parser('quick', help='Repeat quantification by quick mode')
    qucik_parser.add_argument('--BAM', help='input1: *.bam file', default = None)
    qucik_parser.add_argument('--PAT', help='input2: repeat pattern file', default = None)
    qucik_parser.add_argument('--ID', help='output file name', default = None)
    qucik_parser.set_defaults(func=quick)

    ## window
    window_parser = subparsers.add_parser('window', help='Repeat quantification by window scoring. Developing')
    window_parser.add_argument('--BAM', help='input1: *.bam file', default = None)
    window_parser.add_argument('--PAT', help='input2: repeat pattern file', default = None)
    window_parser.add_argument('--ID', help='output file name', default = None)
        
    ## align
    align_parser = subparsers.add_parser('align', help='Repeat quantification by align scoring. Developing')
    align_parser.add_argument('--BAM', help='input1: *.bam file', default = None)
    align_parser.add_argument('--PAT', help='input2: repeat pattern file', default = None)
    align_parser.add_argument('--ID', help='output file name', default = None)

    ## run
    args = parser.parse_args()
    args.func(args)
