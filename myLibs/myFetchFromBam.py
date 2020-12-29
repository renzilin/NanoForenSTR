import os
import sys

def extract_reads_from_bam(bamfile, coordinate, samtools):
    print("## USE samtools to extract the reads covering the REGION;\n"+
    "## the command: %s view -F 2304 %s %s" % (samtools, bamfile, coordinate))
    
    file = os.popen('%s view -F 2304 %s %s' % (samtools, bamfile, coordinate)).readlines()

    reads = []
    for line in file:
        reads.append(line.strip().split())

    if reads == []:
        print("No reads in samfile")
        sys.exit(1)

    return reads

if __name__ == '__main__':
    print(1)