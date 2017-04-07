#!/usr/bin/env python3

import pysam
import numpy as np
from optparse import OptionParser
from multiprocessing import Pool
import sys
import copy

CHECK_VERSION = "0.1.3"
AUTHOR = "youcai"

def parseCommand():
    usage = "\n\tfor Bam QC"
    parser = OptionParser(usage = usage, version = CHECK_VERSION)
    parser.add_option("-b", "--bam", dest = "bam_file",
        help = "file name of bam file, must be sorted and mark duplication, *[requried]")
    parser.add_option("-r", "--bed", dest = "bed_file",
        help = "file name of capture bed file, must be treat with `bedtools merge`, *[requried]")
    parser.add_option("-o", "--output", dest = "output",
        help = "output file name, <stdout>")
    return parser.parse_args()

def Target_region(bed_file):
    '''
    treat bed file
    '''
    with open(bed_file) as f:
        line_1st = f.readline().strip().split('\t')
        target_area = np.array(line_1st)
        for line in f:
            text = line.strip().split('\t')
            tmp_area = np.array(text)
            target_area = np.vstack((target_area, tmp_area))
    return target_area

def Chr_core(samfile, small_bed):
    '''
    for every chromosome, do this analysis
    '''
    chr_coverage = []
    chr_alignments = 0
    for line in small_bed:
        ACGT = samfile.count_coverage(str(line[0]), int(line[1]), int(line[2]))
        chr_coverage.append(sum(np.array(ACGT)))
        chr_alignments += samfile.count(str(line[0]), int(line[1]), int(line[2]))
    return chr_coverage, chr_alignments

def BedBamQC(samfile, target_region, depth_cut = 100):
    '''
    get every samll region's capture_alignments,
    and coverage for every position, then calculate the results
    '''
    starts = list(map(int, target_region[:, 1]))
    ends = list(map(int, target_region[:, 2]))
    target_size = np.sum(np.array(ends) - np.array(starts))
    all_alignments = samfile.mapped + samfile.unmapped
    # 1. parameters
    input_list = []
    chromosome = set(target_region[:, 0])
    for chro in chromosome:
        mask = target_region[:, 0] == chro
        chr_bed = target_region[mask]
        sam = copy.copy(samfile) ###
        tmp_argv = (sam, chr_bed)
        input_list.append(tmp_argv)
    # 2. assignments
    P = Pool(processes = len(chromosome))
    Process = P.starmap(Chr_core, input_list)
    # 3. results collection
    coverage_arr = [] 
    capture_alignments = 0
    for res in Process:
        coverage_arr += res[0] # [] + []
        capture_alignments += res[1] # int + int
    ###for one_region in target_region:
    ###    ACGT = samfile.count_coverage(str(one_region[0]), int(one_region[1]), int(one_region[2]))
    ###    coverage_arr.append(sum(np.array(ACGT)))
    ###    # if one alignment in two or more regions
    ###    # so said
    ###    capture_alignments += samfile.count(str(one_region[0]), int(one_region[1]), int(one_region[2]))
    capture_rate = float(capture_alignments) / all_alignments * 100 ##
    target_base = np.sum(list(map(np.sum, coverage_arr)))
    depth_in_target = float(target_base) / target_size ##
    capture_size = sum(map(lambda x : len(x[x > 0]), coverage_arr))
    target_coverage = float(capture_size) / target_size * 100 ##
    good_capture_size = sum(map(lambda x : len(x[x > depth_cut]), coverage_arr))
    target_100X = float(good_capture_size) / target_size *100 ##
    return capture_rate, depth_in_target, target_coverage, target_100X

def BamQC(samfile, target_area):
    '''
    use bam.bai to calculate mapping_rate, the loop bam file to get mapping_quality,
    duplication_rate, median insert_size
    '''
    mapped_alignments = samfile.mapped
    unmapped_alignments = samfile.unmapped
    all_alignments = mapped_alignments + unmapped_alignments
    mapping_quality = 0
    duplicate_alignments = 0
    insert_arr = []
    for alignment in samfile:
        if alignment.is_unmapped:
            continue
        if alignment.is_duplicate:
            duplicate_alignments += 1
        mapping_quality += alignment.mapping_quality
        insert_arr.append(alignment.template_length)
    # result values
    mapping_rate = float(mapped_alignments) / all_alignments * 100 ##
    mapping_quality = float(mapping_quality) / all_alignments ##
    duplication_rate = float(duplicate_alignments) / all_alignments * 100 ##
    insert_arr = np.array(insert_arr)
    insert_arr = insert_arr[insert_arr>0]
    insert_size = int(np.median(insert_arr)) ##
    return mapping_rate, mapping_quality, insert_size, duplication_rate

def main():
    (options, args) = parseCommand()
    samfile = pysam.AlignmentFile(options.bam_file)
    target_area = Target_region(options.bed_file)

    mapping_rate, mapping_quality, insert_size, duplication_rate = BamQC(samfile, target_area)    
    capture_rate, depth_in_target, target_coverage, target_100X = BedBamQC(samfile, target_area, 100)

    f = sys.stdout if options.output == None else open(options.output, 'w')
    print("Mapping rate: %.3f%%" %mapping_rate, file = f)
    print("Mapping quality: %.3f" %mapping_quality, file = f)
    print("Insert size: %dbp" %insert_size, file = f)
    print("Duplication rate: %.3f%%" %duplication_rate, file = f)
    print("Capture rate: %.3f%%"  %capture_rate, file = f)
    print("Depth in target: %.3fX" %depth_in_target, file = f)
    print("Target coverage: %.3f%%" %target_coverage, file = f)
    print("Target>100X: %.3f%%" %target_100X, file = f)

    samfile.close()
    f.close()

if __name__ == "__main__":
    main()
