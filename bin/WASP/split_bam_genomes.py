# ts220221 - split_bam_genomes.py
#
# Mostly copied from the find_intersecting_snps.py function from WASP, after
# I accidentally deleted my old split_bam script. Oops. 
# 
# All credits go to WASP: https://github.com/bmvdgeijn/WASP

import sys
import os
import gzip
import argparse
import numpy as np
from itertools import product, groupby

import pysam

import util
import snptable

import tables

MAX_SEQS_DEFAULT = 64
MAX_SNPS_DEFAULT = 6

#######################################
### Option parser

def parse_options():
    parser = argparse.ArgumentParser(description="This program will separate "
                                                 "mapped reads (that pass the "
                                                 "WASP criteria) based on the "
                                                 "genomes that are present in "
                                                 "data.")
    
    parser.add_argument("snp_dir", help="directory with SNPs")
    parser.add_argument("bam", help="input BAM with filtered reads")
    parser.add_argument("genomes", help="comma-separated genomes present")
    parser.add_argument("output_dir", help="output directory")
    parser.add_argument("read_stats_file", help="read stats file name")

    return parser.parse_args()
    

#######################################
### Read stats class

class ReadStats(object):
    """Track information about reads and SNPs that they overlap"""

    def __init__(self):
        
        self.total = 0
        
        # number of read matches to reference allele
        self.ref_count = 0
        # number of read matches to alternative allele
        self.alt_count = 0
        # number of reads that overlap SNP but match neither allele
        self.other_count = 0

        # number of reads discarded becaused not mapped
        self.discard_unmapped = 0
        
        ## number of reads discarded because not proper pair
        #self.discard_improper_pair = 0

        ## number of reads discarded because mate unmapped
        #self.discard_mate_unmapped = 0

        ## paired reads map to different chromosomes
        #self.discard_different_chromosome = 0
        
        # no overlapping SNP
        self.nosnp = 0

        # number of reads discarded because overlap an indel
        self.discard_indel = 0

        # number of reads discarded because secondary match
        self.discard_secondary = 0

        # number of chimeric reads discarded
        self.discard_supplementary = 0

        ## number of reads discarded because of too many overlapping SNPs
        #self.discard_excess_snps = 0
        
        ## number of reads discarded because too many allelic combinations
        #self.discard_excess_reads = 0

        ## when read pairs share SNP locations but have different alleles there
        #self.discard_discordant_shared_snp = 0
        
        ## reads where we expected to see other pair, but it was missing
        ## possibly due to read-pairs with different names
        #self.discard_missing_pair = 0
        
        # number of single reads kept
        self.keep_single = 0
        ## number of read pairs kept
        #self.keep_pair = 0

        ## number of single reads that need remapping
        #self.remap_single = 0
        ## number of read pairs kept
        #self.remap_pair = 0
        

    def write(self, file_handle):
        sys.stderr.write("TOTAL reads: \n"
                         "  total: %d\n"
                         "DISCARD reads:\n"
                         "  unmapped: %d\n"
                         "  no SNP: %d\n"
                         "  indel: %d\n"
                         "  secondary alignment: %d\n"
                         "  supplementary alignment: %d\n"
                         "KEEP reads:\n"
                         "  single-end: %d\n\n" %
                         (self.total,
                          self.discard_unmapped,
                          self.nosnp,
                          self.discard_indel,
                          self.discard_secondary,
                          self.discard_supplementary,
                          self.keep_single))

        file_handle.write("read total: %d\n" % self.total)
        file_handle.write("read SNP ref matches: %d\n" % self.ref_count)
        file_handle.write("read SNP alt matches: %d\n" % self.alt_count)
        file_handle.write("read SNP mismatches: %d\n" % self.other_count)
        
        file_handle.write("percentage kept: %.1f%%\n" % 
            (float(self.keep_single) / float(self.total) * 100))
        
        if self.total > 0:
            mismatch_pct = 100.0 * float(self.other_count) / self.total
            if mismatch_pct > 10.0:
                sys.stderr.write("WARNING: many read SNP overlaps do not match "
                                 "either allele (%.1f%%). SNP coordinates "
                                 "in input file may be incorrect.\n" %
                                 mismatch_pct)

    
    
#######################################
### Main function: processing of single-end reads

def process_single_end(bam, snp_dir, genomes, output_dir, read_stats_file):
    """
    Function to process single-end reads, separating them based on overlapping
    SNPs. Output is written to the output_dir in two genome files
    """
    
    #########################
    ## Set parameters
    
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])
    snp_tab = snptable.SNPTable()

    cache_size = 0
    read_count = 0
    read_stats = ReadStats()
    
    
    #########################
    ## Prepare input/output
    
    # 1) Input bam
    input_bam = pysam.Samfile(bam, "r")
    
    # 2 ) Output bams
    genomes_split = genomes.split(",")
    
    basename = os.path.basename(bam).split('.')[0]
    
    output_base = os.path.join(output_dir, basename)
    
    output_file_1 = output_base + "_" + genomes_split[0] + ".bam"
    output_file_2 = output_base + "_" + genomes_split[1] + ".bam"
    
    output_bam_1 = pysam.Samfile(output_file_1, "w",
                                 template=input_bam)
    output_bam_2 = pysam.Samfile(output_file_2, "w",
                                 template=input_bam)
                                 
    
    #########################
    ## Process reads
    
    sys.stderr.write("\nProcessing reads\n")
    
    for read in input_bam:
        read_count += 1
        read_stats.total += 1
        
        if read.tid == -1:
            # unmapped read
            read_stats.discard_unmapped += 1
            continue
        
        if (cur_tid is None) or (read.tid != cur_tid):
            # this is a new chromosome
            cur_chrom = input_bam.getrname(read.tid)
            
            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
                seen_chrom.add(cur_chrom)
            
            cur_tid = read.tid            
            
            read_count = 0
            
            # Read snps
            snp_filename = "%s/%s.snps.txt.gz" % (snp_dir, cur_chrom)
            
            if not os.path.isfile(snp_filename):
                continue
            
            #sys.stderr.write("reading SNPs from file '%s'\n" % snp_filename)
            snp_tab.read_file(snp_filename)
            
            sys.stderr.write("Processing chromosome " + cur_chrom + "\n")
        
        # Determine overlapping SNPs
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos = snp_tab.get_overlapping_snps(read)
        
        if len(indel_idx) > 0:
            # for now discard this read, we want to improve this to handle
            # the indel reads appropriately
            read_stats.discard_indel += 1
            # TODO: add option to handle indels instead of throwing out reads
            continue
            
        if len(snp_idx) == 0:
            # no SNPs, no information - skip read
            read_stats.nosnp += 1
            continue
            
        ref_alleles = snp_tab.snp_allele1[snp_idx]
        alt_alleles = snp_tab.snp_allele2[snp_idx]
        
        count_genome1 = 0
        count_genome2 = 0
        count_other = 0
        
        for i, snp in enumerate(snp_idx):
            
            ref = ref_alleles[i].decode("utf-8")
            alt = alt_alleles[i].decode("utf-8")
            
            if ref == read.query_sequence[snp_read_pos[i]-1]:
                # read matches reference allele
                count_genome1 += 1
            elif alt == read.query_sequence[snp_read_pos[i]-1]:
                # read matches non-reference allele
                count_genome2 += 1
            else:
                # read matches neither ref nor other
                count_other += 1
                
        # Write write to genome is more SNPs overlap than with the other
        # genome
        if count_genome1 > count_genome2:
            output_bam_1.write(read)
            read_stats.ref_count += 1
            read_stats.keep_single += 1
        elif count_genome2 > count_genome1:
            output_bam_2.write(read)
            read_stats.alt_count += 1
            read_stats.keep_single += 1
        else: 
            read_stats.other_count += 1
        
        read_count += 1
        
    
        
    
    #########################
    ## Close files
    
    input_bam.close()
    output_bam_1.close()
    output_bam_2.close()
    
    
    #########################
    ## Write stats
    
    with open(read_stats_file, "w") as file_handle:
        read_stats.write(file_handle)
    



#######################################
### Call script

def main(bam, 
         snp_dir = None,
         genomes = None,
         output_dir = None,
         paired_end = False,
         read_stats_file = None):

    if paired_end:
        raise ValueError("paired end not supported (yet)")
    else:
        process_single_end(bam, 
                           snp_dir = snp_dir,
                           genomes = genomes,
                           output_dir = output_dir,
                           read_stats_file = read_stats_file)



if __name__ == '__main__':
    sys.stderr.write("script 'split_bam_files.py'\n\n")
    sys.stderr.write("command line: %s\n" % " ".join(sys.argv))
    sys.stderr.write("python version: %s\n" % sys.version)
    sys.stderr.write("pysam version: %s\n" % pysam.__version__)
    sys.stderr.write("pytables version: %s\n" % tables.__version__)
    
    options = parse_options()
    bam = options.bam
    
    main(options.bam,
         snp_dir = options.snp_dir,
         genomes = options.genomes,
         output_dir = options.output_dir,
         read_stats_file = options.read_stats_file)
