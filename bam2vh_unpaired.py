#!/usr/bin/env python

import sys
import os
import pysam
import math
import argparse
import gc

MINIMUM_MAPPING_QUALITY = 20

def check_read_pair(read1, read2):
    if read1.qname != read2.qname:
        return "NAME_MISMATCH"
    
    if read1.tid != read2.tid:
        return "CONTIG_MISMATCH"

    if read1.is_duplicate or read2.is_duplicate:
        return "PCR_DUPLICATE"

    if read1.is_qcfail or read2.is_qcfail:
        return "QC_FAIL"

    if(len(read1.qual) == 0 or len(read2.qual) == 0):
        return "NO_QUALITY_STRING"

    min_quality = min(read1.mapping_quality, read2.mapping_quality)

    if min_quality < MINIMUM_MAPPING_QUALITY:
        return "LOW_MAPPING_QUALITY"

    if read1.is_unmapped and read2.is_unmapped:
        return "UNMAPPED_READS"
    elif read1.is_unmapped or read2.is_unmapped:
        return "ONE_UNMAPPED_READ"

    if(read1.pos == read2.pos):
        return "READS_OVERLAP"

    if(len(read1.seq) != len(read2.seq)):
        return "SIZE_MISMATCH"

    if(len(read1.seq) - math.fabs(read1.pos - read2.pos) > 2):
        return "READS_OVERLAP"
    
    return "GOOD"

class VH:
    
    def __init__(self,read1,read2,ins_min,ins_max):
        self.qname = read1.qname
        self.tid = read1.tid
        #print str(samfile.getrname(self.contig))

        self.start1 = read1.pos
        self.start2 = read2.pos
        self.qlen = len(read1.seq)
       
        self.end1 = self.qlen + self.start1
        self.end2 = self.qlen + self.start2

        edit_list = [x[1] for x in read1.tags if x[0] == 'NM']
        self.edit1 = edit_list[0] if len(edit_list) > 0 else 0
        edit_list = [x[1] for x in read2.tags if x[0] == 'NM']
        self.edit2 = edit_list[0] if len(edit_list) > 0 else 0
        
        self.ins1 = sum([x[1] for x in read1.cigar if x[0] == 1]) 
        self.ins2 = sum([x[1] for x in read2.cigar if x[0] == 1])
        
        self.del1 = sum([x[1] for x in read1.cigar if x[0] == 2])
        self.del2 = sum([x[1] for x in read2.cigar if x[0] == 2])

        self.snp1 = max(0,self.edit1 - self.ins1 - self.del1)
        self.snp2 = max(0,self.edit2 - self.ins2 - self.del2)
        self.min_quality = min(read1.mapping_quality, read2.mapping_quality)

        self.phred_avg1 = float(sum([ord(x) - 33 for x in read1.qual]))/len(read1.qual)
        self.phred_avg2 = float(sum([ord(x) - 33 for x in read2.qual]))/len(read2.qual)

        self.phred_avg = float(sum([ord(x) - 33 for x in read1.qual]) + sum([ord(y) - 33 for y in read2.qual]))/(len(read1.qual) + len(read2.qual))

        self.rev1 = read1.is_reverse
        self.rev2 = read2.is_reverse

        self.first1 = self.start1 < self.start2
        self.first2 = self.start2 < self.start1
        
        self.ins_min = ins_min
        self.ins_max = ins_max

        if self.first1:
            self.isize = math.fabs(read2.pos + read2.qend - (read1.pos + read1.qstart))
        else:
            self.isize =  math.fabs(read1.pos + read1.qend - (read2.pos + read2.qstart))
        
    def get_overlap(self):
        return self.qlen - math.fabs(self.start1 - self.start2)

    def get_phred_prob(self):
        return self.__get_phred_prob_read(1) * self.__get_phred_prob_read(2)

    def __get_phred_prob_read(self,readnum):
        if(readnum == 1):
            return self.__calc_phred_prob(self.ins1+self.del1,self.snp1,self.phred_avg1)
        else:
            return self.__calc_phred_prob(self.ins2+self.del2,self.snp2,self.phred_avg2)

    def __calc_phred_prob(self,indel,snp,avg):
        return self.__calc_phred_snp(snp,avg) * self.__calc_phred_indel(indel,avg)
    
    def __calc_phred_snp(self,snp,avg):
        return math.pow((float(1)/1000 + math.pow(10,(float(avg) * -1 )/10)),snp)
    
    def __calc_phred_indel(self,indel,avg):
        return math.pow((float(1)/7000 + math.pow(10,(float(avg) * -1 )/10)),indel)

    def entry(self):
        return "\t".join([self.qname,self.get_str_contig(),str(self.start1),str(self.end1),self.get_orientation(self.rev1),"=",str(self.start2),str(self.end2),self.get_orientation(self.rev2),self.get_event(),str(self.edit1+self.edit2),str(self.phred_avg),"%e" % self.get_phred_prob()])

    def get_orientation(self,rev):
        if rev:
            return 'R'
        else:
            return 'F'

    def get_event(self):
        if(self.check_inversion()):
            return 'V'
        elif(self.check_eversion()):
            return 'E'
        elif(self.check_deletion()):
            return 'D'
        elif(self.check_insertion()):
            return 'I'
        elif(self.check_translocation()): #always returns false, under construction
            return 'T'
        else:
            return 'C'    

    def is_high_quality(self):
        return self.min_quality >= MINIMUM_MAPPING_QUALITY

    def check_inversion(self): #Orientation of the reads are the same 
        return self.rev1 == self.rev2

    def check_eversion(self):
        return self.rev1 != self.rev2 and self.first1 == self.rev1 #Orientations are reversed and either read1 is first and reverse or not first and not reverse 

    def check_deletion(self):
        return self.isize > self.ins_max

    def check_insertion(self):
        return self.isize < self.ins_min

    def check_translocation(self): #always returns false, under construction
        return False

    def get_str_contig(self):
        contig_name = samfile.getrname(self.tid)
        return contig_name if contig_name.startswith('chr') else 'chr' + contig_name

def want_contig(contig):
    return True

def is_empty(file):
    test = pysam.AlignmentFile(file, check_sq=False)
    try:
        read = test.next()
        empty = False
    except StopIteration:
        empty = True
        os.remove(file)
    test.close()
    return empty

def get_next_paired_read(samfile):
    read_a = samfile.next()
    read_b = samfile.next()
    while not read_a.query_name == read_b.query_name:
        read_a = read_b
        read_b = samfile.next()
    return read_a, read_b

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="The sam filename")
    parser.add_argument('manifest', help="Path to sample manifest")
    parser.add_argument('sample_name')
    parser.add_argument('--debug', action='store_true', help='Debug mode with verbose error logging')
    parser.add_argument("--low_qual_reads", default=None, help="Output low quality reads to specified file")
    parser.add_argument("--discordant_reads", default = "/dev/stdout", help="Output discordant reads to specified file (default: %(default)s)")
    parser.add_argument("--discordant_read_format", choices=["bam", "vh"], default = "vh", help="Default: %(default)s")
    parser.add_argument("--check_sq", action="store_true", help="Check SQ flag in input bam header")
    args = parser.parse_args()

    ACCEPTABLE_TYPES = ["GOOD"]

    if args.low_qual_reads is not None:
        ACCEPTABLE_TYPES.extend(["LOW_MAPPING_QUALITY", "ONE_UNMAPPED_READ"])

    with open(args.manifest, 'r') as manfile:
        for line in manfile:
            entries = line.rstrip().split("\t")
            if args.sample_name.startswith(entries[0]):
                min_isize = float(entries[3])
                max_isize = float(entries[4])
                break

    try:
        assert(min_isize < max_isize)
    except NameError:
        sys.stderr.write('Sample name %s not found in manifest %s\n' % (args.sample_name, args.manifest))
        sys.exit(1)
    except AssertionError:
        sys.stderr.write('Min isize %d must be less than max isize %d.\n' % (min_isize, max_isize))
        sys.exit(1)

    samfile = pysam.Samfile(args.filename, check_sq = args.check_sq)

    if args.low_qual_reads is not None:
        lq_file = pysam.Samfile(args.low_qual_reads, 'wb', template=samfile)

    if args.discordant_read_format == "bam":
        discordant_file = pysam.Samfile(args.discordant_reads, 'wb', template=samfile)
    else:
        discordant_file = open(args.discordant_reads, "w")

    while True:
        try:
            read_a, read_b = get_next_paired_read(samfile)
        except StopIteration:
            break

        type = check_read_pair(read_a, read_b)

        if type in ACCEPTABLE_TYPES:
            if type != "GOOD":
                if args.low_qual_reads is not None:
                    lq_file.write(read_a)
                    lq_file.write(read_b)
            else:
                vh_entry = VH(read_a, read_b, min_isize, max_isize)
                if vh_entry.get_event() != "C":
                    if args.discordant_read_format == "bam":
                        discordant_file.write(read_a)
                        discordant_file.write(read_b)
                    else:
                        discordant_file.write(vh_entry.entry() + "\n")
        else:
            if args.debug:
                sys.stderr.write("%s: %s, %s\n" % (type, read_a.qname, read_b.qname))

    samfile.close()
    discordant_file.close()

    lq_empty = False
    if args.low_qual_reads is not None:
        lq_file.close()
        lq_empty = is_empty(args.low_qual_reads)

    if is_empty(discordant_file) or lq_empty:
        sys.exit(1)
