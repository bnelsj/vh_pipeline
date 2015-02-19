import argparse
import pandas as pd

def get_svtype(entries):
    for entry in entries:
        if 'SVtype:' in entry:
            svtype = entry.split(':')[1]
            return svtype

def get_edist(entries):
    for entry in entries:
        if 'AvgEditDits:' in entry:
            edist = float(entry.split(':')[1])
            return edist

def get_max_lib_support(entries):
    max_lib_support = 0
    for entry in entries:
        if 'LibSup:' in entry:
            lib_sup = int(entry.split(':')[1])
            max_lib_support = max(max_lib_support, lib_sup)
    return max_lib_support

def add_dels(entries, dels_per_sample):
    for i, entry in enumerate(entries):
        if entries[i].startswith("Lib:"):
            sn = entries[i].split(':')[1]
            lib_support = int(entries[i+1].split(':')[1])
            if lib_support > 0:
                if sn in dels_per_sample:
                    dels_per_sample[sn] += 1
                else:
                    dels_per_sample[sn] = 1
    return dels_per_sample

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--max_edist', default=2, help='Maximum mean edit distance to count deletion (Default: %(default)s)')
    parser.add_argument('--min_max_lib_support', default=4, help='Minimum max library support to count deletion (Default: %(default)s)')

    args = parser.parse_args()

    dels_per_sample = {}

    max_edist = float(args.max_edist)
    min_max_lib_support = float(args.min_max_lib_support)

    with open(args.infile, 'r') as infile:
        for line in infile:
            entries = line.split()
            if get_svtype(entries) == 'D' and get_edist(entries) <= max_edist and get_max_lib_support(entries) >= min_max_lib_support:
                dels_per_sample = add_dels(entries, dels_per_sample)

    with open(args.outfile, 'w') as outfile:
        for sample, del_count in sorted(dels_per_sample.iteritems()):
            outfile.write("%s\t%s\n" % (sample, del_count))
