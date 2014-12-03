import argparse

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--max_edist', default=2, help='Maximum edit distance to count deletion (Default: %(default)s)')
    parser.add_argument('--min_max_lib_support', default=4, help='Minimum max library support to count deletion (Default: %(default)s)')

    args = parser.parse_args()

    max_edist = float(args.max_edist)
    min_max_lib_support = float(args.min_max_lib_support)

    with open(args.infile, 'r') as infile:
        with open(args.outfile, 'w') as outfile:
            for line in infile:
                entries = line.split()
                if get_svtype(entries) == 'D' and get_edist(entries) <= max_edist and get_max_lib_support(entries) >= min_max_lib_support:
                    outfile.write(line)
