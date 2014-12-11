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

def is_denovo(entries, sample_list, family_member, position_col):
    denovo = False
    for i, entry in enumerate(entries):
        if entries[i].startswith("Lib:"):
            sn = entries[i].split(':')[1]
            lib_support = int(entries[i+1].split(':')[1])
            if lib_support > 0:
                if denovo:
                    return False
                fam_pos = sample_list.ix[sn, position_col]
                if fam_pos == family_member:
                    denovo = True
                else:
                    return False
    if denovo:
        return True
    else:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--max_edist', default=2, help='Maximum edit distance to count deletion (Default: %(default)s)')
    parser.add_argument('--min_max_lib_support', default=4, help='Minimum max library support to count deletion (Default: %(default)s)')
    parser.add_argument('--manifest', required=True, help='Path to sample manifest file')
    parser.add_argument('--family_col_name', default='family', help='Name of manifest column with family info')
    parser.add_argument('--sample_col_name', default='sample', help='Name of manifest column with sample info')
    parser.add_argument('--position_col_name', default='position', help='Name of manifest column with family position (p1, s1, fa, mo) info')
    parser.add_argument('--sex_col_name', default='sex', help='Name of manifest column with sex info')
    parser.add_argument('--family_member', choices=['p1','s1'], default='p1')

    args = parser.parse_args()

    sample_dict = {}

    sample_list = pd.read_csv(args.manifest, sep='\t')
    sample_names = sample_list[args.sample_col_name]
    sample_list.index = sample_names

    max_edist = float(args.max_edist)
    min_max_lib_support = float(args.min_max_lib_support)

    with open(args.infile, 'r') as infile:
        with open(args.outfile, 'w') as outfile:
            for line in infile:
                entries = line.split()
                if get_svtype(entries) == 'D' and get_edist(entries) <= max_edist and get_max_lib_support(entries) >= min_max_lib_support:
                    if is_denovo(entries, sample_list, args.family_member, args.position_col_name):
                        outfile.write(line)
