#!/usr/bin/python

import argparse

def split_manifest_family(manifest, n_groups, outdir, vhdir, contig, family):
    samples = {}
    with open(manifest, 'r') as infile:
        for line in infile:
            sample_dat = line.rstrip().split('\t')
            sample_dat[3] = str(max(0, int(float(sample_dat[3]))))
            sample_dat[4] = str(max(0, int(float(sample_dat[4]))))
            samples[sample_dat[0]] = sample_dat

    families = {}

    with open(family, 'r') as family_info:
        for line in family_info:
            sample_dat = line.rstrip().split('\t')
            family = sample_dat[1].split('.')[0]
            if family not in families:
                families[family] = []
            families[family].append(sample_dat[0])

    batches = [[] for x in range(n_groups)]

    for i, fam in enumerate(sorted(families)):
        batch = i % n_groups
        for sample in families[fam]:
            batches[batch].append(samples[sample])

    extensions = ['txt', 'vh']

    if contig != '':
        suffix = ['%s.%s' % (contig, ext) for ext in extensions]
    else:
        suffix = extensions

    for j, batch in enumerate(batches):
        with open('%s/%d.%s' % (outdir, j, suffix[0]), 'w') as outfile:
            outfile.write(str(len(batch)) + '\n')
            for sample in batch:
                outfile.write('\t'.join([sample[0], sample[0], '%s/%s.%s' % (vhdir, sample[0], suffix[1]), sample[3], sample[4], sample[5]]) + '\n')

def split_manifest(manifest, group_size, n_groups, outdir, vhdir, contig):

    samples = []
    with open(manifest, 'r') as infile:
        for line in infile:
            sample_dat = line.rstrip().split('\t')
            sample_dat[3] = str(max(0, int(float(sample_dat[3]))))
            sample_dat[4] = str(max(0, int(float(sample_dat[4]))))
            samples.append(sample_dat)

    if group_size == 0:
        group_size = len(samples)

    extensions = ['txt', 'vh']

    if contig != '':
        suffix = ['%s.%s' % (contig, ext) for ext in extensions]
    else:
        suffix = extensions

    for group in range(n_groups):
        with open('%s/%d.%s' % (outdir, group, suffix[0]), 'w') as outfile:

            if len(samples) > group_size:
                outfile.write(str(group_size) + '\n')
            else:
                outfile.write(str(len(samples)) + '\n')

            counter = 0
            while counter < group_size and len(samples) > 0:
                dat = samples.pop()
                outfile.write('\t'.join([dat[0], dat[0], '%s/%s.%s' % (vhdir, dat[0], suffix[1]), dat[3], dat[4], dat[5]]) + '\n')
                counter += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--group_size', default=0)
    parser.add_argument('--n_groups', default=1)
    parser.add_argument('--manifest', required=True, help='Manifest file with insert size metrics')
    parser.add_argument('--outdir', default='./')
    parser.add_argument('--contig', default='', help='Chromosome number for analysis (Default: %(default)s)')
    parser.add_argument('--vhdir', required=True, help='Path to discordant read files')
    parser.add_argument('--family', default='', help='Path to manifest file with family info')

    args = parser.parse_args()

    args.group_size = int(args.group_size)
    args.n_groups = int(args.n_groups)
 
    if args.family != '':
        split_manifest_family(args.manifest, args.n_groups, args.outdir, args.vhdir, args.contig, args.family)
    else:
        split_manifest(args.manifest, args.group_size, args.n_groups, args.outdir, args.vhdir, args.contig)
