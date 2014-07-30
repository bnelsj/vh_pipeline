#!/usr/bin/python

from optparse import OptionParser

def split_manifest(manifest, group_size, n_groups, outdir, vhdir):

    samples = []
    with open(manifest, 'r') as infile:
        for line in infile:
            sample_dat = line.rstrip().split('\t')
            sample_dat[3] = str(int(float(sample_dat[3])))
            sample_dat[4] = str(int(float(sample_dat[4])))
            samples.append(sample_dat)

    sample_num = 0
    for group in range(n_groups):
        with open('%s/%d.txt' % (outdir, group), 'w') as outfile:

            if len(samples) > group_size:
                outfile.write(str(group_size) + '\n')
            else:
                outfile.write(str(len(samples)) + '\n')

            counter = 0
            while counter < group_size and len(samples) > 0:
                dat = samples.pop()
                outfile.write('\t'.join([dat[0], dat[0], '%s/%s.vh' % (vhdir, dat[0]), dat[3], dat[4], dat[5]]) + '\n')
                counter += 1

if __name__ == '__main__':
    opts = OptionParser()
    opts.add_option('','--group_size')
    opts.add_option('','--n_groups')
    opts.add_option('','--manifest')
    opts.add_option('','--outdir')
    opts.add_option('','--vhdir')

    (o, args) = opts.parse_args()
    o.group_size = int(o.group_size)
    o.n_groups = int(o.n_groups)
 
    split_manifest(o.manifest, o.group_size, o.n_groups, o.outdir, o.vhdir)
