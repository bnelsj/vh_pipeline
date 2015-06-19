import os

### Snakefile for VariationHunter pipeline
### Assumes you have run Picard's mark duplicates and insert size calculations

### Variables that need to be set
SAMPLE_DIR = "/net/eichler/vol23/projects/human_diversity/nobackups/C_team_bams_nodups"

MANIFEST = "manifest.txt"

NODUPS_DIR = SAMPLE_DIR

READ_LEN = '101'

### Build list of samples and determine how they will be split into batches
### By default, this uses NGROUPS and not VH_GROUP_SIZE
### assigns samples to groups based on family as listed in MANIFEST

SAMPLES = []

with open(MANIFEST, "r") as reader:
    for line in reader:
        sample = line.rstrip().split()[0]
        if sample not in SAMPLES: #and (sample.startswith("OCN") or sample.startswith("AMR") or sample.startswith("SA")):
            SAMPLES.append(sample)

PICARD_ISIZE_PATH = SAMPLE_DIR
PICARD_ISIZE_METRICS = [os.path.basename(file) for file in os.listdir(PICARD_ISIZE_PATH) if file.endswith('insert_size_metrics.txt') if any(map(lambda x: file.startswith(x), SAMPLES))]

VH_GROUP_SIZE = 20
NGROUPS = 8
FAMILY_BATCHES = False

### Manifest file column names (only used if FAMILY_BATCHES = True)
FAMILY_COL_NAME = 'family'
SAMPLE_COL_NAME = 'sample'
POSITION_COL_NAME = 'position'
###

if not FAMILY_BATCHES:
    SIZE = len(SAMPLES)
    NGROUPS = 0
    while NGROUPS * VH_GROUP_SIZE < SIZE:
        NGROUPS += 1

GROUPS = [str(x).zfill(len(str(NGROUPS))) for x in range(NGROUPS)]

### Snakemake variables that probably don't need to be changed

MANIFEST_DIR = 'manifest'
ISIZEFILE = 'isizes.txt'
EDISTFILE = 'edists.txt'
DISCORDANT_READ_DIR = 'discordant_reads'
ALL_DISCO_DIR = 'all_discordant_reads'
VH_OUTDIR = 'vh_analysis'
CALL_DIR = 'calls'
MARKED_DUPS_SUFFIX = 'bam'

CONTIGS = [str(x) for x in range(1,23)] + ['X','Y'] + ['chr' + str(x) for x in range(1,23)] + ['chrX','chrY']
INCLUDE_CHRS = ':'.join(CONTIGS)

### Create directories, load modules

dirs_to_check = ['log', MANIFEST_DIR, NODUPS_DIR, VH_OUTDIR, ALL_DISCO_DIR] + [ALL_DISCO_DIR + x for x in ["/unsorted", "/sorted"]]

for dir in dirs_to_check:
    if not os.path.exists(dir):
        os.makedirs(dir)

### Begin rule definitions

rule all:
    input: expand('%s/{sample}.vh' % ALL_DISCO_DIR, sample = SAMPLES)
    params: sge_opts=""

rule get_vh_files:
    input: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s" % MANIFEST
    output: '%s/{sample}.vh' % ALL_DISCO_DIR
    benchmark: "benchmarks/{sample}_vh.json"
    params: sge_opts="-l mfree=8G -N bam2vh"
    shell:
        """python bam2vh_unpaired.py {input[0]} {input[1]} {wildcards.sample} --discordant_reads {output[0]} --discordant_read_format vh > {output}"""

rule convert_bam_to_fastq:
    input: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s/{sample}.lq.bam" % ALL_DISCO_DIR
    output: '%s/{sample}.fq' % ALL_DISCO_DIR, "%s/{sample}.lq.fq" % ALL_DISCO_DIR
    params: sge_opts="-l mfree=8G -N bam2fq"
    shell:
        """samtools bam2fq {input[0]} > {output[0]}
           samtools bam2fq {input[1]} > {output[1]}"""

rule get_all_discordant_reads:
    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX), '%s' % MANIFEST
    output: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s/{sample}.lq.bam" % ALL_DISCO_DIR
    benchmark: "benchmarks/{sample}.json"
    params: sge_opts="-l mfree=8G -N get_disco_rds -cwd -l disk_free=200G"
    shell:
        "samtools bamshuf -O {input[0]} $TMPDIR/{wildcards.sample} | python bam2vh_unpaired.py /dev/stdin {input[1]} {wildcards.sample} --discordant_reads {output[0]} --discordant_read_format bam --low_qual_reads {output[1]}"

rule get_isize_from_picard:
    input: expand('%s/{picard_isize}' % PICARD_ISIZE_PATH, picard_isize = PICARD_ISIZE_METRICS)
    output: '%s' % MANIFEST
    params: sge_opts="-N isize", edist="4"
    run:
        n_deviations = 3
        outfile = open(output[0], 'w')
        for infile in input:
            sample_name = os.path.basename(infile.replace('.insert_size_metrics.txt', ''))
            sample_path = NODUPS_DIR + '/' + sample_name + '.' + MARKED_DUPS_SUFFIX
            with open(infile, 'r') as reader:
                isize_line = False
                for line in reader:
                    if line.startswith('MEDIAN_INSERT_SIZE'):
                        isize_line = True
                        continue
                    if isize_line:
                        data = line.split()
                        median, stdev = float(data[0]), float(data[5])
                        min_isize = str(median - n_deviations * stdev)
                        max_isize = str(median + n_deviations * stdev)
                        break
            outfile.write('\t'.join([sample_name, sample_path, params.edist, min_isize, max_isize, str(READ_LEN)]) + '\n')
        outfile.close()

#rule make_manifest:
#    input: expand('%s/{sample}.txt' % MANIFEST_DIR, sample = SAMPLES)
#    output: '%s' % MANIFEST
#    params: sge_opts="-l mfree=4G -N make_manifest -cwd"
#    shell:
#        'cat {MANIFEST_DIR}/*.txt > {MANIFEST}'

#rule get_isize_from_stream:
#    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX), '%s/{sample}.%s.bai' % (NODUPS_DIR, MARKED_DUPS_SUFFIX)
#    output: '{MANIFEST_DIR}/{sample}.txt', '{MANIFEST_DIR}/{sample}.nm', '{MANIFEST_DIR}/{sample}.isize'
#    params: sge_opts="-l mfree=12G -N calc_insert_size",  n_samples='1000', n_deviations='4'
#    shell:
#         'python ~bnelsj/src/stream_read_pair/stream_sort_pairs.py --input_bam {input[0]} --n_samples {params.n_samples} --subsample_reads --binary --include_chrs {INCLUDE_CHRS} | python ~bnelsj/stream_read_pair/make_manifest_from_stream.py --input_bam {input[0]} --outdir {MANIFEST_DIR} --read_len {READ_LEN} --deviations {params.n_deviations}'

