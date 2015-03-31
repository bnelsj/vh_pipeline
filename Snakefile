import os

### Snakefile for VariationHunter pipeline
### Assumes you have run Picard's mark duplicates and insert size calculations

### Variables that need to be set
SAMPLE_DIR = "/net/eichler/vol18/genomes/moles/CHM1/UW/"

MANIFEST = "manifest.txt"

NODUPS_DIR = SAMPLE_DIR

READ_LEN = '100'

### Build list of samples and determine how they will be split into batches
### By default, this uses NGROUPS and not VH_GROUP_SIZE
### assigns samples to groups based on family as listed in MANIFEST

SAMPLES = []
with open(MANIFEST, "r") as reader:
    for line in reader:
        sample = line.rstrip().split()[0]
        if sample not in SAMPLES and sample.startswith("WEA"):
            SAMPLES.append(sample)

SAMPLES = ["CHM1_UW_Illumina.original"]

EXCLUDED_SAMPLES = ["WEA_Georgian_mg27_M", "WEA_Bishkek_28439_F"]

SAMPLES = list(set(SAMPLES) - set(EXCLUDED_SAMPLES))

PICARD_ISIZE_PATH = '.'
PICARD_ISIZE_METRICS = [os.path.basename(file) for file in os.listdir(PICARD_ISIZE_PATH) if file.endswith('insert_size_metrics.txt') if any(lambda x: file.startswith(x), SAMPLES)]

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
MARKED_DUPS_SUFFIX = 'sorted.nodups.bam'

CONTIGS = [str(x) for x in range(1,23)] + ['X','Y'] + ['chr' + str(x) for x in range(1,23)] + ['chrX','chrY']
INCLUDE_CHRS = ':'.join(CONTIGS)

### Create directories, load modules

dirs_to_check = ['log', MANIFEST_DIR, NODUPS_DIR, VH_OUTDIR, ALL_DISCO_DIR] + [DISCORDANT_READ_DIR + x for x in ["/unsorted", "/sorted"]]

for dir in dirs_to_check:
    if not os.path.exists(dir):
        os.makedirs(dir)

### Begin rule definitions

rule all:
    input: '%s/ALL.SV' % VH_OUTDIR, expand('%s/ALL.ed{ed}.ls{ls}.{type}' % CALL_DIR, ed = ["1", "2"], ls = ["3","4"], type = ["dels_per_sample", "del", "p1.denovo"])
    params: sge_opts=""

rule get_dels_per_sample:
    input: '%s/ALL.ed{ed}.ls{ls}.del' % CALL_DIR
    output: '%s/ALL.ed{ed}.ls{ls}.dels_per_sample' % CALL_DIR
    params: sge_opts='-l mfree=8G -N dels_per_sample'
    shell:
        'python ~bnelsj/pipelines/VariationHunter/get_deletions_per_sample.py {input} {output}'

rule get_p1_denovo:
    input: '%s/ALL.ed{ed}.ls{ls}.del' % CALL_DIR
    output: '%s/ALL.ed{ed}.ls{ls}.p1.denovo' % CALL_DIR
    params: sge_opts='-l mfree=8G -N denovo'    
    shell:
        'python ~bnelsj/pipelines/VariationHunter/get_denovo.py {input} {output} --manifest {MANIFEST} --family_member p1 --family_col_name {FAMILY_COL_NAME} --sample_col_name {SAMPLE_COL_NAME} --position_col_name {POSITION_COL_NAME}'

rule get_s1_denovo:
    input: '%s/ALL.ed{ed}.ls{ls}.del' % CALL_DIR
    output: '%s/ALL.ed{ed}.ls{ls}.s1.denovo' % CALL_DIR
    params: sge_opts='-l mfree=8G -N denovo'
    shell:
        'python ~bnelsj/pipelines/VariationHunter/get_denovo.py {input} {output} --manifest {MANIFEST} --family_member s1'

rule filter_deletions:
    input: '%s/ALL.SV' % VH_OUTDIR
    output: '%s/ALL.ed{ed}.ls{ls}.del' % CALL_DIR
    params: sge_opts='-l mfree=8G -N get_dels'
    shell:
        'python ~bnelsj/pipelines/VariationHunter/get_deletions.py {input} {output} --max_edist {wildcards.ed} --min_max_lib_support {wildcards.ls}'

rule combine_selections:
    input: expand('%s/{num}.SV' % VH_OUTDIR, num = GROUPS)
    output: '%s/ALL.SV' % VH_OUTDIR
    params: sge_opts='-l mfree=8G -N sel_combo'
    shell:
        'cat {input} > {output}'

rule run_selection:
    input: expand('%s/{{num}}.{ext}' % VH_OUTDIR, ext = ['txt', 'ReadName', 'cluster'])
    output: '%s/{num, \d+}.SV' % VH_OUTDIR
    params: sge_opts='-l mfree=120G -N vhselection', gn = '%s/{num}' % VH_OUTDIR
    shell:
        'module load VariationHunter/0.4; /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/Selection_VH_New2 -l {params.gn}.txt -r {params.gn}.ReadName -c {params.gn}.cluster -t 1000000000 -o {params.gn}.SV'

rule run_vh:
    input: '%s/{num}.txt' % VH_OUTDIR
    output: '%s/{num}.ReadName' % VH_OUTDIR, '%s/{num}.cluster' % VH_OUTDIR
    params: sge_opts="-l mfree=50G -N run_vh", gn = '%s/{num}' % VH_OUTDIR
    shell:
        'module load VariationHunter/0.4; VH -i /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/initInfo -c /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/AllChro -g /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/hg19_Gap.Table.USCS.Clean -r /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/Hg19.Satellite -l {params.gn}.txt -t {params.gn}.ReadName -o {params.gn}.cluster'

rule prep_vh:
    input: '%s' % MANIFEST, expand('%s/{sample}.vh' % ALL_DISCO_DIR, sample = SAMPLES)
    output: '%s/{num}.txt' % VH_OUTDIR
    params: sge_opts='-N make_batches'
    run:
        if FAMILY_BATCHES:
            shell('python ~bnelsj/pipelines/VariationHunter/prep_divet_manifest.py --group_size {VH_GROUP_SIZE} --n_groups {NGROUPS} --manifest {input[0]} --outdir {VH_OUTDIR} --vhdir {ALL_DISCO_DIR} --family {MANIFEST} --family_col_name {FAMILY_COL_NAME} --sample_col_name {SAMPLE_COL_NAME}')
        else:
            shell('python ~bnelsj/pipelines/VariationHunter/prep_divet_manifest.py --group_size {VH_GROUP_SIZE} --n_groups {NGROUPS} --manifest {input[0]} --outdir {VH_OUTDIR} --vhdir {ALL_DISCO_DIR}')

rule prep_mei:
    input: expand('%s/{sample}.{type}' % ALL_DISCO_DIR, sample = SAMPLES, type = ["fq", "vh"])
    params: sge_opts=''

rule get_vh_files:
    input: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s" % MANIFEST
    output: '%s/{sample}.vh' % ALL_DISCO_DIR
    params: sge_opts="-l mfree=8G -N bam2vh"
    shell:
        """python bam2vh_unpaired.py {input[0]} {input[1]} {wildcards.sample} --discordant_reads {output[0]} --discordant_read_format vh
           sort -k 1,1 --buffer-size=500M {ALL_DISCO_DIR}/unsorted/{wildcards.sample}.vh > {output}"""

rule convert_bam_to_fastq:
    input: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s/{sample}.lq.bam" % ALL_DISCO_DIR
    output: '%s/{sample}.fq' % ALL_DISCO_DIR, "%s/{sample}.lq.fq" % ALL_DISCO_DIR
    params: sge_opts="-l mfree=8G -N bam2fq"
    shell:
        """samtools bam2fq {input[0]} > {output[0]}"""
        """samtools bam2fq {input[1]} > {output[1]}"""

rule get_all_discordant_reads:
    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX), '%s' % MANIFEST
    output: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s/{sample}.lq.bam" % ALL_DISCO_DIR
    params: sge_opts="-l mfree=8G -N get_disco_rds -cwd"
    shell:
        "python /net/eichler/vol5/home/bnelsj/src/stream_read_pair/stream_sort_pairs.py --input_bam {input[0]} --binary --include_chrs {INCLUDE_CHRS} | python bam2vh_unpaired.py /dev/stdin {input[1]} {wildcards.sample} --low_qual_reads {output[1]} --discordant_reads {output[0]} --discordant_read_format bam"

rule get_isize_from_picard:
    input: expand('%s/{picard_isize}' % PICARD_ISIZE_PATH, picard_isize = PICARD_ISIZE_METRICS)
    output: '%s' % MANIFEST
    params: sge_opts="-N isize", n_deviations="4", edist="4"
    run:
        outfile = open(output[0], 'w')
        for infile in input:
            sample_name = os.path.basename(infile.replace('_insert_size_metrics.txt', ''))
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
                        min_isize = str(median - 4 * stdev)
                        max_isize = str(median + 4 * stdev)
                        break
            outfile.write('\t'.join([sample_name, sample_path, params.edist, min_isize, max_isize, READ_LEN]) + '\n')
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

