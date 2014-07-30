import os

SAMPLE_LIST = 'samples_good.txt'
SAMPLE_DIR = '/net/eichler/vol23/projects/human_diversity/nobackups/C_team_mappings/bwa_mem_mappings'
MANIFEST_DIR = 'manifest'
OUTFILE = 'manifest.txt'
ISIZEFILE = 'isizes.txt'
EDISTFILE = 'edists.txt'
READ_LEN = 100
VH_INDIR = 'discordant_reads' 
VH_OUTDIR = 'vh_analysis'
VH_GROUP_SIZE = 30
NODUPS_DIR = 'nodups'

#CONTIGS = [str(x) for x in range(1,23)] + ['X','Y']
CONTIGS = ['chr' + str(x) for x in range(1,23)] + ['chrX','chrY']
INCLUDE_CHRS = ':'.join(CONTIGS)

SAMPLES = []

dirs_to_check = ['log', MANIFEST_DIR, NODUPS_DIR]

def makedir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

for dir in dirs_to_check:
    makedir(dir)

with open(SAMPLE_LIST, 'r') as fn:
    for line in fn.readlines():
        SAMPLES.append(os.path.basename(line).split('.')[0])

SIZE = len(SAMPLES)
NGROUPS = 0
while NGROUPS * VH_GROUP_SIZE < SIZE:
    NGROUPS += 1

SUFFIX_LIST = [str(x).zfill(len(str(NGROUPS))) for x in range(NGROUPS)]

rule all:
    input: expand('%s/{num}.SV' % VH_OUTDIR, num = SUFFIX_LIST)
    params: sge_opts=""

rule run_selection:
    input: expand('%s/{num}.{ext}' % VH_OUTDIR, num = SUFFIX_LIST, ext = ['txt', 'ReadName', 'cluster'])
    output: '%s/{num}.SV' % VH_OUTDIR
    params: sge_opts='-l mfree=60G -N vhselection', gn = '%s/{num}' % VH_OUTDIR
    shell:
        'module load VariationHunter/0.4; /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/Selection_VH_New2 -l {params.gn}.txt -r {params.gn}.ReadName -c {params.gn}.cluster -t 1000000000 -o {params.gn}.SV'

rule run_vh:
    input: '%s/{num}.txt' % VH_OUTDIR
    output: '%s/{num}.ReadName' % VH_OUTDIR, '%s/{num}.cluster' % VH_OUTDIR
    params: sge_opts="-l mfree=60G -N run_vh", gn = '%s/{num}' % VH_OUTDIR
    shell:
        'module load VariationHunter/0.4; VH -i /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/initInfo -c /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/AllChro -g /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/hg19_Gap.Table.USCS.Clean -r /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/Hg19.Satellite -l {params.gn}.txt -t {params.gn}.ReadName -o {params.gn}.cluster'

rule get_discordant_reads:
    input: '%s/{sample}/{sample}.sorted.nodups.bam' % NODUPS_DIR, 'manifest.txt'
    output: '%s/{sample}.vh' % VH_INDIR
    params: sge_opts="-l mfree=8G -N get_disco_rds -cwd", sn='{sample}'
    shell:
        'python /net/eichler/vol5/home/bnelsj/src/stream_read_pair/stream_sort_pairs.py --input_bam {input[0]} --binary --include_chrs {INCLUDE_CHRS} | python ~bnelsj/stream_read_pair/bwa_vh_pipeline/bam2vh_unpaired.py /dev/stdin {input[1]} {params.sn} > {VH_INDIR}/{params.sn}.vh 2> {VH_INDIR}/{params.sn}.vh.log'

rule prep_vh:
    input: 'manifest.txt', expand('%s/{sample}.vh' % VH_INDIR, sample = SAMPLES)
    output: '{vhdir}/{num}.txt'.format(num=num, vhdir=VH_OUTDIR) for num in SUFFIX_LIST 
    params: sge_opts='-N make_batches'
    shell:
        'python ~bnelsj/pipelines/VariationHunter/prep_divet_manifest.py --group_size {VH_GROUP_SIZE} --n_groups {NGROUPS} --manifest {input[0]} --outdir {VH_OUTDIR} --vhdir {VH_INDIR}'

rule make_manifest:
    input: expand('%s/{sample}.txt' % MANIFEST_DIR, sample = SAMPLES)
    output: 'manifest.txt'
    params: sge_opts="-l mfree=4G -N make_manifest -cwd"
    shell:
        'cat {MANIFEST_DIR}/*.txt > {OUTFILE}'

rule get_isize_from_stream:
    input: '%s/{sample}/{sample}.sorted.nodups.bam' % NODUPS_DIR, '%s/{sample}/{sample}.sorted.nodups.bam.bai' % NODUPS_DIR
    output: '{MANIFEST_DIR}/{sample}.txt', '{MANIFEST_DIR}/{sample}.nm', '{MANIFEST_DIR}/{sample}.isize'
    params: sge_opts="-l mfree=12G -N calc_insert_size",  n_samples='1000', n_deviations='4'
    shell:
         'python ~bnelsj/src/stream_read_pair/stream_sort_pairs.py --input_bam {input[0]} --n_samples {params.n_samples} --subsample_reads --binary --include_chrs {INCLUDE_CHRS} | python ~bnelsj/stream_read_pair/make_manifest_from_stream.py --input_bam {input[0]} --outdir {MANIFEST_DIR} --read_len {READ_LEN} --deviations {params.n_deviations}'

rule mark_dups:
    input: '%s/{sample}/{sample}.sorted.bam' % SAMPLE_DIR
    output: '%s/{sample}/{sample}.sorted.nodups.bam' % NODUPS_DIR, '%s/{sample}/{sample}.sorted.nodups.bam.bai' % NODUPS_DIR
    params: sge_opts = '-l mfree=8G -l disk_free=180G -N mrkdps', sn = '{sample}'
    shell:
        "java -Xmx8G -jar $PICARD_DIR/MarkDuplicates.jar INPUT={input[0]} OUTPUT={output[0]} METRICS_FILE={NODUPS_DIR}/{params.sn}/{params.sn}.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true COMPRESSION_LEVEL=5 VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=5000000 QUIET=true VERBOSITY=ERROR CREATE_INDEX=true; mv {NODUPS_DIR}/{params.sn}/{params.sn}.sorted.nodups.bai {output[1]}"
