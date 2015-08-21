import os

### Snakefile for VariationHunter pipeline
### Assumes you have run Picard's mark duplicates and insert size calculations

### Variables that need to be set
SAMPLE_DIR = "/net/eichler/vol23/projects/human_diversity/nobackups/bnelsj/vh_hgdp_OCN/vh_pipeline/all_discordant_reads"

MANIFEST = "manifest.txt"

NODUPS_DIR = SAMPLE_DIR

READ_LEN = "100"
EDIST = 4
N_DEV = 4

WHAM_PATH = "/net/eichler/vol5/home/bnelsj/src/wham"
REFERENCE_FASTA = "/net/eichler/vol2/eee_shared/assemblies/human_1kg_v37/human_1kg_v37.fasta"
REFERENCE_GC_PROFILE = "/net/eichler/vol2/eee_shared/assemblies/hg19/gc_profile/ucsc.hg19.gc_profile.bed"

ruleorder: get_isize_from_wham > get_isize_from_picard
#ruleorder: get_isize_from_picard > get_isize_from_wham

### Build list of samples and determine how they will be split into batches
### By default, this uses NGROUPS and not VH_GROUP_SIZE
### assigns samples to groups based on family as listed in MANIFEST

SAMPLES = []
SAMPLE_SUFFIX = "bam"

for file in os.listdir(SAMPLE_DIR):
    if file.endswith(SAMPLE_SUFFIX) and not file.endswith(".lq.bam"):
        SAMPLES.append(file.replace("." + SAMPLE_SUFFIX, ""))

ISIZE_PATH = "/net/eichler/vol23/projects/human_diversity/nobackups/bnelsj/hgdp_remapped_isizes/vh_pipeline/isizes"
PICARD_ISIZE_SUFFIX = "insert_size_metrics.txt"
WHAM_ISIZE_SUFFIX = "wham_isize.txt"

VH_GROUP_SIZE = 22
NGROUPS = 1
FAMILY_BATCHES = False

### Manifest file column names (only used if FAMILY_BATCHES = True)
# FAMILY_MANIFEST must be a tab-delimited file with columns for sample and family
FAMILY_MANIFEST = "/path/to/family/manifest"
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
SUFFIX_LIST = [str(x).zfill(len(str(NGROUPS))) for x in range(NGROUPS)]

### Snakemake variables that probably don't need to be changed

ISIZEFILE = 'isizes.txt'
EDISTFILE = 'edists.txt'
DISCORDANT_READ_DIR = 'discordant_reads'
ALL_DISCO_DIR = 'all_discordant_reads'
VH_OUTDIR = 'vh_analysis'
READ_DEPTH_DIR = "read_depth"
CALL_DIR = 'calls'
MARKED_DUPS_SUFFIX = 'bam'

CHR_CONTIGS = ['chr' + str(x) for x in range(1,23)] + ['chrX','chrY']
CONTIGS = [str(x) for x in range(1,23)] + ['X','Y'] + CHR_CONTIGS
INCLUDE_CHRS = ':'.join(CONTIGS)

### Create directories, load modules

dirs_to_check = ['log', NODUPS_DIR, VH_OUTDIR, ALL_DISCO_DIR, "svs", "calls", "read_depth", "depth", "isizes"]

def _get_family_string(wildcards):
	if FAMILY_BATCHES:
		return "--family %s --family_col_name %s --sample_col_name %s" % (FAMILY_MANIFEST, FAMILY_COL_NAME, SAMPLE_COL_NAME)
	else:
		return ""

for dir in dirs_to_check:
    if not os.path.exists(dir):
        os.makedirs(dir)

### Begin rule definitions

rule all:
    input: "final_calls.bed" #expand("svs/{chr}.SV.DEL.merged.MEI", chr = CHR_CONTIGS)
    params: sge_opts=""

rule genotype_samples: # max SV list size set to 5 million for svs/ALL.SV.DEL.merged.
    input: "samples.txt", "svs/ALL.SV.DEL.merged",  "calls/Alu_L1_SV_Picked.txt", "proband_list.txt", "depth_file_manifest.txt", "calls/VH_calls_gt500bp.txt"
    output: "final_calls.bed"
    params: sge_opts = "-l mfree=64G"
    shell:
        "./bin/genotype_MultipleSamples2 {input} > {output}"

rule get_proband_samples:
    input: "samples.txt"
    output: "proband_list.txt"
    params: sge_opts = ""
    shell: "touch {output}"

rule get_depth_file_manifest:
    input: expand("%s/VH_calls_gt500bp.{sample}.bam.Depth" % "depth", sample = SAMPLES)
    output: "depth_file_manifest.txt"
    params: sge_opts = ""
    run:
        with open(output[0], "w") as outfile:
            for sn in SAMPLES:
                infile = [os.path.abspath(file) for file in input if sn in file][0]
                outfile.write("%s\t%s\n" % (infile, sn))


rule get_gc_corrected_read_depth_per_call:
    input: REFERENCE_GC_PROFILE, "%s/{sample}.bam.Depth" % READ_DEPTH_DIR, "calls/VH_calls_gt500bp.clean"
    output: "%s/VH_calls_gt500bp.{sample}.bam.Depth" % "depth"
    params: sge_opts = "-l mfree=64G"
    shell:
        "./bin/calculateReadDepthFromBAM {input} > {output}"

rule get_depth:
    input: expand("%s/{sample}.bam.Depth" % READ_DEPTH_DIR, sample = SAMPLES)
    params: sge_opts = ""

rule get_read_depth:
    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX)
    output: "%s/{sample}.bam.Depth" % READ_DEPTH_DIR
    params: sge_opts = "", tmpfile = "$TMPDIR/{sample}.bam.Depth"
    shell:
        "samtools depth {input} > {params.tmpfile}; "
        "rsync --bwlimit 10000 {params.tmpfile} {output}"

rule get_clean_gt500bp_file:
    input: "calls/VH_calls_gt500bp.tab"
    output: "calls/VH_calls_gt500bp.clean"
    params: sge_opts = ""
    shell: "cut -f 1-4 {input} > {output}"

rule get_gt500bp_call_names:
    input: "calls/VH_calls_gt500bp.tab"
    output: "calls/VH_calls_gt500bp.txt"
    params: sge_opts = ""
    shell: "cut -f 4 {input} > {output}"

rule get_calls_over_500bp:
    input: "svs/ALL.SV.DEL.merged"
    output: "calls/VH_calls_gt500bp.tab"
    params: sge_opts = ""
    shell:
        "cat {input} | awk '{{if ($3-$2>500) print;}}' > {output}"

rule combine_all_calls:
    input: expand("svs/{chr}.SV.DEL.merged", chr = CHR_CONTIGS)
    output: "svs/ALL.SV.DEL.merged"
    params: sge_opts = ""
    shell:
        "cat {input} > {output}"

rule get_picked_mei:
    input: expand("svs/{chr}.SV.DEL.merged.MEI", chr = CHR_CONTIGS)
    output: "calls/Alu_L1_SV_Picked.txt"
    params: sge_opts = ""
    shell:
        "cat {input} | cut -f 4 > {output}"

rule get_mei_svs: # Currently set for max of 160 samples
    input: "svs/{chr}.SV.DEL.merged"
    output: "svs/{chr}.SV.DEL.merged.MEI"
    params: sge_opts = "-l mfree=8G", mei = "/net/eichler/vol1/home/fhormozd/1000Genome_SV/Alu_L1_Del/Stewart_Alu_L1_Del.hg19"
    shell: "bin/spansKnownME {params.mei} {input} > {output}"

rule get_merged:
    input: expand("svs/{chr}.SV.DEL.merged", chr = CHR_CONTIGS)
    params: sge_opts = "-l mfree=8G"

rule merge_samples: # Need to test mergeSamples with different number of samples
    input: "svs/ALL.SV.{chr}.DEL", "samples.txt"
    output: "svs/{chr}.SV.DEL.merged"
    params: sge_opts = "-l mfree=8G", nsamples = str(len(SAMPLES))
    run:
        with open(input[0]) as f:
            for i, l in enumerate(f):
                pass
            num_sv = str(i + 1)
        shell("bin/mergeSamples {input} {num_sv} {params.nsamples} > {output}")

rule make_sample_list_file:
    input: MANIFEST
    output: "samples.txt"
    params: sge_opts = ""
    shell: "cut -f 1 {input} > {output}"

rule filter_deletions_by_chr:
    input: "%s/ALL.SV.{chr}" % VH_OUTDIR
    output: "%s/ALL.SV.{chr}.DEL" % "svs"
    params: sge_opts = ""
    shell:
        "grep SVtype:D {input[0]} > {output}"

rule combine_sv:
    input: expand("%s/{num}.SV.{chr}" % VH_OUTDIR, num = SUFFIX_LIST, chr = CHR_CONTIGS)
    output: "%s/ALL.SV.{chr}" % VH_OUTDIR
    params: sge_opts = "", chr = "{chr}"
    run:
        files = []
        for file in input:
            if file.endswith(params.chr):
                files.append(file)
        shell("cat {files} > {output}")

rule run_selection:
    input: expand('%s/{{num}}.{ext}' % VH_OUTDIR, ext = ['txt', 'ReadName']), '%s/{num}.cluster.{chr}' % VH_OUTDIR
    output: '%s/{num}.SV.{chr}' % VH_OUTDIR
    params: sge_opts='-l mfree=80G -N vhselection', gn = '%s/{num}' % VH_OUTDIR
    shell:
        'module load VariationHunter/0.4; /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/Selection_VH_New2 -l {params.gn}.txt -r {params.gn}.ReadName -c {params.gn}.cluster.{wildcards.chr} -t 1000000000 -o {output}'

rule split_clust_by_chr:    
    input: "%s/{num}.cluster" % VH_OUTDIR
    output: "%s/{num}.cluster.{chr}" % VH_OUTDIR
    params: sge_opts="-l mfree=8G", chr = "{chr}"
    shell:
        "grep -w {params.chr} {input} > {output}"

rule run_vh:
    input: '%s/{num}.txt' % VH_OUTDIR
    output: '%s/{num}.ReadName' % VH_OUTDIR, '%s/{num}.cluster' % VH_OUTDIR
    params: sge_opts="-l mfree=80G -N run_vh", gn = '%s/{num}' % VH_OUTDIR
    shell:
        'module load VariationHunter/0.4; VH -i /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/initInfo -c /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/AllChro -g /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/hg19_Gap.Table.USCS.Clean -r /net/eichler/vol5/home/bknelson/src/Hg19_NecessaryFiles/Hg19.Satellite -l {params.gn}.txt -t {params.gn}.ReadName -o {params.gn}.cluster'

rule prep_vh:
    input: 'manifest.txt', expand('%s/{sample}.vh' % ALL_DISCO_DIR, sample = SAMPLES)
    output: '{vhdir}/{num}.txt'.format(num=num, vhdir=VH_OUTDIR) for num in SUFFIX_LIST
    params: sge_opts='-N make_batches', family_string = _get_family_string
    shell:
        'python prep_divet_manifest.py --group_size {VH_GROUP_SIZE} --n_groups {NGROUPS} --manifest {input[0]} --outdir {VH_OUTDIR} --vhdir {ALL_DISCO_DIR} {params.family_string}'

rule do_get_vh_files:
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
    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX), MANIFEST
    output: '%s/{sample}.bam' % ALL_DISCO_DIR, "%s/{sample}.lq.bam" % ALL_DISCO_DIR
    benchmark: "benchmarks/{sample}.json"
    params: sge_opts="-l mfree=8G -N get_disco_rds -cwd -l disk_free=200G"
    shell:
        "samtools bamshuf -O {input[0]} $TMPDIR/{wildcards.sample} | python bam2vh_unpaired.py /dev/stdin {input[1]} {wildcards.sample} --discordant_reads {output[0]} --discordant_read_format bam --low_qual_reads {output[1]}"

rule get_isize_from_wham:
    input: expand("%s/{sample}.%s" % (ISIZE_PATH, WHAM_ISIZE_SUFFIX), sample = SAMPLES)
    output: "%s" % MANIFEST
    params: sge_opts = "-N isize"
    run:
        outfile = open(output[0], "w")
        for infile in input:
            sample_name = os.path.basename(infile.replace("." + WHAM_ISIZE_SUFFIX, ''))
            sample_path = NODUPS_DIR + '/' + sample_name + '.' + MARKED_DUPS_SUFFIX
            with open(infile, 'r') as reader:
                median_isize, sd_isize = None, None
                for line in reader:
                    if "median insert length" in line:
                        median_isize = float(line.rstrip().split()[-1])
                    elif "sd insert length" in line:
                        sd_isize = float(line.rstrip().split()[-1])
                    if median_isize is not None and sd_isize is not None:
                        min_isize = str(int(median_isize - N_DEV * sd_isize))
                        max_isize = str(int(median_isize + N_DEV * sd_isize))
                        break
            outfile.write('\t'.join([sample_name, sample_path, str(EDIST), min_isize, max_isize, str(READ_LEN)]) + '\n')
        outfile.close()
       

rule get_isize_from_picard:
    input: expand("%s/{sample}.%s" % (ISIZE_PATH, PICARD_ISIZE_SUFFIX), sample = SAMPLES)
    output: '%s' % MANIFEST
    params: sge_opts="-N isize"
    run:
        outfile = open(output[0], 'w')
        for infile in input:
            sample_name = os.path.basename(infile.replace("." + PICARD_ISIZE_SUFFIX, ''))
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
                        min_isize = str(median - N_DEV * stdev)
                        max_isize = str(median + N_DEV * stdev)
                        break
            outfile.write('\t'.join([sample_name, sample_path, str(EDIST), min_isize, max_isize, str(READ_LEN)]) + '\n')
        outfile.close()

rule calc_isize_picard:
    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX)
    output: "%s/{sample}.%s" % (ISIZE_PATH, PICARD_ISIZE_SUFFIX)
    params: sge_opts = "-l mfree=8G -N isize_picard"
    shell:
        """module load java/7u17 picard/1.111; """
        """java -Xmx8G -jar $PICARD_DIR/CollectInsertSizeMetrics.jar I={input} O={output} H={output}.hist"""

rule calc_isize_wham:
    input: '%s/{sample}.%s' % (NODUPS_DIR, MARKED_DUPS_SUFFIX), REFERENCE_FASTA
    output: "%s/{sample}.%s" % (ISIZE_PATH, WHAM_ISIZE_SUFFIX)
    params: sge_opts = "-l mfree=8G -N isize_wham"
    shell:
        """{WHAM_PATH}/bin/WHAM-GRAPHENING -s -f {input[0]} -a {input[1]} 2> {output[0]}"""
