"""Snakemake wrapper for bwa bwasw."""

from os.path import basename
from snakemake.shell import shell
import re

# params
readgroup = snakemake.params.get('readgroup', '')
library = snakemake.params.get('library', '')
sample = snakemake.params.get('sample', '')
binary = snakemake.params.get('binary', 'bwa')
extra = snakemake.params.get('extra', '')
paired_end = snakemake.params.get('paired_end', False)
tmp_dir = snakemake.params.get('tmp_dir', '.')
reference = snakemake.params.get('reference', '')
sort_by = snakemake.params.get('sort_by', 'coordinate')
modules = snakemake.params.get("modules",dict())

# basename for temporary files
base_name = re.sub('(_R\d)$','',re.sub('\.\S+$','',basename(snakemake.input[0])))

# read group, library and sample defaults if necessary
if not readgroup:
    readgroup = base_name
if not library:
    library = readgroup
if not sample:
    sample = readgroup

# build full bwa sw command piece by piece
bwa_sw_cmd = "{binary} bwasw -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} - 2) {extra}  {snakemake.params.index_directory}/{snakemake.params.index}"

# bam input
if re.search('\.bam$',snakemake.input[0]):
    if not reference:
        raise ValueError("Need a reference path given by params.reference when using BAM files as input.")

    if paired_end==True:
        bwa_sw_cmd = bwa_sw_cmd + " <(samtools fastq -f 0x40 {snakemake.input[0]}) <(samtools fastq -f 0x80 {snakemake.input[0]})"
    else:
        bwa_sw_cmd = bwa_sw_cmd + " <(samtools fastq {snakemake.input[0]})"

    bwa_mem_cmd = bwa_sw_cmd + " | picardtools MergeBamAlignment UNMAPPED={snakemake.input[0]} ALIGNED=/dev/stdin O={snakemake.output[0]} R="+reference+" MAX_GAPS=-1 ALIGNER_PROPER_PAIR_FLAGS=true CLIP_OVERLAPPING_READS=false ADD_MATE_CIGAR=true TMP_DIR={tmp_dir}"     
    if sort_by == "none":
        bwa_sw_cmd = bwa_sw_cmd + " SO=unsorted"
    elif sort_by == "coordinate":
        bwa_sw_cmd = bwa_sw_cmd + " SO=coordinate"
    elif sort_by == "queryname":
        bwa_sw_cmd = bwa_sw_cmd + " SO=queryname"
    else:
        raise ValueError("Unexpected value for params.sort_by ({})".format(sort_by))
# fastq input
else:
    if paired_end==True:
        bwa_sw_cmd = bwa_sw_cmd + " {snakemake.input[0]} {snakemake.input[1]}" 
    else:
        bwa_sw_cmd = bwa_sw_cmd + " {snakemake.input[0]}"

    bwa_sw_cmd = bwa_sw_cmd + " | samtools addreplacerg -r \"@RG\\tID:"+readgroup+"\\tSM:"+sample+"\\tLB:"+library+"\\tCN:DeepSeqDresden\" -"
    if sort_by == "none":
        bwa_sw_cmd = bwa_sw_cmd + " | samtools view -Sbh -o {snakemake.output[0]} -"
    elif sort_by == "coordinate":
        bwa_sw_cmd = bwa_sw_cmd + " | samtools sort -O bam -@ 2 -T {tmp_dir}/"+base_name+" -o {snakemake.output[0]} -"
    elif sort_by == "queryname":
        bwa_sw_cmd = bwa_sw_cmd + " | samtools sort -O bam -n -@ 2 -T {tmp_dir}/"+base_name+" -o {snakemake.output[0]} -"
    else:
        raise ValueError("Unexpected value for params.sort_by ({})".format(sort_by))

# module load command
if "bwa" in modules and modules["bwa"]:
    modules["bwa"] = "/" + modules["bwa"]
else:
    modules["bwa"] = ""
if "samtools" in modules and modules["samtools"]:
    modules["samtools"] = "/" + modules["samtools"]
else:
    modules["samtools"] = ""
if "picardtools" in modules and modules["picardtools"]:
    modules["picardtools"] = "/" + modules["picardtools"]
else:
    modules["picardtools"] = ""
bwa_sw_cmd = "(module load apps/bwa"+modules["bwa"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true; trap \"kill 0\" SIGINT;" + bwa_sw_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+bwa_sw_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(bwa_sw_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))
