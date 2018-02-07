"""Snakemake wrapper for bwa aln/samse/sampe."""


from os.path import basename
from snakemake.shell import shell
import re

# params
readgroup = snakemake.params.get('readgroup', '')
library = snakemake.params.get('library', '')
sample = snakemake.params.get('sample', '')
binary = snakemake.params.get('binary', 'bwa')
extra_aln = snakemake.params.get('extra_aln', '')
extra_samspe = snakemake.params.get('extra_samspe', '')
paired_end = snakemake.params.get('paired_end', False)
tmp_dir = snakemake.params.get('tmp_dir', '.')
reference = snakemake.params.get('reference', '')
sort_by = snakemake.params.get('sort_by', 'coordinate')
default_threads = snakemake.threads
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

# build full bwa aln command piece by piece

# bam input
if re.search('\.bam$',snakemake.input[0]):
    if not reference:
        raise ValueError("Need a reference path given by params.reference when using BAM files as input.")
    if paired_end:
        bwa_aln_cmd = "{binary} sampe {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}"
        bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -b1 -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} / 2 - 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]})"
        bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -b2 -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} / 2 - 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]})"
        bwa_aln_cmd = bwa_aln_cmd + " {snakemake.input[0]} {snakemake.input[0]}"
    else:
        bwa_aln_cmd = "{binary} samse {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}" 
        bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -b -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} - 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]}) {snakemake.input[0]}"

    bwa_aln_cmd = bwa_aln_cmd + " | picardtools MergeBamAlignment UNMAPPED={snakemake.input[0]} ALIGNED=/dev/stdin O={snakemake.output[0]} R="+reference+" MAX_GAPS=-1 ALIGNER_PROPER_PAIR_FLAGS=true CLIP_OVERLAPPING_READS=false ADD_MATE_CIGAR=true TMP_DIR={tmp_dir}"
    if sort_by == "none":
        bwa_mem_cmd = bwa_mem_cmd + " SO=unsorted"
    elif sort_by == "coordinate":
        bwa_mem_cmd = bwa_mem_cmd + " SO=coordinate"
    elif sort_by == "queryname":
        bwa_mem_cmd = bwa_mem_cmd + " SO=queryname"
    else:
        raise ValueError("Unexpected value for params.sort_by ({})".format(sort_by))
# fastq input
else:
    if paired_end:
        bwa_aln_cmd = "{binary} sampe -r \"@RG\\tID:"+readgroup+"\\tSM:"+sample+"\\tLB:"+library+"\\tCN:DeepSeqDresden\" {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}"
        bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} / 2 - 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]})"
        bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} / 2 - 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[1]})"
        bwa_aln_cmd = bwa_aln_cmd + " {snakemake.input[0]} {snakemake.input[1]}"
    else:
        bwa_aln_cmd = "{binary} samse -r \"@RG\\tID:"+readgroup+"\\tSM:"+sample+"\\tLB:"+library+"\\tCN:DeepSeqDresden\" {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}" 
        bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} - 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]}) {snakemake.input[0]}"
    
    if sort_by == "none":
        bwa_aln_cmd = bwa_aln_cmd + " | samtools view -Sbh -o {snakemake.output[0]} -"
    elif sort_by == "coordinate":
            bwa_aln_cmd = bwa_aln_cmd + " | samtools sort -O bam -@ 2 -T {tmp_dir}/"+base_name+" -o {snakemake.output[0]} -"
    elif sort_by == "queryname":
            bwa_aln_cmd = bwa_aln_cmd + " | samtools sort -O bam -n -@ 2 -T {tmp_dir}/"+base_name+" -o {snakemake.output[0]} -"
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
# note the trap is neccessary to kill the subshells <(..) when the parent shell dies
bwa_aln_cmd = "(module load apps/bwa"+modules["bwa"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true; trap \"kill 0\" SIGINT;" + bwa_aln_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+bwa_aln_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(bwa_aln_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))


