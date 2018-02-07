"""Snakemake wrapper for gsnap."""


from os.path import basename
from snakemake.shell import shell
import re

# params
readgroup = snakemake.params.get('readgroup', '')
library = snakemake.params.get('library', '')
sample = snakemake.params.get('sample', '')
binary = snakemake.params.get('binary', 'gsnap')
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

# build full gsnap command piece by piece
gsnap_cmd = "{binary} -D {snakemake.params.index_directory} -d {snakemake.params.index} -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} - 2) -A sam {extra} --use-sarray=1 --input-buffer-size=500000 --output-buffer-size=500000 -B 5"

# bam input
if re.search('\.bam$',snakemake.input[0]):
    if paired_end==True:
        gsnap_cmd = gsnap_cmd + " <(samtools fastq -f 0x40 {snakemake.input[0]}) <(samtools fastq -f 0x80 {snakemake.input[0]})"
    else:
        gsnap_cmd = gsnap_cmd + " <(samtools fastq {snakemake.input[0]})"

    gsnap_cmd = gsnap_cmd + " | picardtools MergeBamAlignment UNMAPPED={snakemake.input[0]} ALIGNED=/dev/stdin O={snakemake.output[0]} R="+reference+" MAX_GAPS=-1 ALIGNER_PROPER_PAIR_FLAGS=true CLIP_OVERLAPPING_READS=false ADD_MATE_CIGAR=true TMP_DIR={tmp_dir}"     
    if sort_by == "none":
        gsnap_cmd = gsnap_cmd + " SO=unsorted"
    elif sort_by == "coordinate":
        gsnap_cmd = gsnap_cmd + " SO=coordinate"
    elif sort_by == "queryname":
        gsnap_cmd = gsnap_cmd + " SO=queryname"
    else:
        raise ValueError("Unexpected value for params.sort_by ({})".format(sort_by))
else:
    gsnap_cmd = gsnap_cmd + " --read-group-id=\""+readgroup+"\" --read-group-name=\""+sample+"\" --read-group-library=\""+library+"\" "+("--gunzip" if re.search('\.gz$',snakemake.input[0]) else "") 

    if paired_end==True:
        gsnap_cmd = gsnap_cmd + " {snakemake.input[0]} {snakemake.input[1]}" 
    else:
        gsnap_cmd = gsnap_cmd + " {snakemake.input[0]}"

    if sort_by == "none":
        gsnap_cmd = gsnap_cmd + " | samtools view -Sbh -o {snakemake.output[0]} -"
    elif sort_by == "coordinate":
        gsnap_cmd = gsnap_cmd + " | samtools sort -O bam -@ 2 -T {tmp_dir}/"+base_name+" -o {snakemake.output[0]} -"
    elif sort_by == "queryname":
        gsnap_cmd = gsnap_cmd + " | samtools sort -O bam -n -@ 2 -T {tmp_dir}/"+base_name+" -o {snakemake.output[0]} -"
    else:
        raise ValueError("Unexpected value for params.sort_by ({})".format(sort_by))

# module load command
if "gmap" in modules and modules["gmap"]:
    modules["gmap"] = "/" + modules["gmap"]
else:
    modules["gmap"] = ""
if "samtools" in modules and modules["samtools"]:
    modules["samtools"] = "/" + modules["samtools"]
else:
    modules["samtools"] = ""
if "picardtools" in modules and modules["picardtools"]:
    modules["picardtools"] = "/" + modules["picardtools"]
else:
    modules["picardtools"] = ""
gsnap_cmd = "(module load apps/gmap"+modules["gmap"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true;" + gsnap_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+gsnap_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(gsnap_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))
