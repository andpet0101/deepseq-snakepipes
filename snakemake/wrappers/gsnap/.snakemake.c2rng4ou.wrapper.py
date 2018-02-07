
######## Snakemake header ########
import sys; sys.path.insert(0, "/share/apps/python/3.4.2/lib/python3.4/site-packages/snakemake-3.12.0-py3.4.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x07\x00\x00\x00threadsq\x03K\x03X\t\x00\x00\x00resourcesq\x04csnakemake.io\nResources\nq\x05)\x81q\x06(K\x01K\x03e}q\x07(X\x06\x00\x00\x00_coresq\x08K\x03X\x06\x00\x00\x00_nodesq\tK\x01X\x06\x00\x00\x00_namesq\n}q\x0b(h\tK\x00N\x86q\x0ch\x08K\x01N\x86q\ruubX\x04\x00\x00\x00ruleq\x0eX\t\x00\x00\x00run_gsnapq\x0fX\x05\x00\x00\x00inputq\x10csnakemake.io\nInputFiles\nq\x11)\x81q\x12X\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzq\x13a}q\x14h\n}q\x15sbX\x06\x00\x00\x00paramsq\x16csnakemake.io\nParams\nq\x17)\x81q\x18(X\x0e\x00\x00\x003PrimeCaoTest2q\x19X&\x00\x00\x00/projects/seq-work/user/pipeline/gmap/q\x1aX\x0c\x00\x00\x00${TMPDIR:-.}q\x1bX\x06\x00\x00\x00L25037q\x1cX\x05\x00\x00\x00gsnapq\x1dX\x15\x00\x00\x00L25037_3PrimeCaoTest2q\x1eX\x04\x00\x00\x00hg38q\x1f}q X\n\x00\x00\x00coordinateq!X\x00\x00\x00\x00q"e}q#(X\x06\x00\x00\x00sampleq$h\x19X\x07\x00\x00\x00modulesq%h X\x0f\x00\x00\x00index_directoryq&h\x1aX\x07\x00\x00\x00tmp_dirq\'h\x1bX\x05\x00\x00\x00indexq(h\x1fX\x06\x00\x00\x00binaryq)h\x1dX\t\x00\x00\x00readgroupq*h\x1eh\n}q+(h$K\x00N\x86q,h%K\x07N\x86q-h&K\x01N\x86q.h\'K\x02N\x86q/h)K\x04N\x86q0h*K\x05N\x86q1h(K\x06N\x86q2X\x07\x00\x00\x00libraryq3K\x03N\x86q4X\x07\x00\x00\x00sort_byq5K\x08N\x86q6X\x05\x00\x00\x00extraq7K\tN\x86q8uh3h\x1ch5h!h7h"ubX\t\x00\x00\x00wildcardsq9csnakemake.io\nWildcards\nq:)\x81q;X\x15\x00\x00\x00L25037_3PrimeCaoTest2q<a}q=(X\x08\x00\x00\x00basenameq>h<h\n}q?X\x08\x00\x00\x00basenameq@K\x00N\x86qAsubX\x03\x00\x00\x00logqBcsnakemake.io\nLog\nqC)\x81qDX#\x00\x00\x00gsnap/log/L25037_3PrimeCaoTest2.logqEa}qFh\n}qGsbX\x06\x00\x00\x00configqH}qI(X\x07\x00\x00\x00speciesqJX\x0c\x00\x00\x00homo_sapiensqKX\x07\x00\x00\x00projectqLX\x06\x00\x00\x00bfx123qMX\t\x00\x00\x00referenceqN}qO(X\x07\x00\x00\x00versionqPX\x04\x00\x00\x00hg38qQX\x04\x00\x00\x00pathqRXA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38qSuX\x07\x00\x00\x00mappingqT}qU(X\x05\x00\x00\x00indexqVh\x1fX\x07\x00\x00\x00programqWh\x1duuX\x06\x00\x00\x00outputqXcsnakemake.io\nOutputFiles\nqY)\x81qZX\x1f\x00\x00\x00gsnap/L25037_3PrimeCaoTest2.bamq[a}q\\(X\x03\x00\x00\x00bamq]h[h\n}q^h]K\x00N\x86q_subub.')
######## Original script #########
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
gsnap_cmd = "{binary} -D {snakemake.params.index_directory} -d {snakemake.params.index} -t $(expr ${{NSLOTS:-"+str(snakemake.threads)+"}} - 2) -A sam {extra}"

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
    gsnap_cmd = bwa_mem_cmd + " --read-group-id=\""+readgroup+"\" --read-group-name=\""+sample+"\" --read-group-library=\""+library+"\" "+("--gunzip" if re.search('\.gz$',snakemake.input[0]) else "") 

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
gsnap_cmd = "(module load apps/bwa"+modules["bwa"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true;" + gsnap_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+gsnap_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(gsnap_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))
