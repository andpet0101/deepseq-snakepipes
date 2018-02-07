
######## Snakemake header ########
import sys; sys.path.insert(0, "/share/apps/python/3.4.2/lib/python3.4/site-packages/snakemake-3.12.0-py3.4.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x03\x00\x00\x00logq\x03csnakemake.io\nLog\nq\x04)\x81q\x05X#\x00\x00\x00gsnap/log/L25037_3PrimeCaoTest2.logq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x07\x00\x00\x00threadsq\nK\x03X\x06\x00\x00\x00paramsq\x0bcsnakemake.io\nParams\nq\x0c)\x81q\r(X\x05\x00\x00\x00gsnapq\x0eX&\x00\x00\x00/projects/seq-work/user/pipeline/gmap/q\x0fX\x06\x00\x00\x00L25037q\x10X\x04\x00\x00\x00hg38q\x11X\x0e\x00\x00\x003PrimeCaoTest2q\x12X\n\x00\x00\x00coordinateq\x13X\x0c\x00\x00\x00${TMPDIR:-.}q\x14X\x00\x00\x00\x00q\x15X\x15\x00\x00\x00L25037_3PrimeCaoTest2q\x16}q\x17e}q\x18(X\x0f\x00\x00\x00index_directoryq\x19h\x0fX\x07\x00\x00\x00libraryq\x1ah\x10X\t\x00\x00\x00readgroupq\x1bh\x16X\x06\x00\x00\x00sampleq\x1ch\x12h\x08}q\x1d(h\x19K\x01N\x86q\x1eX\x06\x00\x00\x00binaryq\x1fK\x00N\x86q h\x1bK\x08N\x86q!X\x07\x00\x00\x00modulesq"K\tN\x86q#h\x1cK\x04N\x86q$X\x07\x00\x00\x00sort_byq%K\x05N\x86q&X\x07\x00\x00\x00tmp_dirq\'K\x06N\x86q(X\x05\x00\x00\x00extraq)K\x07N\x86q*X\x05\x00\x00\x00indexq+K\x03N\x86q,h\x1aK\x02N\x86q-uh%h\x13h\'h\x14h)h\x15h+h\x11h"h\x17h\x1fh\x0eubX\x05\x00\x00\x00inputq.csnakemake.io\nInputFiles\nq/)\x81q0X\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzq1a}q2h\x08}q3sbX\x04\x00\x00\x00ruleq4X\t\x00\x00\x00run_gsnapq5X\x06\x00\x00\x00outputq6csnakemake.io\nOutputFiles\nq7)\x81q8X\x1f\x00\x00\x00gsnap/L25037_3PrimeCaoTest2.bamq9a}q:(h\x08}q;X\x03\x00\x00\x00bamq<K\x00N\x86q=sh<h9ubX\t\x00\x00\x00resourcesq>csnakemake.io\nResources\nq?)\x81q@(K\x03K\x01e}qA(X\x06\x00\x00\x00_coresqBK\x03h\x08}qC(hBK\x00N\x86qDX\x06\x00\x00\x00_nodesqEK\x01N\x86qFuhEK\x01ubX\x06\x00\x00\x00configqG}qH(X\x07\x00\x00\x00speciesqIX\x0c\x00\x00\x00homo_sapiensqJX\x07\x00\x00\x00projectqKX\x06\x00\x00\x00bfx123qLX\x07\x00\x00\x00mappingqM}qN(X\x07\x00\x00\x00programqOh\x0eX\x05\x00\x00\x00indexqPh\x11uX\t\x00\x00\x00referenceqQ}qR(X\x07\x00\x00\x00versionqSX\x04\x00\x00\x00hg38qTX\x04\x00\x00\x00pathqUXA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38qVuuX\t\x00\x00\x00wildcardsqWcsnakemake.io\nWildcards\nqX)\x81qYX\x15\x00\x00\x00L25037_3PrimeCaoTest2qZa}q[(h\x08}q\\X\x08\x00\x00\x00basenameq]K\x00N\x86q^sX\x08\x00\x00\x00basenameq_hZubub.')
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
gsnap_cmd = "(module load apps/bwa"+modules["bwa"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true;" + gsnap_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+gsnap_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(gsnap_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))
