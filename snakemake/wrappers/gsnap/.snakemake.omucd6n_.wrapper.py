
######## Snakemake header ########
import sys; sys.path.insert(0, "/share/apps/python/3.4.2/lib/python3.4/site-packages/snakemake-3.12.0-py3.4.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\t\x00\x00\x00resourcesq\x03csnakemake.io\nResources\nq\x04)\x81q\x05(K\x03K\x01e}q\x06(X\x06\x00\x00\x00_namesq\x07}q\x08(X\x06\x00\x00\x00_nodesq\tK\x01N\x86q\nX\x06\x00\x00\x00_coresq\x0bK\x00N\x86q\x0cuh\tK\x01h\x0bK\x03ubX\x06\x00\x00\x00configq\r}q\x0e(X\x07\x00\x00\x00speciesq\x0fX\x0c\x00\x00\x00homo_sapiensq\x10X\x07\x00\x00\x00projectq\x11X\x06\x00\x00\x00bfx123q\x12X\x07\x00\x00\x00mappingq\x13}q\x14(X\x07\x00\x00\x00programq\x15X\x05\x00\x00\x00gsnapq\x16X\x05\x00\x00\x00indexq\x17X\x04\x00\x00\x00hg38q\x18uX\t\x00\x00\x00referenceq\x19}q\x1a(X\x07\x00\x00\x00versionq\x1bX\x04\x00\x00\x00hg38q\x1cX\x04\x00\x00\x00pathq\x1dXA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38q\x1euuX\t\x00\x00\x00wildcardsq\x1fcsnakemake.io\nWildcards\nq )\x81q!X\x15\x00\x00\x00L25037_3PrimeCaoTest2q"a}q#(h\x07}q$X\x08\x00\x00\x00basenameq%K\x00N\x86q&sX\x08\x00\x00\x00basenameq\'h"ubX\x03\x00\x00\x00logq(csnakemake.io\nLog\nq))\x81q*X#\x00\x00\x00gsnap/log/L25037_3PrimeCaoTest2.logq+a}q,h\x07}q-sbX\x06\x00\x00\x00paramsq.csnakemake.io\nParams\nq/)\x81q0(X\x00\x00\x00\x00q1X\x0e\x00\x00\x003PrimeCaoTest2q2h\x18X\x06\x00\x00\x00L25037q3X\x15\x00\x00\x00L25037_3PrimeCaoTest2q4X&\x00\x00\x00/projects/seq-work/user/pipeline/gmap/q5}q6X\x0c\x00\x00\x00${TMPDIR:-.}q7X\n\x00\x00\x00coordinateq8h\x16e}q9(X\x07\x00\x00\x00modulesq:h6X\x06\x00\x00\x00sampleq;h2X\t\x00\x00\x00readgroupq<h4h\x07}q=(h:K\x06N\x86q>h;K\x01N\x86q?X\x0f\x00\x00\x00index_directoryq@K\x05N\x86qAX\x07\x00\x00\x00libraryqBK\x03N\x86qCX\x06\x00\x00\x00binaryqDK\tN\x86qEh<K\x04N\x86qFX\x05\x00\x00\x00indexqGK\x02N\x86qHX\x05\x00\x00\x00extraqIK\x00N\x86qJX\x07\x00\x00\x00sort_byqKK\x08N\x86qLX\x07\x00\x00\x00tmp_dirqMK\x07N\x86qNuhBh3hIh1h@h5hGh\x18hMh7hKh8hDh\x16ubX\x06\x00\x00\x00outputqOcsnakemake.io\nOutputFiles\nqP)\x81qQX\x1f\x00\x00\x00gsnap/L25037_3PrimeCaoTest2.bamqRa}qS(h\x07}qTX\x03\x00\x00\x00bamqUK\x00N\x86qVshUhRubX\x07\x00\x00\x00threadsqWK\x03X\x04\x00\x00\x00ruleqXX\t\x00\x00\x00run_gsnapqYX\x05\x00\x00\x00inputqZcsnakemake.io\nInputFiles\nq[)\x81q\\X\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzq]a}q^h\x07}q_sbub.')
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
gsnap_cmd = "(module load apps/gmap"+modules["gmap"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true;" + gsnap_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+gsnap_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(gsnap_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))
