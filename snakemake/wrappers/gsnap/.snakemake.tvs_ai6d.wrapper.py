
######## Snakemake header ########
import sys; sys.path.insert(0, "/share/apps/python/3.4.2/lib/python3.4/site-packages/snakemake-3.12.0-py3.4.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x04\x00\x00\x00ruleq\nX\t\x00\x00\x00run_gsnapq\x0bX\x07\x00\x00\x00threadsq\x0cK\x03X\t\x00\x00\x00resourcesq\rcsnakemake.io\nResources\nq\x0e)\x81q\x0f(K\x03K\x01e}q\x10(h\x08}q\x11(X\x06\x00\x00\x00_coresq\x12K\x00N\x86q\x13X\x06\x00\x00\x00_nodesq\x14K\x01N\x86q\x15uh\x12K\x03h\x14K\x01ubX\x06\x00\x00\x00paramsq\x16csnakemake.io\nParams\nq\x17)\x81q\x18(X&\x00\x00\x00/projects/seq-work/user/pipeline/gmap/q\x19}q\x1aX\x06\x00\x00\x00L25037q\x1bX\x05\x00\x00\x00gsnapq\x1cX\x00\x00\x00\x00q\x1dX\x0e\x00\x00\x003PrimeCaoTest2q\x1eX\n\x00\x00\x00coordinateq\x1fX\x04\x00\x00\x00hg38q X\x15\x00\x00\x00L25037_3PrimeCaoTest2q!X\x0c\x00\x00\x00${TMPDIR:-.}q"e}q#(X\x0f\x00\x00\x00index_directoryq$h\x19X\x07\x00\x00\x00modulesq%h\x1aX\x07\x00\x00\x00libraryq&h\x1bX\x06\x00\x00\x00binaryq\'h\x1cX\x05\x00\x00\x00extraq(h\x1dX\x06\x00\x00\x00sampleq)h\x1eX\x07\x00\x00\x00sort_byq*h\x1fh\x08}q+(h$K\x00N\x86q,h%K\x01N\x86q-h&K\x02N\x86q.h\'K\x03N\x86q/h(K\x04N\x86q0h)K\x05N\x86q1h*K\x06N\x86q2X\x05\x00\x00\x00indexq3K\x07N\x86q4X\x07\x00\x00\x00tmp_dirq5K\tN\x86q6X\t\x00\x00\x00readgroupq7K\x08N\x86q8uh3h h5h"h7h!ubX\x06\x00\x00\x00outputq9csnakemake.io\nOutputFiles\nq:)\x81q;X\x1f\x00\x00\x00gsnap/L25037_3PrimeCaoTest2.bamq<a}q=(h\x08}q>X\x03\x00\x00\x00bamq?K\x00N\x86q@sh?h<ubX\x06\x00\x00\x00configqA}qB(X\t\x00\x00\x00referenceqC}qD(X\x04\x00\x00\x00pathqEXA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38qFX\x07\x00\x00\x00versionqGX\x04\x00\x00\x00hg38qHuX\x07\x00\x00\x00speciesqIX\x0c\x00\x00\x00homo_sapiensqJX\x07\x00\x00\x00projectqKX\x06\x00\x00\x00bfx123qLX\x07\x00\x00\x00mappingqM}qN(X\x07\x00\x00\x00programqOh\x1cX\x05\x00\x00\x00indexqPh uuX\t\x00\x00\x00wildcardsqQcsnakemake.io\nWildcards\nqR)\x81qSX\x15\x00\x00\x00L25037_3PrimeCaoTest2qTa}qU(h\x08}qVX\x08\x00\x00\x00basenameqWK\x00N\x86qXsX\x08\x00\x00\x00basenameqYhTubX\x03\x00\x00\x00logqZcsnakemake.io\nLog\nq[)\x81q\\X#\x00\x00\x00gsnap/log/L25037_3PrimeCaoTest2.logq]a}q^h\x08}q_sbub.')
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
