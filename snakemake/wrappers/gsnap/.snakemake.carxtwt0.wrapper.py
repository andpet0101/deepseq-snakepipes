
######## Snakemake header ########
import sys; sys.path.insert(0, "/share/apps/python/3.4.2/lib/python3.4/site-packages/snakemake-3.12.0-py3.4.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x07\x00\x00\x00threadsq\x03K\x03X\x06\x00\x00\x00outputq\x04csnakemake.io\nOutputFiles\nq\x05)\x81q\x06X\x1f\x00\x00\x00gsnap/L25037_3PrimeCaoTest2.bamq\x07a}q\x08(X\x06\x00\x00\x00_namesq\t}q\nX\x03\x00\x00\x00bamq\x0bK\x00N\x86q\x0csh\x0bh\x07ubX\x06\x00\x00\x00paramsq\rcsnakemake.io\nParams\nq\x0e)\x81q\x0f(X\x00\x00\x00\x00q\x10X\x06\x00\x00\x00L25037q\x11X&\x00\x00\x00/projects/seq-work/user/pipeline/gmap/q\x12}q\x13X\x0c\x00\x00\x00${TMPDIR:-.}q\x14X\x15\x00\x00\x00L25037_3PrimeCaoTest2q\x15X\x04\x00\x00\x00hg38q\x16X\x0e\x00\x00\x003PrimeCaoTest2q\x17X\x05\x00\x00\x00gsnapq\x18X\n\x00\x00\x00coordinateq\x19e}q\x1a(X\x07\x00\x00\x00tmp_dirq\x1bh\x14X\x07\x00\x00\x00libraryq\x1ch\x11X\x06\x00\x00\x00binaryq\x1dh\x18X\x07\x00\x00\x00modulesq\x1eh\x13X\t\x00\x00\x00readgroupq\x1fh\x15X\x05\x00\x00\x00extraq h\x10h\t}q!(h\x1bK\x04N\x86q"h\x1cK\x01N\x86q#h\x1dK\x08N\x86q$h\x1eK\x03N\x86q%h K\x00N\x86q&h\x1fK\x05N\x86q\'X\x05\x00\x00\x00indexq(K\x06N\x86q)X\x07\x00\x00\x00sort_byq*K\tN\x86q+X\x0f\x00\x00\x00index_directoryq,K\x02N\x86q-X\x06\x00\x00\x00sampleq.K\x07N\x86q/uh(h\x16h*h\x19h,h\x12h.h\x17ubX\x03\x00\x00\x00logq0csnakemake.io\nLog\nq1)\x81q2X#\x00\x00\x00gsnap/log/L25037_3PrimeCaoTest2.logq3a}q4h\t}q5sbX\x06\x00\x00\x00configq6}q7(X\x07\x00\x00\x00mappingq8}q9(X\x07\x00\x00\x00programq:h\x18X\x05\x00\x00\x00indexq;h\x16uX\t\x00\x00\x00referenceq<}q=(X\x04\x00\x00\x00pathq>XA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38q?X\x07\x00\x00\x00versionq@X\x04\x00\x00\x00hg38qAuX\x07\x00\x00\x00projectqBX\x06\x00\x00\x00bfx123qCX\x07\x00\x00\x00speciesqDX\x0c\x00\x00\x00homo_sapiensqEuX\t\x00\x00\x00wildcardsqFcsnakemake.io\nWildcards\nqG)\x81qHX\x15\x00\x00\x00L25037_3PrimeCaoTest2qIa}qJ(h\t}qKX\x08\x00\x00\x00basenameqLK\x00N\x86qMsX\x08\x00\x00\x00basenameqNhIubX\x05\x00\x00\x00inputqOcsnakemake.io\nInputFiles\nqP)\x81qQX\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzqRa}qSh\t}qTsbX\x04\x00\x00\x00ruleqUX\t\x00\x00\x00run_gsnapqVX\t\x00\x00\x00resourcesqWcsnakemake.io\nResources\nqX)\x81qY(K\x01K\x03e}qZ(h\t}q[(X\x06\x00\x00\x00_nodesq\\K\x00N\x86q]X\x06\x00\x00\x00_coresq^K\x01N\x86q_uh\\K\x01h^K\x03ubub.')
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
    gsnap_cmd = bwa_mem_cmd + " --read-group-id=\""+readgroup+"\" --read-group-name=\""+sample+"\" --read-group-library=\""+library+"\" "+("--gunzip" if re.search('\.gz$',snakemake.input[0]) : "") 

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
