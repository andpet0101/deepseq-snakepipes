
######## Snakemake header ########
import sys; sys.path.insert(0, "/share/apps/python/3.4.2/lib/python3.4/site-packages/snakemake-3.12.0-py3.4.egg"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x04\x00\x00\x00ruleq\x03X\t\x00\x00\x00run_gsnapq\x04X\x05\x00\x00\x00inputq\x05csnakemake.io\nInputFiles\nq\x06)\x81q\x07X\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzq\x08a}q\tX\x06\x00\x00\x00_namesq\n}q\x0bsbX\t\x00\x00\x00wildcardsq\x0ccsnakemake.io\nWildcards\nq\r)\x81q\x0eX\x15\x00\x00\x00L25037_3PrimeCaoTest2q\x0fa}q\x10(X\x08\x00\x00\x00basenameq\x11h\x0fh\n}q\x12X\x08\x00\x00\x00basenameq\x13K\x00N\x86q\x14subX\x03\x00\x00\x00logq\x15csnakemake.io\nLog\nq\x16)\x81q\x17X#\x00\x00\x00gsnap/log/L25037_3PrimeCaoTest2.logq\x18a}q\x19h\n}q\x1asbX\x06\x00\x00\x00configq\x1b}q\x1c(X\t\x00\x00\x00referenceq\x1d}q\x1e(X\x07\x00\x00\x00versionq\x1fX\x04\x00\x00\x00hg38q X\x04\x00\x00\x00pathq!XA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38q"uX\x07\x00\x00\x00projectq#X\x06\x00\x00\x00bfx123q$X\x07\x00\x00\x00speciesq%X\x0c\x00\x00\x00homo_sapiensq&X\x07\x00\x00\x00mappingq\'}q((X\x07\x00\x00\x00programq)X\x05\x00\x00\x00gsnapq*X\x05\x00\x00\x00indexq+X\x04\x00\x00\x00hg38q,uuX\t\x00\x00\x00resourcesq-csnakemake.io\nResources\nq.)\x81q/(K\x03K\x01e}q0(X\x06\x00\x00\x00_coresq1K\x03X\x06\x00\x00\x00_nodesq2K\x01h\n}q3(h1K\x00N\x86q4h2K\x01N\x86q5uubX\x06\x00\x00\x00paramsq6csnakemake.io\nParams\nq7)\x81q8(h,X\x06\x00\x00\x00L25037q9X\x0e\x00\x00\x003PrimeCaoTest2q:X&\x00\x00\x00/projects/seq-work/user/pipeline/gmap/q;X\n\x00\x00\x00coordinateq<X\x00\x00\x00\x00q=X\x0c\x00\x00\x00${TMPDIR:-.}q>}q?X\x15\x00\x00\x00L25037_3PrimeCaoTest2q@h*e}qA(X\x07\x00\x00\x00modulesqBh?h\n}qC(X\x06\x00\x00\x00binaryqDK\tN\x86qEX\x05\x00\x00\x00extraqFK\x05N\x86qGX\x07\x00\x00\x00libraryqHK\x01N\x86qIX\x06\x00\x00\x00sampleqJK\x02N\x86qKX\x0f\x00\x00\x00index_directoryqLK\x03N\x86qMX\x07\x00\x00\x00sort_byqNK\x04N\x86qOX\x05\x00\x00\x00indexqPK\x00N\x86qQX\x07\x00\x00\x00tmp_dirqRK\x06N\x86qSX\t\x00\x00\x00readgroupqTK\x08N\x86qUhBK\x07N\x86qVuhDh*hHh9hJh:hLh;hNh<hPh,hRh>hTh@hFh=ubX\x07\x00\x00\x00threadsqWK\x03X\x06\x00\x00\x00outputqXcsnakemake.io\nOutputFiles\nqY)\x81qZX\x1f\x00\x00\x00gsnap/L25037_3PrimeCaoTest2.bamq[a}q\\(X\x03\x00\x00\x00bamq]h[h\n}q^h]K\x00N\x86q_subub.')
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
