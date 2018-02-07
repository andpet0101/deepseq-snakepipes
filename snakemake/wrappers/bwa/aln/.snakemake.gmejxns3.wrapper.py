
######## Snakemake header ########
import sys; sys.path.insert(0, "/home/sequencing/andpetzo/tools/miniconda3/lib/python3.5/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x06\x00\x00\x00outputq\x03csnakemake.io\nOutputFiles\nq\x04)\x81q\x05X \x00\x00\x00bwaaln/L25037_3PrimeCaoTest2.bamq\x06a}q\x07(X\x03\x00\x00\x00bamq\x08h\x06X\x06\x00\x00\x00_namesq\t}q\nh\x08K\x00N\x86q\x0bsubX\x05\x00\x00\x00inputq\x0ccsnakemake.io\nInputFiles\nq\r)\x81q\x0eX\'\x00\x00\x00fastq/L25037_3PrimeCaoTest2_R1.fastq.gzq\x0fa}q\x10h\t}q\x11sbX\x07\x00\x00\x00threadsq\x12K\x01X\x03\x00\x00\x00logq\x13csnakemake.io\nLog\nq\x14)\x81q\x15X$\x00\x00\x00bwaaln/log/L25037_3PrimeCaoTest2.logq\x16a}q\x17h\t}q\x18sbX\t\x00\x00\x00resourcesq\x19csnakemake.io\nResources\nq\x1a)\x81q\x1b(K\x01K\x01e}q\x1c(X\x06\x00\x00\x00_nodesq\x1dK\x01X\x06\x00\x00\x00_coresq\x1eK\x01h\t}q\x1f(h\x1dK\x00N\x86q h\x1eK\x01N\x86q!uubX\t\x00\x00\x00wildcardsq"csnakemake.io\nWildcards\nq#)\x81q$X\x15\x00\x00\x00L25037_3PrimeCaoTest2q%a}q&(X\x08\x00\x00\x00basenameq\'h%h\t}q(X\x08\x00\x00\x00basenameq)K\x00N\x86q*subX\x06\x00\x00\x00configq+}q,(X\t\x00\x00\x00referenceq-}q.(X\x07\x00\x00\x00versionq/X\x04\x00\x00\x00hg38q0X\x04\x00\x00\x00pathq1XA\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38q2uX\x07\x00\x00\x00mappingq3}q4(X\x05\x00\x00\x00indexq5X\x04\x00\x00\x00hg38q6X\x07\x00\x00\x00programq7X\x06\x00\x00\x00bwaalnq8uX\x07\x00\x00\x00speciesq9X\x0c\x00\x00\x00homo_sapiensq:X\x07\x00\x00\x00projectq;X\x06\x00\x00\x00bfx123q<uX\x06\x00\x00\x00paramsq=csnakemake.io\nParams\nq>)\x81q?(X\x00\x00\x00\x00q@X\x15\x00\x00\x00L25037_3PrimeCaoTest2qAX%\x00\x00\x00/projects/seq-work/user/pipeline/bwa/qBX\n\x00\x00\x00coordinateqCh6X\x0e\x00\x00\x003PrimeCaoTest2qDh@\x89X\x06\x00\x00\x00L25037qEX\x03\x00\x00\x00bwaqFXD\x00\x00\x00/projects/seq-work/user/pipeline/reference/homo_sapiens/hg38/hg38.faqGX\x0c\x00\x00\x00${TMPDIR:-.}qH}qIe}qJ(X\x0c\x00\x00\x00extra_samspeqKh@X\n\x00\x00\x00paired_endqL\x89X\x0f\x00\x00\x00index_directoryqMhBX\x07\x00\x00\x00sort_byqNhCX\x05\x00\x00\x00indexqOh6X\x06\x00\x00\x00sampleqPhDX\t\x00\x00\x00extra_alnqQh@X\t\x00\x00\x00readgroupqRhAX\x07\x00\x00\x00libraryqShEX\x06\x00\x00\x00binaryqThFX\t\x00\x00\x00referenceqUhGX\x07\x00\x00\x00tmp_dirqVhHh\t}qW(hKK\x00N\x86qXhLK\x07N\x86qYhMK\x02N\x86qZhNK\x03N\x86q[hOK\x04N\x86q\\hPK\x05N\x86q]hQK\x06N\x86q^hRK\x01N\x86q_hSK\x08N\x86q`hTK\tN\x86qahUK\nN\x86qbhVK\x0bN\x86qcX\x07\x00\x00\x00modulesqdK\x0cN\x86qeuhdhIubX\x04\x00\x00\x00ruleqfX\x0b\x00\x00\x00run_bwa_alnqgub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
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
bwa_aln_cmd = "{binary} aln -t $(expr ${{NSLOTS:-{default_threads}}} - 2) {extra_aln}"

# bam input
if re.search('\.bam$',snakemake.input[0]):
	if not reference:
		raise ValueError("Need a reference path given by params.reference when using BAM files as input.")
	if paired_end:
		bwa_aln_cmd = "{binary} sampe {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}"
		bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -b1 -t $(expr ${{NSLOTS:-{default_threads}}} / 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]})"
		bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -b2 -t $(expr ${{NSLOTS:-{default_threads}}} / 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]})"
		bwa_aln_cmd = bwa_aln_cmd + " {snakemake.input[0]} {snakemake.input[0]}"
	else:
		bwa_aln_cmd = "{binary} samse {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}" 
		bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -b0 -t ${{NSLOTS:-{default_threads}}} {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]}) {snakemake.input[0]}"

	bwa_aln_cmd = bwa_aln_cmd + " | picardtools MergeBamAlignment UNMAPPED={snakemake.input[0]} ALIGNED=/dev/stdin O={snakemake.output[0]} R="+reference+" MAX_GAPS=-1 ALIGNER_PROPER_PAIR_FLAGS=true CLIP_OVERLAPPING_READS=false ADD_MATE_CIGAR=true TMP_DIR={tmp_dir}"
# fastq input
else:
	if paired_end:
		bwa_aln_cmd = "{binary} sampe -r \"@RG\\tID:"+readgroup+"\\tSM:"+sample+"\\tLB:"+library+"\\tPL:ILLUMINA\\tCN:DeepSeqDresden\" {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}"
		bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -t $(expr ${{NSLOTS:-{default_threads}}} / 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]})"
		bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -t $(expr ${{NSLOTS:-{default_threads}}} / 2) {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[1]})"
		bwa_aln_cmd = bwa_aln_cmd + " {snakemake.input[0]} {snakemake.input[1]}"
	else:
		bwa_aln_cmd = "{binary} samse -r \"@RG\\tID:"+readgroup+"\\tSM:"+sample+"\\tLB:"+library+"\\tPL:ILLUMINA\\tCN:DeepSeqDresden\" {extra_samspe} {snakemake.params.index_directory}/{snakemake.params.index}" 
		bwa_aln_cmd = bwa_aln_cmd + " <({binary} aln -t ${{NSLOTS:-{default_threads}}} {extra_aln} {snakemake.params.index_directory}/{snakemake.params.index} {snakemake.input[0]}) {snakemake.input[0]}"
	
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
bwa_aln_cmd = "(module load apps/bwa"+modules["bwa"]+" apps/samtools"+modules["samtools"]+" apps/picardtools/"+modules["picardtools"]+" || true;module list -l || true; trap 'kill 0' SIGINT;" + bwa_aln_cmd + ")"


# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+bwa_aln_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(bwa_aln_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))


