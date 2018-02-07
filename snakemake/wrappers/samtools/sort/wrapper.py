"""Snakemake wrapper for samtools sort."""


from os.path import basename
from snakemake.shell import shell
import re

# optional params
binary = snakemake.params.get('binary', 'samtools')
extra = snakemake.params.get('extra', '')
tmp_dir = snakemake.params.get("tmp_dir", ".")
sort_by = snakemake.params.get("sort_by", "coordinate")
default_threads = snakemake.threads
modules = snakemake.params.get("modules",dict())

base_name = re.sub('\.[^\.]+$','',basename(snakemake.input[0]))

sort_cmd = "{binary} sort "+("-n" if sort_by=="queryname" else "")+" -T {tmp_dir}/"+base_name+" {extra} -o {snakemake.output[0]} {snakemake.input[0]}"
	
# module load command
if "samtools" in modules and modules["samtools"]:
    modules["samtools"] = "/" + modules["samtools"]
else:
    modules["samtools"] = ""
sort_cmd = "(module load apps/samtools"+modules["samtools"]+" || true;module list -l || true;" + sort_cmd + ")"

# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+sort_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(sort_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))

