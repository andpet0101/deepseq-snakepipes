"""Snakemake wrapper for samtools index."""


from os.path import basename
from snakemake.shell import shell

# optional params
binary = snakemake.params.get('binary', 'samtools')
extra = snakemake.params.get('extra', '')
modules = snakemake.params.get("modules",dict())

index_cmd = "{binary} index {snakemake.input[0]} {snakemake.output[0]}"
	
# module load command
if "samtools" in modules and modules["samtools"]:
    modules["samtools"] = "/" + modules["samtools"]
else:
    modules["samtools"] = ""
index_cmd = "(module load apps/samtools"+modules["samtools"]+" || true;module list -l || true;" + index_cmd + ")"

# first print shell command, then run
shell("(echo NODE: $(hostname);echo STARTED AT: $(date); echo COMMAND: '"+index_cmd+"'; echo) "+snakemake.log_fmt_shell(append=False))
shell(index_cmd + ' '+snakemake.log_fmt_shell(append=True))
shell("(echo;echo FINISHED AT: $(date)) "+snakemake.log_fmt_shell(append=True))
