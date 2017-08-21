#!/bin/sh
# properties = {properties}

# the job runs on
hostname

# do not do core dumps (you cannot imagine how fast 5Tb of data are consumed)
ulimit -c 0

# set maximum number of open files to hard limit (ulimit -Hn
ulimit -n 4096

# fix group write
umask g+w

# load snakemake module
module load apps/snakemake

# run job
{exec_job}
