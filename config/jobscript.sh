#!/bin/sh
# properties = {properties}

# works with qsub --notify
#trap "echo 'Caught SIGUSR1/SIGURS2 - going down';kill -- -$$;kill -9 -- -$$" SIGUSR1 SIGUSR2

hostname
if [ ! -d "/project/seq-work" ] && [ ! -d "/home/sequencing" ]
then
  echo "Error: Either /project/seq-work or /home/sequencing is not mounted"
  exit 1
fi

module load apps/snakemake
{exec_job}
