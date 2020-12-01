#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-1:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p guenette         # Partition to submit to
#SBATCH --mem=16000           # Memory pool for all cores (see also --mem-per-cpu)

PATH=$PATH:$HOME/bin
. /n/home02/jhaefner/miniconda/etc/profile.d/conda.sh
export ICTDIR=/n/home02/jhaefner/IC_modifications/IC
export ICDIR=$ICTDIR/invisible_cities
export PATH=$ICTDIR/bin:$PATH
export PYTHONPATH=$ICTDIR:$PYTHONPATH
conda activate IC-3.7-2020-06-16

city irene_unrolled /n/holystore01/LABS/guenette_lab/Users/jhaefner/irene_sliding_window_out/7487/config/irene_standard_FILENO.conf 
