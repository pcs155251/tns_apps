#!/bin/bash
#PBS -N conv14
#PBS -P MST107027
#PBS -q dc20190021
#PBS -l select=1:ncpus=40:mem=180gb
#PBS -M pcs155251@gmail.com
#PBS -m bea
#allow checking error and output at real time (home directory)
#PBS -W sandbox=PRIVATE

cd $PBS_O_WORKDIR

module load intel/2017_u4

./en.e input.rc  -dimDstart 14 -dimDend 15 -saveFolder ctm_convP/ -enviFolder ctm_convP/

