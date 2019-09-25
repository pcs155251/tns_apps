#!/bin/bash
#PBS -N testProj
#PBS -P MST107027
#PBS -q cf160
#PBS -l select=1:ncpus=10
#PBS -M pcs155251@gmail.com
#PBS -m bea
#allow checking error and output at real time (home directory)
#PBS -k oe  

cd $PBS_O_WORKDIR

module load intel/2017_u4

./xxzdm.e input_xxzdm.rc
