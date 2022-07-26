#!/bin/bash -l
#SBATCH -J "GTs"
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 23-23:59:00
#SBATCH --mem=8000
#SBATCH -p standard
#SBATCH	-o test.r1.log

nextflow \
	run \
	/path/to/main.nf \
	-profile mfn
