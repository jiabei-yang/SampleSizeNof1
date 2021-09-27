#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/output/RandEffPtsParams_%j.out
#SBATCH -e ../Data/error/RandEffPtsParams_%j.err
#SBATCH --job-name=RandEffPtsParams
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainPtsParams.R \
	--psblseq "Pairwise Randomization" \
	--corrstr "AR-1" \
	--alpha 0.05 \
	--beta 0.2 \
	--sigmaepsilonsq 4 \
	--rho 0.4 \
	--sigmamusq $1 \
	--sigmatausq $2 \
	--sigmamutau $3 \
	--mindiff 1 \
	--filename $4