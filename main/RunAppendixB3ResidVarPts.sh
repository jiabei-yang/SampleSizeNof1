#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/output/ResidVarPtsParams_%j.out
#SBATCH -e ../Data/error/ResidVarPtsParams_%j.err
#SBATCH --job-name=ResidVarPtsParams
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainPtsParams.R \
	--psblseq "Pairwise Randomization" \
	--corrstr $1 \
	--alpha 0.05 \
	--beta 0.2 \
	--sigmaepsilonsq $2 \
	--rho $3 \
	--sigmamusq 4 \
	--sigmatausq 1 \
	--sigmamutau 1 \
	--mindiff 1 \
	--filename $4