#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/output/AlphBtPpobsParams_%j.out
#SBATCH -e ../Data/error/AlphBtPpobsParams_%j.err
#SBATCH --job-name=AlphBtPpobsParams
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainPpobsParams.R \
	--psblseq "Pairwise Randomization" \
	--corrstr "AR-1" \
	--alpha $1 \
	--beta $2 \
	--sigmaepsilonsq 4 \
	--rho 0.4 \
	--sigmamusq 4 \
	--sigmatausq 1 \
	--sigmamutau 1 \
	--mindiff 1 \
	--filename $3