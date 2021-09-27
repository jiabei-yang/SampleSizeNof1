#!/bin/bash

#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH --nodes 1
#SBATCH -o ../Data/output/MainPpobsParams_%j.out
#SBATCH -e ../Data/error/MainPpobsParams_%j.err
#SBATCH --job-name=MainPpobsParams
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jiabei_yang@brown.edu

Rscript MainPpobsParams.R \
	--psblseq "Pairwise Randomization" \
	--corrstr "AR-1" \
	--alpha 0.05 \
	--beta 0.2 \
	--sigmaepsilonsq 4 \
	--rho 0.4 \
	--sigmamusq 4 \
	--sigmatausq 1 \
	--sigmamutau 1 \
	--mindiff 1 \
	--filename "main/mainPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData"