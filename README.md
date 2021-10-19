# Sample size calculations for n-of-1 trials

Code for paper [Yang, J., Steingrimsson, J. A., & Schmid, C. H. Sample size calculations for n-of-1 trials.](https://arxiv.org/abs/2110.08970)

* `main/`: code for illustrations in the main manuscript and appendices.
  * Simulations:
    + `MainPpobsParams.R`: run the illustration in Figure 2 (left) and appendices which require appropriate input from the shell scripts.
    + `MainPtsParams.R`: run the illustration in Figure 2 (right) and appendices which require appropriate input from the shell scripts.
    + `MainIndSe.R`: calculate individual standard error and produce Figure 3.
  * `Run*.sh`: shell scripts to run associated R scripts on `slurm`.

* `Functions/`: functions to calculate the power for estimating the population average treatment effect and standard errors for individual-specific treatment effects.

* `Data/`: results from running the code in under `main/`.

* `Results/`: code to generate figures and tables in the paper.

* `Shiny/`: code for the Shiny app.
  * `TestDesgin.csv`: example input file if one wants to specify sequences of interest.