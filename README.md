# Refining capture-recapture methods to estimate case counts in a finite population setting

This github repository contains the code for replicating the simulation studies and real-data analysis in the paper, "Refining capture-recapture methods to estimate case counts in a finite population setting."  `01_functions.R` contains helper functions for running the simulation study, including functions for calculating:


- Chapman estimator:  **N̂<sub>Chap</sub>**
- Random sample estimator:  **N̂<sub>RS</sub>**
- 5-cell estimator:  **N̂<sub>MLE</sub>**
- FPC1 and FPC2 Variance Estimators
-Wald Intervals
-Adapted Bayesian Credible Intervals

All functions assume the input data comes in the form **(n₁₅, n₂, n₄, n₆, n₃₇)<sup>T</sup>**.

`02_simulations.R` contains functions that implement the simulation studies in the paper and automatically convert output to LaTeX.  

`CRISP Data Analysis.RMD` is an RMarkdown file that replicates the real-data analysis on the CRISP data.

`sim_results1b_arXiv.RData` contains the simulation results for Table 4 in the mansucript.

`sim_results2b_arXiv.RData` contains the simulation results for Table 6 in the mansucript.

`sim_resultsA3b_arXiv.RData` contains the simulation results for Table S1 in the mansucript.

`sim_resultsA4b_arXiv.RData` contains the simulation results for Table S2 in the mansucript.

`sim_results_CRISPb_arXiv.RData` contains the simulation results for Table S3 in the mansucript.



🧑‍💻 Author

Michael Doerfler
Ph.D. Candidate, Biostatistics, Emory University
📧 michael.doerfler@emory.edu
