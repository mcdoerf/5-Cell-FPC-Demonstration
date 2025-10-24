# 5-Cell-FPC-Demonstration
This github repository contains the code for replicating the simulation studies and real-data analysis in the paper, "Refining capture-recapture methods to estimate case counts in a finite population setting."  `Script1.R` contains helper functions for running the simulation study, including functions for calculating:


- Chapman estimator:  **N̂<sub>Chap</sub>**
- Random sample estimator:  **N̂<sub>RS</sub>**
- 5-cell estimator:  **N̂<sub>MLE</sub>**

All functions assume the input data comes in the form **(n₁₅, n₂, n₄, n₆, n₃₇)<sup>T</sup>**.
