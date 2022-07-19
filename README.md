# ARE, AIE and MIG estimators for ITR evaluation from observation data


This repository reproduces results from the paper *A Comprehensive Framework for the Evaluation of Individual Treatment Rules From Observational Data*.
The code implements the average rule effect (ARE), average implementation effect (AIE) and maximal implementation gain (MIG) estimators as proposed in the original paper.

### Reproducibility

- The **new_itr_situation** folder contains the following files.

 >`figure_2.R` implements the toy example given in the paper and repoduces Figure 2.
 >
 >`boot_func_new_itr.R` contains the bootstrap functions used for the new ITR situation application.
 >
 >`mimic_new_itr.R` reproduces Figure 5 for the new ITR situation application.

- The **partially_implemented_itr_situation** folder contains the following files.

>`algo1.R` implements the EM algorithm from the paper and returns ARE, AIE and MIG estimates. 
>ARE, AIE and MIG estimates along their bootstrap standard errors can be obtained in one line of code. 
>An example is given at the end of the file.
>
>`simulations.R` reproduces the simulations given in the paper.
>
>`plot_results.R` plots the results of the simulations and reproduces Figure 4 from the paper.
