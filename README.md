# RHDSI: A Novel Dimensionality Reduction Based Algorithm on High Dimensional Feature Selection with Interactions 
The algorithm provides a unique approach to select features in high dimensional settings. It reduces the dimensions of feature before performing three stage feature selection. Currently, algorithm only handles continuous features and outcome. This tutorial explains the process of using the algorithm.

## Install Libraries
This algorithm uses multiple libraries. Those libraries need to be installed before this algorithm could be used.

```
# Versions used in the current algorithm are:
# memoise (>= 1.1.0), dplyr (>= 0.8.4), plyr (>= 1.8.6), stringr (>= 1.4.0), stats (>= 3.6.2), scales (>= 1.1.0), 
# caret (>= 6.0-86), magrittr (>= 1.5), janitor (>= 1.2.1), naturalsort (>= 0.1.3), plsdepot (>= 0.1.17), 
# spatstat.utils (>= 1.17-0), GA (>= 3.2), future (>= 1.17.0), future.apply (>= 1.6.0), purrr (>= 0.3.3), 
# MLmetrics (>= 1.1.1), reshape2 (>= 1.4.4), ggplot2 (>= 3.3.2), readxl (>= 1.3.1)

install.packages("pacman")
library(pacman)
pacman::p_load(memoise, dplyr, plyr, stringr, stats, scales, caret, magrittr, janitor, naturalsort, plsdepot, 
               spatstat.utils, GA, future, future.apply, purrr, MLmetrics, reshape2, ggplot2, readxl) 
```
## Perform hyperparameter optimization on simulated datset
This step will generate a simulated dataset and select the best RHDSI hyperparameters for the dataset. In the current example, dataset contains 15 marginal features with 19 target features (10 marginal features and 9 interaction features). All 15 marginal features are correlated with pariwise correlation of 0.5. The dataset has only 100 samples in the training set and 500 samples in the test set.

```
best_para = cv_hyperpara(datatype=c("simulated"), param = list(seed =1, test_train_ratio=NA, varnum = 15, 
                         setting= "Correlation", main_var=10, var_effect=c(0.5, -0.5), correlation_var=15, 
                         correlation_val=5, high_dim=T, train_sample=100, var = "Mar"), 
                         seeder=1, effectsize=13, optimtype = "GA", cv = T, predict_method = "reg")
best_para = as.numeric(best_para)
print(best_para)
```
```param``` takes the various arguments needed to generate simulated dataset. The four hyperparameters namely number of features used in sample (q) [```k```], number of latent factors (Lf) [```ncomp```], confidence level for selecting features in Stage I (CL) [```summary_ci```] and number of features to remove from the target cluster (Af) in Stage II [```varmax```]. 

## Perform RHDSI
Use the best values of hyperparameters to run RHDSI algorithm to get final results.

```
res= fit_function(datatype=c("simulated"), param = list(seed =1, varnum = 15, setting= "Correlation", 
                  main_var=10, var_effect=c(0.5, -0.5), correlation_var=15, correlation_val=5, high_dim=T,
                  train_sample=100, var = "Mar"), seeder=i, ncomp=best_para[2], k=best_para[1], effectsize=13, 
                  summary_ci = best_para[3], varmax=best_para[4], predict_method = "reg")
res
```

The output ```res``` is a list. The first element in the list provide the performance of the algorithm in the test dataset and second element provides the list of features selected by the algorithm.
