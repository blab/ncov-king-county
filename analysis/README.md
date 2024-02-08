# Processing MASCOT-GLM results 

## About the Analysis

Here we work on taking the output from MASCOT-GLM (which can be found in `../mascot_glm/results/`) and analyzing the logs and posterior set of trees. 

## Organization

1. `data-files/` contains both the input and output files for all the posterior processing analyses. Here you can find MCC trees that are then fed into `../figures/` to make the final figures

2. `mascot predictors/` contains all the notebooks to make clean and format the empirical predictors used in our MASCOT-GLM analysis. All relevant predictors have already been included in the template xml under `../mascot_glm/` so no further modification is needed to reproduce the results. 

3. `post-process` contains all the notebooks that extract and analyze information from the posterior set of trees from MASCOT-GLM. `post-process/beast/` contains all the notebooks directly related to MASCOT-GLM files. All other notebooks analyze mobility patterns and demographic information 

### Final figures found in the manuscript can be reproducted in `../figures/`

