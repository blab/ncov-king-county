# Phylodynamic Analysis 

## About the Analysis

To understand SARS-CoV-2 transmission within and between regions studied, we employ an approximate structured coalescent approach called MASCOT-GLM which also uses generalized linear models to inform estimates of effective population size and migration rates. More information can found in [Müller et al](https://academic.oup.com/ve/article/5/2/vez030/5549805?login=false)

## Running the Model

1. Make sure to run the script `KCparsimonyclusters.m` in `nextstrain_build/cluster_assignment/` beforehand 

2. navigate to `mascot_glm/data_for_tsv/` to combine all the cluster assignments from the four variant builds into one composite .tsv
2.1. first run `combining_metadata_and_cluster_assignments` followed by `combining_variant_files`

3. Produce the xmls for the three subsampling schemes using `mascotGLM.m` which uses `template_glm_from_multicoal_nomig_kc_clusters_combined_new.xml` as the template. Empirical predictors are already included in this template. The notebooks that produced and formatted these predictors can be found in `../analysis/mascot predictors/`

4. To run MASCOT GLM, while in the `mascot_glm/`folder in your terminal, type `java -jar Nab.jar xml/{xml_name}` for each separate subsampling scheme.

