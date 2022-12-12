# Local-Scale phylodynamics reveal differential community impact of SARS-CoV-2 in metropolitan US county
Miguel I. Paredes1,2*, Amanda C. Perofsky3,4, Lauren Frisbie5, Louise H. Moncla6, Pavitra Roychoudhury2,7, Hong Xie7, Shah A. Mohamed Bakhash7, Kevin Kong7, Isabel Arnould7, Tien V. Nguyen7, Seffir T. Wendm7, Pooneh Hajian7, Sean Ellis7, Patrick C. Mathias7, Alexander Greninger2,7, Lea M. Starita3,8, Chris D. Frazar8, Erica Ryke8, Weizhi Zhong3, Luis Gamboa3, Machiko Threlkeld8, Jover Lee2, Jeremy Stone3, Evan McDermot3, Melissa Truong8, Jay Shendure3,8,9, Hanna Oltean5, Cécile Viboud4, Helen Chu10, Nicola F. Müller2† , Trevor Bedford1,2,3,8,9†


Affiliations
1 Department of Epidemiology, University of Washington, Seattle, WA, USA
2 Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, Washington, USA
3 Brotman Baty Institute for Precision Medicine, University of Washington, Seattle, WA USA
4 Fogarty International Center, National Institutes of Health, Bethesda, MD, USA
5 Washington State Department of Health, Shoreline, WA USA
6 The University of Pennsylvania, Department of Pathobiology, Philadelphia, PA
7 Department of Laboratory Medicine and Pathology, University of Washington, Seattle, WA, USA
8 Department of Genome Sciences, University of Washington, Seattle, WA, USA
9 Howard Hughes Medical Institute, Seattle, WA, USA
10 Department of Medicine, Division of Allergy and Infectious Diseases, University of Washington, Seattle, WA
†These authors jointly supervised this work.

## Abstract

SARS-CoV-2 transmission is largely driven by heterogeneous dynamics at a local scale, leaving local health departments to design interventions with limited information. We analyzed SARS-CoV-2 genomes sampled between February 2020 and March 2022 jointly with epidemiological and cell phone mobility data to investigate fine scale spatiotemporal SARS-CoV-2 transmission dynamics in King County, Washington, a diverse, metropolitan US county. We applied an approximate structured coalescent approach to model transmission within and between North King County and South King County alongside the rate of outside introductions into the county. Our phylodynamic analyses reveal that following stay-at-home orders, the epidemic trajectories of North and South King County began to diverge. We find that South King County consistently had more reported and estimated cases, COVID-19 hospitalizations, and longer persistence of local viral transmission when compared to North King County, where viral importations from outside drove a larger proportion of new cases. Using mobility and demographic data, we also find that South King County experienced a more modest and less sustained reduction in mobility following stay-at-home orders than North King County, while also bearing more socioeconomic inequities that might contribute to a disproportionate burden of SARS-CoV-2 transmission. Overall, our findings suggest a role for local-scale phylodynamics in understanding the heterogeneous transmission landscape.

## Organization

This repository contains the analytic code needed to reproduce the results from the above paper. To start, begin with the folder `nextstrain_build` to run the maximum likelihood analysis and create the temporally resolved phylogeny of SARS-CoV-2 in King County.

Then `nextstrain_build/cluster_assignment/` contains the script `KCparsimonyclusters.m` needed to assign local outbreak clusters needed to conduct the MASCOT-GLM analysis

After cluster assignment, navigate to `mascot_glm/data_for_tsv/` to combine all the cluster assignments from the four variant builds into one composite .tsv

`mascot_glm/mascotGLM.m` contains the code to produce the xmls for the three subsampling schemes that are inputed into MASCOT GLM

To run MASCOT GLM, while in the `mascot_glm/`folder in your terminal, type `java -jar Nab.jar xml/{xml_name}` for each separate subsampling scheme.

To analyze the posterior set of trees and log files produced by MASCOT-GLM, go to `analysis/post-process/` for all scripts

Finally, to produce the manuscript figures, run the notebooks found in `analysis/figures/`

All sequence accession IDs and associated North/South King County divisions used in this manuscript and provided by the Washington State Department of Health can be found in `final_acknowledgements_gisaid.csv`
