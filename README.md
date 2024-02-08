# Local-Scale phylodynamics reveal differential community impact of SARS-CoV-2 in metropolitan US county
Miguel I. Paredes<sup>1,2,*</sup>, Amanda C. Perofsky<sup>3,4</sup>, Lauren Frisbie<sup>5</sup>, Louise H. Moncla<sup>6</sup>, Pavitra Roychoudhury<sup>2,7</sup>, Hong Xie<sup>7</sup>, Shah A. Mohamed Bakhash<sup>7</sup>, Kevin Kong<sup>7</sup>, Isabel Arnould<sup>7</sup>, Tien V. Nguyen<sup>7</sup>, Seffir T. Wendm<sup>7</sup>, Pooneh Hajian<sup>7</sup>, Sean Ellis<sup>7</sup>, Patrick C. Mathias<sup>7</sup>, Alexander Greninger<sup>2,7</sup>, Lea M. Starita<sup>3,8</sup>, Chris D. Frazar<sup>8</sup>, Erica Ryke<sup>8</sup>, Weizhi Zhong<sup>3</sup>, Luis Gamboa<sup>3</sup>, Machiko Threlkeld<sup>8</sup>, Jover Lee<sup>2</sup>, Jeremy Stone<sup>3</sup>, Evan McDermot<sup>3</sup>, Melissa Truong<sup>8</sup>, Jay Shendure<sup>3,8,9</sup>, Hanna Oltean<sup>5</sup>, Cécile Viboud<sup>4</sup>, Helen Chu<sup>10</sup>, Nicola F. Müller<sup>2,†</sup>, Trevor Bedford<sup>1,2,3,8,9,†</sup>

<sup>1</sup>Department of Epidemiology, University of Washington, Seattle, WA, USA
<sup>2</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, Washington, USA
<sup>3</sup>Brotman Baty Institute for Precision Medicine, University of Washington, Seattle, WA USA
<sup>4</sup>Fogarty International Center, National Institutes of Health, Bethesda, MD, USA
<sup>5</sup>Washington State Department of Health, Shoreline, WA USA
<sup>6</sup>The University of Pennsylvania, Department of Pathobiology, Philadelphia, PA
<sup>7</sup>Department of Laboratory Medicine and Pathology, University of Washington, Seattle, WA, USA
<sup>8</sup>Department of Genome Sciences, University of Washington, Seattle, WA, USA
<sup>9</sup>Howard Hughes Medical Institute, Seattle, WA, USA
<sup>10</sup>Department of Medicine, Division of Allergy and Infectious Diseases, University of Washington, Seattle, WA
<sup>*</sup>To whom correspondence should be addressed.
<sup>†</sup>These authors jointly supervised this work.

## Abstract

SARS-CoV-2 transmission is largely driven by heterogeneous dynamics at a local scale, leaving local health departments to design interventions with limited information. We analyzed SARS-CoV-2 genomes sampled between February 2020 and March 2022 jointly with epidemiological and cell phone mobility data to investigate fine scale spatiotemporal SARS-CoV-2 transmission dynamics in King County, Washington, a diverse, metropolitan US county. We applied an approximate structured coalescent approach to model transmission within and between North King County and South King County alongside the rate of outside introductions into the county. Our phylodynamic analyses reveal that following stay-at-home orders, the epidemic trajectories of North and South King County began to diverge. We find that South King County consistently had more reported and estimated cases, COVID-19 hospitalizations, and longer persistence of local viral transmission when compared to North King County, where viral importations from outside drove a larger proportion of new cases. Using mobility and demographic data, we also find that South King County experienced a more modest and less sustained reduction in mobility following stay-at-home orders than North King County, while also bearing more socioeconomic inequities that might contribute to a disproportionate burden of SARS-CoV-2 transmission. Overall, our findings suggest a role for local-scale phylodynamics in understanding the heterogeneous transmission landscape.

## Organization

This repository contains the analytic code needed to reproduce the results from the above paper. To start, begin with the folder `nextstrain_build` to run the maximum likelihood analysis and create the temporally resolved phylogenies of SARS-CoV-2 in King County.

Then `nextstrain_build/cluster_assignment/` contains the script `KCparsimonyclusters.m` needed to assign local outbreak clusters needed to conduct the MASCOT-GLM analysis

After cluster assignment, navigate to `mascot_glm/data_for_tsv/` to combine all the cluster assignments from the four variant builds into one composite .tsv

`mascot_glm/mascotGLM.m` contains the code to produce the xmls for the three subsampling schemes that are inputed into MASCOT GLM

To run MASCOT GLM, while in the `mascot_glm/`folder in your terminal, type `java -jar Nab.jar xml/{xml_name}` for each separate subsampling scheme.

To analyze the posterior set of trees and log files produced by MASCOT-GLM, go to `analysis/post-process/` for all scripts

Finally, to produce the manuscript figures, run the notebooks found in `analysis/figures/`

All sequence accession IDs and associated North/South King County divisions used in this manuscript and provided by the Washington State Department of Health can be found in `final_acknowledgements_gisaid.csv`

Raw Beast Results and MASCOT-GLM XMLs can be found in `mascot_glm/results` and `mascot_glm/xmls`. MCC trees and any other results that require posterior processing can be found under `analysis/data-files/`

