# Lake sturgeon *p*CO<sub>2</sub> and temperature
A repository for scripts related to a lake sturgeon *p*CO<sub>2</sub> and temperature exposure-based experiment. The following sections have explanations of repo directories and files.


## R_scripts:
These R scripts are divided into those used for looking at RNAseq data and those used for modeling different physiological and behavioural variables.
- `RNAseq`
  - [`edgeR.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/RNAseq/edgeR.R) contains code for running differential gene expression analyses with [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), making a PCA and heatmap, and writing out results tables.
  - [`enrichR_nonunique.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/RNAseq/enrichR_nonunique.R) contains code for using edgeR results with [enrichR](https://maayanlab.cloud/Enrichr/) to identify gene ontology terms specific, but not unique to each experimental treatment.
  - [`enrichR_overall.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/RNAseq/enrichR_overall.R) contains code for using edgeR results with enrichR to identify gene ontology terms unique to each experimental treatment.

- `stats`
  -  [`Ammonia.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/Ammonia.R) contains brms-based Bayesian analyses of ammonia excretion among the different experimental groups.
  -  [`NKA_corr.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/NKA_corr.R) contains Bayesian models of Na<sup>+</sup>/K<sup>+</sup> ATPase activity and transcript abundance for different subunits.
  -  [`behaviour.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/behaviour.R) contains models of responses to cues at 0 and 24h after a transient increase in *p*CO<sub>2</sub> of 10,000 μatm.
  -  [`boldness.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/boldness.R) contains models of activity and time near a [novel object](https://img.bricklink.com/ItemImage/MN/0/sw0036.png) during the middle 3 minutes of a 5 minute experiment.
  -  [`brms_ATPase.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/brms_ATPase.R) contains models of Na<sup>+</sup>/K<sup>+</sup> ATPase activity among different experimental groups.
  -  [`hematocrit.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/hematocrit.R) contains models of hematocrit levels among different groups.
  -  [`metabolic_rate.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/metabolic_rate.R) contains models of routine and maximum metabolic rates among the different groups.
  -  [`network_plot.R`](https://github.com/BioMatt/lkst_pCO2/blob/main/R_scripts/stats/network_plot.R) contains visualizations of gene betweenness scores in networks, found using [OmicsNet 2.0](https://www.omicsnet.ca/). 
 

## shell_scripts:
These shell scripts were run on the [Cedar](https://docs.alliancecan.ca/wiki/Cedar) cluster of the [Digital Research Alliance of Canada](https://alliancecan.ca/en).
  - [`fastp.sh`](https://github.com/BioMatt/lkst_pCO2/blob/main/shell_scripts/fastp.sh) uses the program [fastp](https://github.com/OpenGene/fastp) to perform quality control on the raw RNAseq data.
  - [`raw_reads.txt`](https://github.com/BioMatt/lkst_pCO2/blob/main/shell_scripts/raw_reads.txt) is not a script, but a text file with explicit file paths to raw data useful for running array jobs. This file can be easily created with a line like `printf '%s\n' "$PWD"/*_R1.fastq.gz > forward_reads.txt` and another for reverse reads, then putting them together.
  - [`salmon_index.sh`](https://github.com/BioMatt/lkst_pCO2/blob/main/shell_scripts/salmon_index.sh) uses the program [Salmon](https://combine-lab.github.io/salmon/) to index the lake sturgeon gill transcriptome for subsequent transcript quantification.
  - [`salmon_quant.sh`](https://github.com/BioMatt/lkst_pCO2/blob/main/shell_scripts/salmon_quant.sh) uses Salmon to quantify transcript abundance. The output files of this script would then be used with the `edgeR.R` script under the `R_scripts` folder.
