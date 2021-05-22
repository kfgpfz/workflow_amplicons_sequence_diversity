# workflow_amplicons_sequence_diversity

THe folder contains the scripts and workflow for the analysis of amplicon sequencing data published in the paper XXX.

Workflow:


First, for every amplicon pool the bash-script pipeline_SE.sh needs to be executed. Please adapt the the different path variables to your environment. A specific conda environment for this project is highly recommended. The script and the positions folder need to be in the same subdirectory.

bash pipeline_SE.sh

Second, for every locus the R scripts algG_analysis_preparation_200.R and algG_analysis_coverage.R are executed based on the data of the first step. The name of the locus, the orientation and filepath need to be adapted. The SNP statistics is performed using SNPeff (http://pcingola.github.io/SnpEff/). The command is, where ** is the name of the investigated locus.

snpEff -ud 0 PA14 **_SNPeff_input.vcf > **_SNPeff_results.vcf 


Finally, for every locus the R script Zusammenfassung_algG_200_20_envs_2021.Rmd needs the be executed. For other loci the name of the locus and file path need to be adapted.

The script Zusammenfassung_inter_gene_comparison_200_20_envs_2021.Rmd compares all the loci for specific SNP characteristics.
