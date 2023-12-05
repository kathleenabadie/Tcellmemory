# Tcellmemory_Abadie-et-al-2024
This repository includes code for analysis of genomics, imaging, and flow cytometry datasets for Abadie et al., 2024, " Reversible, tunable epigenetic silencing of TCF1 generates flexibility in the T cell memory decision." These scripts are minimally commented to be functional and reproducible, but they are not tutorials. Please see the Methods section in the manuscript for additional details. 

## scifate_analysis
This scripts are related to Figure 2 and Figure S2.
scifate_clustering.Rmd contains umap dimensional reduction and clustering of sci-fate-seq dataset (R).
scifate_DEG.Rmd contains differential expression analyses between umap clusters and between pseudotime trajectories (R).
Dynamo_Tcell_newRNA_Tcell_workflow.ipynb contains RNA velocity analyses using Dynamo (Python). 

## image_analysis_matlab
These scripts are related to Figure 3 and Figure S3, as well as Figure 4E-F. 
2022_ictrack_general_public.tar.gz contains scripts for image pre-processing, cell segmentation, and cell tracking. The ictrack software is adapted from Kueh et al., Nature Immunology, 2016 and Ng et al., eLife, 2018.
2022_longtracks_general_public.tar.gz contains scripts used to generate YFP intensity tracks over time that are used as imput to image_analysis_downstream_R analyses.
calc_YFP_fractions.m, plot_time_YFP_hist_1D.m, time_YFP_hist.m are used to convert time-lapse imaging output to YFP histograms and YFP+ fraction line plots in Figure 4E-F.

## image_analysis_downstream_R
These scripts are related to Figure 3 and Figure S3.
HMM_promoter_states_analysis_manuscript.Rmd implements a hidden markov model to call Tcf7 promoter states from the live imaging data and generate subsequent Tcf7 silencing analyses (R).

## hom_het_analysis


## bulk_genomics_analysis


## modeling_simulations_matlab


