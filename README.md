This repository contains scripts associated with the paper:  

## Stability in fecal metabolites amid a diverse gut microbiome composition: a one-month longitudinal study of variability in healthy individuals. 
Published in *Gut Microbes*, 12 November 2024.
[Link to the paper](URL_of_your_paper)
[DOI: 10.1080/19490976.2024.2427878](https://doi.org/10.1080/19490976.2024.2427878)

#### Authors
**Matteo Sangermani**<sup>1,2</sup>, Indri Desiati<sup>3</sup>, Solveig M. Jørgensen<sup>3</sup>, Jia V. Li<sup>4</sup>, Trygve Andreassen<sup>3,5</sup>, Tone F. Bathen<sup>3</sup>, Guro F. Giskeødegård<sup>1,2</sup>

#### Affiliations
<h6>
1. Dep. of Public Health and Nursing NTNU, Trondheim
2. Dep of Surgery, St. Olavs University Hospital, Trondheim
3. Dep. of Circulation and Medical Imaging, NTNU, Trondheim.
4. Department of Metabolism, Digestion and Reproduction, Imperial College London, London
5. Central Staff, St. Olavs Hospital HF, Trondheim.
</h6>


## Abstract of the paper
<h6l>*An extensive network of microbial-host interactions exists in the gut, making the gut microbiome a complex ecosystem to untangle. The microbial composition and the fecal metabolites are important readouts to investigate intricate microbiota-diet-host interplay. However, this ecosystem is dynamic, and it is of interest to understand the degree and timescale of changes occurring in the gut microbiota, during disease as well as in healthy individuals. Cross-sectional study design is often used to investigate the microbiome, but this design provides a static snapshot and cannot provide evidence on the dynamic nature of the gut microbiome. Longitudinal studies are better suited to extrapolate causation in a study or assess changes over time.
This study investigates longitudinal change in the gut microbiome and fecal metabolites in 14 healthy individuals with weekly sampling over a period of one-month (four time points), to elucidate the temporal changes occurring in the gut microbiome composition and fecal metabolites. Utilizing 16S rRNA amplicon sequencing for microbiome analysis and NMR spectroscopy for fecal metabolite characterization, we assessed the stability of these two types of measurable parameters in fecal samples during the period of one month. Our results show that the gut microbiome display large variations between healthy individuals, but relatively lower within-individual variations, which makes it possible to uniquely identify individuals. The fecal metabolites showed higher stability over time compared to the microbiome and exhibited consistently smaller variations both within and between individuals. This relative higher stability of the fecal metabolites suggests a balanced, consistent output even amid individual’s differences in microbial composition and they can provide a viable complementary readout to better understand the microbiome activity.*</h6>


## Description of the Scripts

This repository includes several scripts used in the study. Below is a brief overview of each:

- **Figure_1_Composition.R**: &nbsp;&nbsp;&nbsp;*Generates composition plot and tSNE clustering.*
- **Figure_2_CoV.R**: &nbsp;&nbsp;&nbsp;*Analyzes the Coefficient of Variation (CoV).*
- **Figure_3_Concentration_Plots.R**: &nbsp;&nbsp;&nbsp;*Tracks time evolution of specific microbiome features and fecal metabolites.*
- **Figure_3_PartialSpearman.R**: &nbsp;&nbsp;&nbsp;*Correlates microbiome and fecal metabolites.*
- **Figure_Suppl_1.R**: &nbsp;&nbsp;&nbsp;*Analyzes alpha and beta diversity.*
- **Figure_Suppl_2.R**: &nbsp;&nbsp;&nbsp;*Calculates CoV of microbial features at taxonomic ranks: class, family, genus.*
- **Figure_Suppl_3_ICC.R**: &nbsp;&nbsp;&nbsp;*Interclass correlation analysis for microbial and fecal metabolites.*
- **Figure_Suppl_4.R**: &nbsp;&nbsp;&nbsp;*CoV of microbial taxa and metabolites grouped by individuals.*
- **Figure_Suppl_5.R**: &nbsp;&nbsp;&nbsp;*Scatter plot of IQR vs. Median from CoV values.*
- **Figure_Suppl_6.R**: &nbsp;&nbsp;&nbsp;*Concentration plot analysis on a different feature set.*

Each script is designed to be run independently to generate the figures presented in the published paper.

## Description of the Data
*Provide an overview of the datasets used, generated, or required by the scripts. Mention if there are specific data preprocessing steps or data formats expected.*

## Main Figures
The following figures showcase the types of analyses or results generated by the scripts in this repository:

| ![Figure 1](Final_Figures/Figure_1.png) | ![Figure 2](Final_Figures/Figure_2.png) | ![Figure 3](Final_Figures/Figure_3.png) |
|---------------------------------------|---------------------------------------|---------------------------------------|
| Description of Figure 1               | Description of Figure 2               | Description of Figure 3               |

## Usage

To run these scripts, the following R packages are required:

#### Visualization
- `ggplot2`
- `cowplot`
- `tiff`
- `ggrepel`
- `ggtext`
- `showtext`

#### Data Manipulation
- `readr`
- `tidyr`
- `tibble`
- `dplyr`

#### Statistical Analysis
- `Rtsne`
- `nlme`
- `nortest`
- `PResiduals` (for Partial Spearman correlation)

#### Setup Instructions
In each script, change the variable `cDir_Rscript` to the **absolute path** of the repository on the machine running the code. All other paths are relative to this and follow Unix/MacOS format.


## Data Tables and Format:
Table of the data used for analysis are included in the repository. Folder *Data_NMR* contains the quantification of fecal metabolites measured from NMR spectras, whereas  *Data_Seq* contains the count data and taxonomic information of microbiome features obtained with 16S rRNA sequencing. Below is an overview of the files:
- *Data_NMR*
	- *Metabolites_Concentration.txt*: quantification of identified peaks.
	- *Metabolites_MatchFactor.txt*: peak fitting evaluation from Chenomx NMR Suite.
	- *Metabolites_Ranking.txt* and *Metabolites_Ranking.xls*: a classification and ranking of all the metabolites identified in the spectra (using Chenomx NMR Suite software) with additional information. NOTE: column F_Type is our custom. classification, based on HMDB classification.
	- *Ctrl_Repeated_Measurments.txt*: results of NMR repeated measurements from the same samples, to as control.
	- *LUT_Metabolites*: LUT color table; each metabolite in the cohort is assigned a range of tonalities from the same base color.
	
- *Data_Seq*
	- *ASV_longit.tsv*: read count table of the cohort, where each column is a microbial ASV.
	- *TAXO_longit.tsv*: taxonomic information for each ASV present in the cohort.
	- *META_longit*: The metadata information of the cohort, such as "Sample_ID", "Exp_Group", "Exp_Name", "Exp_TimePoint", "Sex", etc.
	- *RefList_Class.txt* and *RefList_Phyla.txt*: a simple list of all the classes and phyla present in the cohort listed by abundance and with a reference color to use to highlight plots by tese specific ranks.
	- *Drop_Features_Named.txt*: species taxonomy in 16S rRNS sequencing is often incomplete. THese is a list of names (such as "unkown", "unclassified genome", etc) that should not be identified as "species".
	- *LUT_Microbes*: LUT color table; each class in the cohort is assigned a range of tonalities from the same base color.
















