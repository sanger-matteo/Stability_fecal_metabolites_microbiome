# +++++ DESCRIPTION - Figure_Suppl_3_ICC.R ++++++++++++++++++++++++++++++++++++++
#
# Calculate the IntraClass Correlation (ICC).
# This Script calculate and plot the the IntraClass Correlation (ICC) for both
# the microbiome taxa and metabolites
#
# DEF: ICC quantifies the degree to which individuals with a fixed degree of 
# relatedness (assumption) resemble each other in terms of a quantitative trait.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2023
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(nlme)
library(nortest)
library(ggplot2)

# Define the ICC function: very simple, no package needed
icc.lme <- function(lmeobj)
{
  varcorr <- VarCorr(lmeobj)
  return( as.numeric(varcorr[[1]]) / (as.numeric(varcorr[[1]])+as.numeric(varcorr[[2]]) ) )
}  

# General plotting setting
figFontSize <- 12
theme_BoxPlot <- theme(
      plot.margin = margin(t=0.3, r=0.5, b=0.5, l=0.3, "cm"),
      panel.background = element_blank() ,
      panel.grid.minor = element_blank() ,
      panel.grid.major.x = element_line( size = 0.5, colour = "#cbcbcb", linetype="dashed") ,
      panel.grid.major.y = element_line( size = 0.5, colour = "#cbcbcb", linetype="solid") ,
      plot.title  = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
      axis.title  = element_text( color="#444444", size=figFontSize+1 ) ,
      axis.text   = element_text( color="#444444", size=figFontSize   ) ,
      axis.ticks  = element_blank(),
      legend.position = "none",
      legend.background = element_blank() ,
      legend.title = element_text(size = figFontSize), 
      legend.text  = element_text(size = figFontSize-2),
      axis.line.x = element_line(size = 1, colour = "#4e4e4e"),
      axis.text.x = element_text(angle = 45, hjust=1) ,
      axis.text.y  = element_text( color="#444444", size=figFontSize ) 
)

cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures/")



# ---- Load DATA tables --------------------------------------------------------------

# Load ASV table(s) and Metadata
ASV  <- read_tsv( paste0( path_2_data,"Data_Seq/ASV_longit.tsv"),
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) )
TAXO <- read_tsv( paste0( path_2_data,"Data_Seq/TAXO_longit.tsv"),
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) ) 

MetaData <- read_tsv( paste0(path_2_data, "Data_Seq/META_longit.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

metabols <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Concentration.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 



# --- Pruning MICROBIOME data ------------------------------------------------------------------- 
# Prune the TAXO table removing descriptive columns and simplifying names  
# (e.g., strip off the "__" separator and "DX" prefixes)
TAXO <- TAXO[4:10, -1]
for( ii in seq(1, dim(TAXO)[1]) ){
  names <- unlist(TAXO[ ii ,], use.names=FALSE)
  temp  <- as.data.frame( sapply( strsplit( names, "__" ), head, 2) )
  names <- as.character(as.vector( temp[2,] ))
  TAXO[ ii ,] <- as.list(names)
}

# Convert from count feature table to relative abundances (normalized)
ASV <- as.data.frame( dplyr::select( ASV, -c(Sample_ID, Sample_TotReads) ) )
ASV <- ASV/rowSums(ASV)
# Filter taxonomic features with a Rel. abundance less than 0.03%, and in less than 10% of the cohort
thres_abundance  <- 0.0003
thres_freqSample <- 0.010
mask <- unlist( colSums(ASV > thres_abundance) / dim(ASV)[1] >= thres_freqSample, use.names=FALSE )
ASV  <- ASV[,mask]
TAXO <- TAXO[,mask]

# AGGREGATE and reduce ASVs at a specific taxa level ("aggreg_at_rank"), 
taxa_rank <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
aggreg_at_rank <- "Genus"
reduc_ASVs <- ASV[ ]
reduc_TAXO <- as.data.frame(t(TAXO))
colnames(reduc_TAXO) <- taxa_rank
reduc_TAXO$Sample <- rownames(reduc_TAXO)

# Determine unique_taxa and remove those that are unclassified
drop_feat_named <- read_tsv( paste0(path_2_data, "Data_Seq/Drop_Features_Named.txt") )
drop_feat_named <- drop_feat_named[[1]]
unique_taxa <- unique( unlist( reduc_TAXO[aggreg_at_rank] , use.names=FALSE) )
unique_taxa <- unique_taxa[ ! unlist( lapply( unique_taxa , function(x) any( stringr::str_detect( x ,drop_feat_named)) )) ]

temp <- list()
for (rr in seq(1, length(unique_taxa)) ) {
  selecols <- intersect( colnames(reduc_ASVs) , reduc_TAXO[reduc_TAXO[aggreg_at_rank] == unique_taxa[rr], "Sample"] )
  temp[[rr]] <- rowSums( data.frame(reduc_ASVs[ , selecols ]) , na.rm=F)
}
reduc_ASVs <- as.data.frame(temp)
colnames(reduc_ASVs) <- unique_taxa


# --- Pruning METABOLITES data ------------------------------------------------------------------- 
# Reduce the sample names in MetaData to the essentials and match metabolites tables
MetaData$External_ID <- lapply( MetaData$External_ID, function(x) sub( "Exp1_DNA_", "", x ) )

# Subset the samples used for analysis
mask <- metabols$External_ID %in%  MetaData$External_ID 
metabols <- metabols[ mask, ]

# Reorganize metabols dataset rows to match the order of MetaData df
metabols$External_ID <- factor( metabols$External_ID, levels=MetaData$External_ID )
metabols <- metabols[ order(metabols$External_ID), ]

# Normalized metabols areas of entire dataset
MB_normal <- metabols[, -c(1)] 
MB_normal <- sapply( as.data.frame(t( MB_normal )), function(x) x/sum(x))
MB_normal <- as.data.frame( t( MB_normal ) )
colnames(MB_normal) <- colnames(metabols[, -c(1)] )
rownames(MB_normal) <- seq(1,dim(MB_normal)[1])

# Subset MetaData to the useful response variables and factorize them
sele_cols <- c("Exp_TimePoint", "Exp_Name", "BMI_cat", "Age_Rank", "Sex" )
sub_MetaData <- MetaData %>% dplyr::select( all_of(sele_cols) )
for (ii in seq(1,length(sele_cols))) {
  ii_Col <- sele_cols[ii]
  sub_MetaData[[ii_Col]] <- MetaData[[ii_Col]] %>% factor()
  levels(sub_MetaData[[ii_Col]]) <- seq( 1, length(unique(MetaData[[ii_Col]])))
  sub_MetaData[[ii_Col]] <- as.numeric(sub_MetaData[[ii_Col]])
}
# Transform the Exp_Name into string
sub_MetaData$Exp_Name <- sprintf( sub_MetaData$Exp_Name, fmt = '%#.2d') 
sub_MetaData$Exp_Name <- sapply( sub_MetaData$Exp_Name, function(x) paste0("A", x) )



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Perfor ICC on MICROBIOME -------------------------------------------------------
Data <- reduc_ASVs *100
Individuals <- MetaData$Exp_Name
n_features <- dim(Data)[2]
ICC_res <- vector('numeric',n_features)

Data <- compositions::clr( Data )

for (ii in 1:n_features) {
  y_Feature <- Data[ , ii]
  resultat  <- lme( y_Feature ~1, random = ~1 | as.factor(Individuals) , na.action=na.omit ) 
  resultatY <- summary(resultat)
  temp <- resultatY[["tTable"]]
  ICC_res[ii] <- icc.lme(resultat)
}
df = data.frame( Feature_Name = colnames(Data), ICC = ICC_res)
avg_icc_ubiome = mean(df$ICC) 

dsPlo_A <- df
colnames(dsPlo_A)[1] <- "Feature"

# Import color LUT for microbiome ranks and set order to match color palette
upper_rank <- "Class"
ref_ranks  <- read_tsv(paste0(path_2_data, "Data_Seq/RefList_",upper_rank,".txt"),
                       col_names = FALSE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                       guess_max  = min( 10 , 300 ) )

# The columns names from aggreg_at_rank and upper_rank, those that appear in dsPlo_A dataset
mask      <- reduc_TAXO[[aggreg_at_rank]] %in% dsPlo_A$Feature
sub_TAXO  <- reduc_TAXO[ mask , c(upper_rank, aggreg_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")
# First, merge reduc_TAXO with ref_rank to add the "Color" column (matching upper_rank ID) to reduc_TAXO
sub_TAXO <- merge( sub_TAXO, ref_ranks, by=upper_rank)
# Second, select only unique elements of "aggreg_at_rank" (since it is in long-form)
# and then merge with the dsPlo_A, to add Colors and upper_rank to final plotting dataset
temp <- unique( sub_TAXO[, c(aggreg_at_rank,"Colour",upper_rank)])
colnames(temp)[1] <- "Feature"
dsPlo_A <- merge( temp, dsPlo_A, by="Feature", all=TRUE)
# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsPlo_A$Feature  <- sapply( dsPlo_A$Feature , function(x)  gsub("\\[", "", x) )
dsPlo_A$Feature  <- sapply( dsPlo_A$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" and "F_Type" 
dsPlo_A$Feature <- factor( dsPlo_A$Feature, levels=rev(ordered_Features))  
dsPlo_A[upper_rank] <- factor( dsPlo_A[[upper_rank]], levels=ref_ranks[[upper_rank]] )
# Now subset to take only the "top" taxa
dsPlo_A <- dsPlo_A[dsPlo_A$Feature %in% ordered_Features[1:51], ]

# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlo_A %>% distinct( .[[upper_rank]], .keep_all=TRUE)
reference_Factors[upper_rank] <- factor( reference_Factors[[upper_rank]], levels=ref_ranks[[upper_rank]] ) 
reference_Factors <- reference_Factors[ order(reference_Factors[[upper_rank]]), ]


# Set Plotting theme and colors
my_palette_A <- reference_Factors$Colour

hnd_A <- ggplot( dsPlo_A, aes( x=Feature, y=ICC) )  +
        geom_point( aes(fill=.data[[upper_rank]]), color="#666666", size=3, alpha = 0.8, shape=23) +  
        geom_hline(yintercept=avg_icc_ubiome, linetype="dashed", alpha = 0.6) +
        coord_flip() + 
        scale_fill_manual( values=my_palette_A ) +
        ylim(0,1) + 
        labs(x = "", colour="Metabolite Group:") +
        theme_BoxPlot
hnd_A



# --- Perform ICC on Metabolites -------------------------------------------------------
Data <- MB_normal *100
Individuals <- MetaData$Exp_Name
n_features <- dim(Data)[2]
ICC_res_B <- vector('numeric',n_features)
hnd_B <- list()

for (ii in 1:n_features) {
  colnames(Data)[[ii]]
  y_Feature <- Data[ , ii]
  resultat  <- lme( y_Feature ~1, random = ~1 | as.factor(Individuals) , na.action=na.omit ) 
  resultatY <- summary(resultat)
  temp <- resultatY[["tTable"]]
  ICC_res_B[ii] <- icc.lme(resultat)
  qqnorm( residuals(resultat) ) 
  qqline( residuals(resultat) ) 
}
df = data.frame( Feature_Name = colnames(Data), ICC = ICC_res_B)
avg_icc_metabolites = mean(df$ICC)

metab_classes <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Ranking.txt"),
                           col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                           guess_max  = min( 10 , 300 ) )
metab_classes <- as.data.frame(metab_classes)
metab_classes <- metab_classes %>% dplyr::select( c("Feature_Name","F_Type","Colour") )
dsPlot_B <- merge( df, metab_classes, by="Feature_Name", all=TRUE) 
dsPlot_B$Feature_Name <- factor( dsPlot_B$Feature_Name, levels=rev( unique(metab_classes$Feature_Name)) )
# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot_B %>% distinct( F_Type, .keep_all=TRUE)
reference_Factors <- reference_Factors[ order(reference_Factors$F_Type) , ]

my_palette_B <- reference_Factors$Colour

hnd_B <- ggplot( dsPlot_B, aes( x=Feature_Name, y=ICC) )  +
          geom_point( aes(fill=F_Type), color="#666666", size=3, alpha = 0.8, shape=23) +   
          geom_hline( yintercept=avg_icc_metabolites, linetype="dashed", alpha = 0.6) +
          coord_flip() + 
          scale_fill_manual( values=my_palette_B ) +
          ylim(0,1) + 
          labs(x = "", colour="Metabolite Group:") +
          theme_BoxPlot
hnd_B



# --- COMBINE sub-Figures ------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")
ggdraw()

ggdraw() +
  draw_plot(hnd_A, x = 0, y = .28, width = .55, height = .71) +
  draw_plot(hnd_B, x = 0.55, y = 0, width = .45, height = .99) +
  draw_plot_label(label = c("A", "B"), size = 17, x = c(0.07, 0.58), y = c(1, 1))

ggsave( paste0(path_saveFig, "Suppl_Figure_3_ICC.tiff"), width=10, height=12.5, dpi=300)












