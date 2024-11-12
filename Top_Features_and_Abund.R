##### Description --- Top_Features_and_Abund ###################################
#
# The script create simple summary tables for the relative abundance and 
# frequency of features, for the entire cohort, of grouping by individuals.
# This is done for both microbiome and metabolites data.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2023
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Load Libraries -----------------------------------------------------------
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)


# General plotting setting
figFontSize <- 13
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



# ---- Load DATA tables --------------------------------------------------------
cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )

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
aggreg_at_rank <- "Class"
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
sele_cols <- c("Sex", "Exp_TimePoint", "BMI_cat", "Age_Rank", "Height_cat", "Exp_Name")
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


# Gather the relative abundance of all samples (in %)
RelAbu_by_Samples <-cbind( "Sample_ID"=unlist(MetaData$External_ID,use.names=F), reduc_ASVs*100 )

write_tsv( RelAbu_by_Samples, paste0( path_2_data, "Res_Tables/uBiome_RelAbu_",aggreg_at_rank,"_bySamples.txt") )

# Gather the relative abundance (in %) grouping results by individual, thus, 
# calculate mean abundance for all the samples belonging to one individual
RelAbu_by_Indiv  <- data.frame()
for(ii in unique(sub_MetaData$Exp_Name)) {
  temp  <- reduc_ASVs[ sub_MetaData$Exp_Name==ii, ]
  avg   <- as.data.frame( lapply( temp, function(x) mean(x)) ) *100
  RelAbu_by_Indiv <- rbind( RelAbu_by_Indiv, c( "Exp_Name"=ii, avg) )
}
RelAbu_by_Indiv <- replace(RelAbu_by_Indiv, is.na(RelAbu_by_Indiv), 0)
colnames(RelAbu_by_Indiv) <- c( "Exp_Name", colnames(reduc_ASVs))

write_tsv( RelAbu_by_Indiv, paste0( path_2_data, "Res_Tables/uBiome_RelAbu_",aggreg_at_rank,"_byIndividual.txt") )



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# Summarize the data for entire cohort, by calculating frequency of a feature 
# in the cohort, average abundance and StDev (all samples in cohort)
# Calculate the frequency of Taxa
sub_ASVs  <- reduc_ASVs[ MetaData$External_ID != "A12_28" , ]   # Remove the sample with problematic sequencing 
RelAbu_by_Cohort <- as.data.frame( (colSums(sub_ASVs>0) / nrow(sub_ASVs)) *100)
colnames(RelAbu_by_Cohort) <- "Freq_Cohort"
# Calculate the average abundance and StDev
RelAbu_by_Cohort$Avg_abundance   <- unlist( lapply( RelAbu_by_Indiv[,-1], function(x) mean(x) ), use.names = F )
RelAbu_by_Cohort$StDev_abundance <- unlist( lapply( RelAbu_by_Indiv[,-1], function(x) sd(x) ), use.names = F )
RelAbu_by_Cohort$Median          <- unlist( lapply( RelAbu_by_Indiv[,-1], function(x) median(x) ), use.names = F )
RelAbu_by_Cohort$InterQuartRange <- unlist( lapply( RelAbu_by_Indiv[,-1], function(x) IQR(x) ), use.names = F )
RelAbu_by_Cohort[aggreg_at_rank] <- rownames(RelAbu_by_Cohort)
# Format and sort the final table
mask <- reduc_TAXO[[aggreg_at_rank]] %in% RelAbu_by_Cohort[[aggreg_at_rank]]
sub_TAXO  <- reduc_TAXO[ mask , ] 
temp <- unique( sub_TAXO[ , taxa_rank[ 2:which( taxa_rank==aggreg_at_rank ) ] ])
if( aggreg_at_rank != "Phyla") {
  RelAbu_by_Cohort <- merge( RelAbu_by_Cohort, temp,  by=aggreg_at_rank, all=TRUE)  }
RelAbu_by_Cohort <- RelAbu_by_Cohort[ order(RelAbu_by_Cohort$Avg_abundance, decreasing=T) , ]
# Rename the columns to pubblish ready format
if( aggreg_at_rank == "Phyla") {
  colnames(RelAbu_by_Cohort)[1:5] <- c( "Frequency",	"Average",	"St. Dev.", "Median", "Interq. Range" ) 
} else {
  colnames(RelAbu_by_Cohort)[2:6] <- c( "Frequency",	"Average",	"St. Dev.", "Median", "Interq. Range" ) }

write_tsv( RelAbu_by_Cohort, paste0( path_2_data, "Res_Tables/uBiome_RelAbu_",aggreg_at_rank,"_byCohort.txt") )





# # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # --- Plotting CoV of METABOLITES  -----------------------------------------------------
# Save the table with the measured metabolites levels
metab_Avg_SD <- cbind( unlist(MetaData$External_ID), MB_normal )
colnames( metabols )[1] <- "Sample_ID"
write_tsv( metabols, paste0( path_2_data, "Res_Tables/Concentrations_Metabolites.txt") )

# Calculate the average abundance and StDev of metabolites for the entire cohort
calculate_stats <- function(column) {
  # Calc. mean and standard deviation for each column of a data.frame
  mean_value <- mean(column)
  sd_value <- sd(column)
  median_value <- median(column)
  iqr_value <- IQR(column)
  return(c(mean_value, sd_value, median_value, iqr_value))
}

# Apply the function to each column of the data frame
metab_Avg_SD <- data.frame( t(apply( MB_normal, 2, calculate_stats)) )
colnames(metab_Avg_SD) <- c("Average", "St_Dev", "Median", "Interq. Range" )
metab_Avg_SD$Feature <- rownames(metab_Avg_SD)

# Import metabolites ranking to add information on F_type and sort it by average abundance
metabo_Ranking <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Ranking.txt"),
                            col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                            guess_max  = min( 10 , 300 ) )
metabo_Ranking <- as.data.frame(metabo_Ranking)
metabo_Ranking <- metabo_Ranking %>% dplyr::select( c("Feature_Name","F_Type") )
colnames(metabo_Ranking) <- c("Feature","F_Type")
metabo_Ranking <- metabo_Ranking[ order(metabo_Ranking$F_Type), ]
metab_Avg_SD <- merge( metab_Avg_SD, metabo_Ranking, by="Feature", all=TRUE)
metab_Avg_SD <- metab_Avg_SD[order(metab_Avg_SD[["Average"]], decreasing = TRUE), ]

colnames(metab_Avg_SD) <- c("Metabolite","Average","St_Dev", "Median", "Interq. Range", "Class")

write_tsv( metab_Avg_SD, paste0( path_2_data, "Res_Tables/Table_Avg_Metabolites.txt") )
 














