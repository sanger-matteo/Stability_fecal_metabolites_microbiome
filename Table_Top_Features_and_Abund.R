##### Description --- Top_Features_and_Abund ###################################
#
# The script create summary tables for the microbiome and metabolite features,
# specifically the relative abundance and frequency in the cohort, grouping by 
# individuals, sample, or cohort. These tables generate the Table 1 and Table 2
# of the paper 
# 
# Output:
# The tables are saved in folder Res_Tables/ (inside path_2_data)
# - Table_Avg_Metabolites.txt
# - Concentrations_Metabolites.txt
# - uBiome_RelAbu_Class_bySamples.txt
# - uBiome_RelAbu_Class_byCohort.txt
# - uBiome_RelAbu_Class_byIndividual.txt
# The last three can be generated for any taxonomic rank chose (here "Class")
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


calculate_stats <- function(column) {
  # Calc. and return the mean, standard deviation, median and interquartile
  # range for one column of a data.frame 
  mean_value <- mean(column)
  sd_value <- sd(column)
  median_value <- median(column)
  iqr_value <- IQR(column)
  return(c(mean_value, sd_value, median_value, iqr_value))
}

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


# --- Pruning the two datasets
aggreg_at_rank <- "Genus"
source(paste0( path_2_data, "pruning_2Datasets.R"))



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # --- Tables of MICROBIOME dataset  --------------------------------------------------
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
# # --- Tables of METABOLITES dataset ----------------------------------------------------
# Save the table with the measured metabolites levels
metab_Avg_SD <- cbind( unlist(MetaData$External_ID), MB_normal )
colnames( metabols )[1] <- "Sample_ID"
write_tsv( metabols, paste0( path_2_data, "Res_Tables/Concentrations_Metabolites.txt") )

# Apply the function "calculate_stats" to each column of the data frame
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
 














