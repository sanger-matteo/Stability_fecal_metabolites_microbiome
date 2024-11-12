# +++++ DESCRIPTION - Figure_Suppl_4.R ++++++++++++++++++++++++++++++++++++++
#
# Create Supplementary Figure 4, which are CoV of microbial taxa and metabolites
# grouped by individuals. There are two subfigures:
# - SF 4 A - CoV of microbial (genera)       --> hnd_3
# - SF 4 B - CoV of metabolites              --> hnd_4
#
# Essentially, this scripts is based on:
# - CoV_analysis.R
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2023
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Load Libraries -----------------------------------------------------------
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)

# Define function to perform division by zero and return zero (not NaN)
"%/0%" <- function(x,y) { res <- x/y;   res[is.na(res)] <-0;  return(res) }
# General plotting setting
figFontSize <- 12
theme_BoxPlot <- theme(
      plot.margin = margin(t=0.3, r=0.5, b=0.5, l=0.3, "cm"),
      panel.background = element_blank() ,
      panel.grid.minor = element_blank() ,
      panel.grid.major.x = element_line( size = 0.5, colour = "#cbcbcb", linetype="dashed") ,
      panel.grid.major.y = element_line( size = 0.5, colour = "#cbcbcb", linetype="solid") ,
      plot.title  = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
      axis.title  = element_text( color="#444444", size=figFontSize ) ,
      axis.text   = element_text( color="#444444", size=figFontSize   ) ,
      axis.ticks  = element_blank(),
      legend.position = "none",
      legend.background = element_blank() ,
      legend.title = element_text(size = figFontSize), 
      legend.text  = element_text(size = figFontSize),
      axis.line.x = element_line(size = 1, colour = "#444444"),
      axis.text.x = element_text(angle = 45, hjust=1) ,
      axis.text.y  = element_text( color="#444444", size=figFontSize-2 ) 
)



aggreg_at_rank <- "Family"


# ---- Load DATA tables --------------------------------------------------------------
cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures/")

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
# Prune the TAXO table removing descriptive columns and simplyfying names  
# (e.g., strip off the "__" separator and "DX" prefixes)
TAXO <- TAXO[4:10, -1]
for( ii in seq(1, dim(TAXO)[1]) ){
  names <- unlist(TAXO[ ii ,], use.names=FALSE)
  temp  <- as.data.frame( sapply( strsplit( names, "__" ), head, 2) )
  names <- as.character(as.vector( temp[2,] ))
  TAXO[ ii ,] <- as.list(names)
}

# Convert from count feature table to proportion (normalized)
ASV <- as.data.frame( dplyr::select( ASV, -c(Sample_ID, Sample_TotReads) ) )
ASV <- ASV/rowSums(ASV)
# Filtered taxonomic features with a Rel. abundance less than (0.03%) in greater than 10% of all samples
thres_abundance  <- 0.0003
thres_freqSample <- 0.010
mask <- unlist( colSums(ASV > thres_abundance) / dim(ASV)[1] >= thres_freqSample, use.names=FALSE )
ASV <- ASV[,mask]
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
# MB_normal <- cbind( metabols[, 1] , MB_normal)

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



# --- Plotting CoV by INDIVIDUALS - Metabolites  -----------------------------------------
# Group by individual and calculate the coefficient of variation
res <- data.frame()
for(ii in unique(sub_MetaData$Exp_Name)) {
  temp  <- MB_normal[ sub_MetaData$Exp_Name==ii, ]
  stdev <- as.data.frame( lapply( temp, function(x) sd(x)  ) )
  avg   <- as.data.frame( lapply( temp, function(x) mean(x)) )
  cov   <- stdev / avg
  res   <- rbind( res, c( "Exp_Name"=ii, cov) )
}
res <- replace(res, is.na(res), 0)
colnames(res) <- c( "Exp_Name", colnames(MB_normal))
cov_val <- res

dsPlot <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                    names_to  = "Feature",     # selected column for a new variable ("groupings")
                                    values_to = "CoV")         # values all into one new variables

# Import metabolites general classification and set order to match color palette
metabo_Ranking <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Ranking.txt"),
                            col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                            guess_max  = min( 10 , 300 ) )
metabo_Ranking <- as.data.frame(metabo_Ranking)
metabo_Ranking <- metabo_Ranking %>% dplyr::select( c("Feature_Name","F_Type","Colour") )
colnames(metabo_Ranking) <- c("Feature","F_Type","Colour")
metabo_Ranking <- metabo_Ranking[ order(metabo_Ranking$F_Type), ]
dsPlot <- merge( dsPlot, metabo_Ranking, by="Feature", all=TRUE)
# Factorize the "Feature" and "F_Type" 
dsPlot$Feature <- factor( dsPlot$Feature, levels=metabo_Ranking$Feature )
dsPlot$F_Type  <- factor( dsPlot$F_Type, levels= unique(sort(dsPlot$F_Type)) ) 
# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot %>% distinct( F_Type, .keep_all=TRUE)
reference_Factors <- reference_Factors[ order(reference_Factors$F_Type) , ]

# ----- 
# Plot
dsPlot_3 <- dsPlot[ !is.na(dsPlot$CoV) , ]
hnd_3 <- ggplot( dsPlot_3, aes(x=forcats::fct_rev(Exp_Name), y=CoV)) + 
          geom_boxplot( fill="#DCA32A", color="#666666", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
          geom_jitter( color="#666666", alpha=0.5, size=0.75, width=0.2) +
          coord_flip() +
          ylim( -0.1, 2.1) +
          labs(fill="Metabolites Type:") +
          xlab("") + ylab("CoV") +
          theme_BoxPlot +
          theme(legend.position = "none")
hnd_3


# --- Plotting CoV by INDIVIDUALS - Microbiome  -----------------------------------------
# Group by individual and calculate the coefficient of variation
cov_val <- cbind(sub_MetaData,  reduc_ASVs) %>% 
          group_by(Exp_Name) %>%
          reframe(  across( colnames(reduc_ASVs), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
          as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)

dsPlot <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                    names_to  = "Feature",     # selected column for a new variable ("groupings")
                                    values_to = "CoV")         # values all into one new variables

# Import color LUT for microbiome ranks and set order to match color palette
upper_rank <- "Class"
ref_ranks  <- read_tsv(paste0(path_2_data, "Data_Seq/RefList_",upper_rank,".txt"),
                       col_names = FALSE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                       guess_max  = min( 10 , 300 ) )

# The columns names from aggreg_at_rank and upper_rank, those that appear in dsPlot dataset
mask      <- reduc_TAXO[[aggreg_at_rank]] %in% dsPlot$Feature
reduc_TAXO  <- reduc_TAXO[ mask , c(upper_rank, aggreg_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")
# First, merge reduc_TAXO with ref_rank to add the "Color" column (matching upper_rank ID) to reduc_TAXO
reduc_TAXO <- merge( reduc_TAXO, ref_ranks, by=upper_rank)
# Second, select only unique elements of "aggreg_at_rank" (since it is in long-form)
# and then merge with the dsPlot, to add Colors and upper_rank to final plotting dataset
temp <- unique( reduc_TAXO[, c(aggreg_at_rank,"Colour",upper_rank)])
colnames(temp)[1] <- "Feature"
dsPlot <- merge( temp, dsPlot, by="Feature", all=TRUE)
# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsPlot$Feature <- sapply( dsPlot$Feature , function(x)  gsub("\\[", "", x) )
dsPlot$Feature <- sapply( dsPlot$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" and "F_Type" 
dsPlot$Feature <- factor( dsPlot$Feature, levels=ordered_Features)
dsPlot[upper_rank] <- factor( dsPlot[[upper_rank]], levels=ref_ranks[[upper_rank]] )

# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot %>% distinct( .[[upper_rank]], .keep_all=TRUE)

# Remove rows that had empty calculation 0/0 ( there are only a handful)
dsPlot <- dsPlot[ !is.na(dsPlot$Exp_Name) ,]

# ----- 
dsPlot_4 <- dsPlot
dsPlot_4[dsPlot_4$CoV == 0.0 , "CoV"] <- NaN
hnd_4 <- ggplot( dsPlot_4, aes(x=forcats::fct_rev(Exp_Name), y=CoV)) + 
          geom_boxplot( fill="#2593BD", color="#666666", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
          geom_jitter( color="#666666", alpha=0.5, size=0.8, width=0.2) +
          coord_flip() +
          ylim( -0.1, 2.1) +
          xlab("") + ylab("CoV") +
          theme_BoxPlot +
          theme(legend.position = "none")

hnd_4



# --- FIGURE - combine sub-plots -------------------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")

ggdraw() +
  draw_plot(hnd_4, x = 0.15, y = 0, width = .3, height = .92) +
  draw_plot(hnd_3, x = 0.55, y = 0, width = .3, height = .92) +
  draw_plot_label(label = c("A", "B"), size = 17, x = c(0.12, 0.52), y = c(1,1))

ggsave( paste0(path_saveFig, "Suppl_Figure_4.tiff"), width=10, height=6, dpi=300)










 