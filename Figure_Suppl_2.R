# +++++ DESCRIPTION - Figure_Suppl_2.R ++++++++++++++++++++++++++++++++++++++
#
# Create Supplementary Figure 2, which are CoV of microbial taxa grouped at 
# different taxonomic levels. There are three subfigures:
# - SF 2 A - CoV of microbial Classes                   --> hnd_1
# - SF 2 B - CoV of microbial Families                  --> hnd_2
# - SF 2 C - CoV of double-sorted microbial Genera      --> hnd_3
#
# This figures summarise the results from the analysis of the Coefficient of
# Variation (CoV) of microbiome aggregating features at different taxonomic
# ranks.
#
# All those plots are sorted by the respective taxa frequency in the cohort
# and color coded by Class; thus, no legend is needed and simply the color used 
# in sub-figure A act as a legend for B and C.
# In SF 2 C we sort by class, and each class is then itself sorted by frequency,
# so that genera are sorted by both hierarchy + frequency
#
# Essentially this scripts combines:
# - CoV_analysis.R
# - CoV_analysis_Control.R
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
      plot.title  = element_text( color="#222222", size=figFontSize, face = "bold" ) ,
      axis.title  = element_text( color="#444444", size=figFontSize +4 ) ,
      axis.text   = element_text( color="#444444", size=figFontSize   ) ,
      axis.ticks  = element_blank(),
      #legend.position = "none",
      legend.background = element_blank() ,
      legend.title = element_text(size = figFontSize), 
      legend.text  = element_text(size = figFontSize),
      axis.line.x = element_line(size = 1, colour = "#444444"),
      axis.text.x = element_text(angle = 45, hjust=1) ,
      axis.text.y  = element_text( color="#444444", size=figFontSize-2 ) 
)

cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures/")

aggreg_at_rank <- "Family"

# --- Fig. A - microbiome FAMILIES -----------------------------------------------
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
ASV  <- ASV[,mask]
TAXO <- TAXO[,mask]

# AGGREGATE and reduce ASVs at a specific taxa level ("aggreg_at_rank"), 
taxa_rank <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
aggreg_at_rank <- "Family"
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




# --- Plotting CoV of MICROBIOME -------------------------------------------------------
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
sub_TAXO  <- reduc_TAXO[ mask , c(upper_rank, aggreg_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")
# First, merge reduc_TAXO with ref_rank to add the "Color" column (matching upper_rank ID) to reduc_TAXO
sub_TAXO <- merge( sub_TAXO, ref_ranks, by=upper_rank)
# Second, select only unique elements of "aggreg_at_rank" (since it is in long-form)
# and then merge with the dsPlot, to add Colors and upper_rank to final plotting dataset
temp <- unique( sub_TAXO[, c(aggreg_at_rank,"Colour",upper_rank)])
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
dsPlot$Feature <- factor( dsPlot$Feature, levels=rev(ordered_Features))                      # rev(ordered_Features)
dsPlot[upper_rank] <- factor( dsPlot[[upper_rank]], levels=ref_ranks[[upper_rank]] )
# Now subset to take only the "top" taxa
dsPlot <- dsPlot[dsPlot$Feature %in% ordered_Features[1:25], ]

# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot %>% distinct( .[[upper_rank]], .keep_all=TRUE)
reference_Factors[upper_rank] <- factor( reference_Factors[[upper_rank]], levels=ref_ranks[[upper_rank]] ) 
reference_Factors <- reference_Factors[ order(reference_Factors[[upper_rank]]), ]

# Plot sub-figure A - first layer
my_palette_1 <- reference_Factors$Colour
dsPlot_1 <- dsPlot
hnd_1 <- ggplot( dsPlot_1, aes(x=Feature, y=CoV) ) +  
        geom_boxplot( aes(fill=.data[[upper_rank]]), color="#777777", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
        geom_jitter(  color="#555555", alpha=0.5, size=0.6, width=0.2) +
        coord_flip() + 
        ylim( -0.1, 2.1) + 
        labs(fill="Class:") +
        xlab("") + ylab("CoV") +        ggtitle("CoV - by Family") +
        theme_BoxPlot +
        scale_colour_manual( values = my_palette_1) +
        scale_fill_manual( values = my_palette_1)
hnd_1

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create Legend for Figure 2 and 3
# Remove the features that have CoV at zeros
my_palette_2 <- reference_Factors$Colour
dsPlot_X2 <- dsPlot[ dsPlot$Exp_Name == "A10" | dsPlot$Exp_Name == "A12", ]
hnd_X <- ggplot( dsPlot_X2, aes(x=Feature, y=CoV, fill=.data[[upper_rank]]) ) +  
  geom_col( position="fill", color="#CCCCCC", size=0.25, just=1, width =1) + 
  coord_flip() + 
  theme_BoxPlot + 
  theme( panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_blank(),
         axis.line.x = element_blank() ) + 
  scale_colour_manual( values = my_palette_2) +
  scale_fill_manual( values = my_palette_2)
hnd_X
ggsave( paste0(path_saveFig, "Legend_uBiome_Class2Family.tiff"), width=9, height=15, dpi=300)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Prepare QC of METABOLITES 
# Now lets include the repeated measurements control (CoV_analysis_Control.R) 
# Load ASV table(s) and Metadata
ASV  <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/Reg_All/comb_ASV_Freq.txt",
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) )
TAXO <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/Reg_All/comb_TAXO.txt",
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) ) 
MetaData <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/comb_Metadata.txt",
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

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

# AGGREGATE and reduce ASVs at a specific taxa level ("aggreg_at_rank"), 
taxa_rank <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
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

# Take only features EXP2 datasets
mask <- str_detect( MetaData$External_ID , "Exp2_DNA_" )
mask <- mask + str_detect( MetaData$External_ID , "A09" )    # A01, A06, A09, A16, A20
mask <- mask + str_detect( MetaData$External_ID , "C[123]" )     # C1 to C6    
mask <- mask == 3
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]

sub_MetaData <- MetaData


# Group by individual and calculate the coefficient of variation
cov_val <- cbind(sub_MetaData,  reduc_ASVs) %>% 
            group_by(Exp_Name) %>%
            reframe(  across( colnames(reduc_ASVs), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
            as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)

dsCTRL_1 <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                    names_to  = "Feature",     # selected column for a new variable ("groupings")
                                    values_to = "CoV")         # values all into one new variables

# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsCTRL_1$Feature <- sapply( dsCTRL_1$Feature , function(x)  gsub("\\[", "", x) )
dsCTRL_1$Feature <- sapply( dsCTRL_1$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" 
dsCTRL_1$Feature <- factor( dsCTRL_1$Feature, levels=rev( sort(unique(dsPlot_1$Feature), decreasing=T)) )         
dsCTRL_1 <- dsCTRL_1[ !is.na(dsCTRL_1$Feature) , ]


# Plot sub-figure A - overlay layer 2, controls
hnd_1 <- hnd_1 + #ggplot(NULL) + 
          geom_jitter(  data=dsCTRL_1, aes(x=Feature, y=CoV),
                        color="#8a4949", fill="#963939", shape=23, size=1.5, width=0) +
          coord_flip() + 
          ylim( -0.1, 2.1)+
          theme_BoxPlot +
          scale_colour_manual( values = my_palette_1) +
          scale_fill_manual( values = my_palette_1)
hnd_1





aggreg_at_rank <- "Class"
# --- Fig. B - microbiome CLASSES -----------------------------------------------

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
ASV  <- ASV[,mask]
TAXO <- TAXO[,mask]

# AGGREGATE and reduce ASVs at a specific taxa level ("aggreg_at_rank"), 
taxa_rank <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
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




# --- Plotting CoV of MICROBIOME -------------------------------------------------------
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
sub_TAXO  <- reduc_TAXO[ mask , c(upper_rank, aggreg_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")
# First, merge reduc_TAXO with ref_rank to add the "Color" column (matching upper_rank ID) to reduc_TAXO
sub_TAXO <- merge( sub_TAXO, ref_ranks, by=upper_rank)
# Second, select only unique elements of "aggreg_at_rank" (since it is in long-form)
# and then merge with the dsPlot, to add Colors and upper_rank to final plotting dataset
temp <- unique( sub_TAXO[, c(aggreg_at_rank,"Colour",upper_rank)])
colnames(temp)[1] <- "Feature"
dsPlot <- merge( temp, dsPlot, by="Feature", all=TRUE)
dsPlot <- rename( dsPlot, "Class"="Class.1")
# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsPlot$Feature <- sapply( dsPlot$Feature , function(x)  gsub("\\[", "", x) )
dsPlot$Feature <- sapply( dsPlot$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" and "F_Type" 
dsPlot$Feature <- factor( dsPlot$Feature, levels=rev(ordered_Features))                      # rev(ordered_Features)
dsPlot[upper_rank] <- factor( dsPlot[[upper_rank]], levels=ref_ranks[[upper_rank]] )
# Now subset to take only the "top" taxa
dsPlot <- dsPlot[dsPlot$Feature %in% ordered_Features[1:25], ]

# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot %>% distinct( .[[upper_rank]], .keep_all=TRUE)
reference_Factors[upper_rank] <- factor( reference_Factors[[upper_rank]], levels=ref_ranks[[upper_rank]] ) 
reference_Factors <- reference_Factors[ order(reference_Factors[[upper_rank]]), ]

# Filter out clsses that have a median = 0
temp <- dsPlot %>% group_by(Feature) %>% summarise(median = median(CoV), IQR = IQR(CoV) )
mask <- temp[ temp$median != 0 , "Feature"]
dsPlot <- dsPlot[ dsPlot$Feature %in% unlist(mask) , ]

# Plot sub-figure B - first layer
my_palette_2 <- reference_Factors$Colour
dsPlot_2 <- dsPlot
hnd_2 <- ggplot( dsPlot_2, aes(x=Feature, y=CoV) ) +  
          geom_boxplot( aes(fill=.data[[upper_rank]]), color="#777777", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
          geom_jitter(  color="#555555", alpha=0.5, size=0.6, width=0.2) +
          coord_flip() + 
          ylim( -0.1, 2.1) + 
          labs(fill="Class:") +
          xlab("") + ylab("CoV") +        # ggtitle("CoV - Microbiome") +
          ggtitle("CoV - by Class") +
          theme_BoxPlot +
          scale_colour_manual( values = my_palette_2) +
          scale_fill_manual( values = my_palette_2)
hnd_2


# --- Prepare QC of METABOLITES 
# Now lets include the repeated measurements control (CoV_analysis_Control.R)
ASV  <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/Reg_All/comb_ASV_Freq.txt",
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) )
TAXO <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/Reg_All/comb_TAXO.txt",
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) ) 
MetaData <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/comb_Metadata.txt",
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

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

# AGGREGATE and reduce ASVs at a specific taxa level ("aggreg_at_rank"), 
taxa_rank <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
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

# Take only features EXP2 datasets
mask <- str_detect( MetaData$External_ID , "Exp2_DNA_" )
mask <- mask + str_detect( MetaData$External_ID , "A09" )    # A01, A06, A09, A16, A20
mask <- mask + str_detect( MetaData$External_ID , "C[123]" )     # C1 to C6    
mask <- mask == 3
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]

sub_MetaData <- MetaData


# Group by individual and calculate the coefficient of variation
cov_val <- cbind(sub_MetaData,  reduc_ASVs) %>% 
            group_by(Exp_Name) %>%
            reframe(  across( colnames(reduc_ASVs), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
            as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)

dsCTRL_2 <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                    names_to  = "Feature",     # selected column for a new variable ("groupings")
                                    values_to = "CoV")         # values all into one new variables

# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsCTRL_2$Feature <- sapply( dsCTRL_2$Feature , function(x)  gsub("\\[", "", x) )
dsCTRL_2$Feature <- sapply( dsCTRL_2$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" 
dsCTRL_2$Feature <- factor( dsCTRL_2$Feature, levels=rev( sort(unique(dsPlot_2$Feature), decreasing=T)) )         
dsCTRL_2 <- dsCTRL_2[ !is.na(dsCTRL_2$Feature) , ]

# Plot sub-figure B - overlay layer 2, controls
hnd_2 <- hnd_2 + #ggplot(NULL) + 
  geom_jitter(  data=dsCTRL_2, aes(x=Feature, y=CoV),
                color="#8a4949", fill="#963939", shape=23, size=1.5, width=0) +
  coord_flip() + 
  ylim( -0.1, 2.1)+
  theme_BoxPlot +
  scale_colour_manual( values = my_palette_2) +
  scale_fill_manual( values = my_palette_2)

hnd_2





aggreg_at_rank <- "Genus"
# --- Fig. B - microbiome GENUS -----------------------------------------------
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
ASV  <- ASV[,mask]
TAXO <- TAXO[,mask]

# AGGREGATE and reduce ASVs at a specific taxa level ("aggreg_at_rank"), 
taxa_rank <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
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



sub_MetaData <- MetaData

# Prepare CoV values
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
sub_TAXO  <- reduc_TAXO[ mask , c(upper_rank, aggreg_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")
# First, merge reduc_TAXO with ref_rank to add the "Color" column (matching upper_rank ID) to reduc_TAXO
sub_TAXO <- merge( sub_TAXO, ref_ranks, by=upper_rank)
# Second, select only unique elements of "aggreg_at_rank" (since it is in long-form)
# and then merge with the dsPlot, to add Colors and upper_rank to final plotting dataset
temp <- unique( sub_TAXO[, c(aggreg_at_rank,"Colour",upper_rank)])
colnames(temp)[1] <- "Feature"
dsPlot <- merge( temp, dsPlot, by="Feature", all=TRUE)
# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsPlot$Feature <- sapply( dsPlot$Feature , function(x)  gsub("\\[", "", x) )
dsPlot$Feature <- sapply( dsPlot$Feature , function(x)  gsub("\\]", "", x) )


# We try to also merge the information regarding all the taxonomic levels onto dsPlot
taxonomics <- t(TAXO)
colnames(taxonomics) <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Feature", "Species")
rownames(taxonomics) <- seq(1, length(rownames(taxonomics) ))
taxonomics <- as.data.frame(taxonomics)
# Select only taxonomy that that is present as "Feature" in dsPlot
unique_features <- unique(dsPlot$Feature)   
taxonomics <- taxonomics[ unique_features %in% taxonomics$Feature, ]
taxonomics <- taxonomics %>% distinct( Feature, .keep_all = TRUE)
# Merge the taxonomy to dsPlot
dsPlot <- merge( dsPlot, taxonomics, by="Feature")   # inner_join( dsPlot, temp, by=c("Feature"="Genus", "Phyla"="Phyla")) 
# dsPlot <- select( dsPlot, -c("Class.x") )
colnames(dsPlot)[colnames(dsPlot) == "Class.x"] = "Class"
# Select only features that belong to one specific Class
unique(dsPlot$Class)


# --- plot all Classes together -----------------------------------------------
# We will plot all the top 20 genera for each class, all in one plot. However, 
# we will also sort taxa to be ordered first by class, and then by abundancy.
list_uppper_taxa <- c( "Clostridia", "Bacteroidia","Erysipelotrichia","Actinobacteria",
                       "Coriobacteriia","Negativicutes","Verrucomicrobiae","Mollicutes",
                       "Deltaproteobacteria","Gammaproteobacteria","Bacilli") 
sub_dsPlot <- data.frame()                    #Store values of all feature to include in chart
reference_Factors <- data.frame()             # store the class (upper_taxa) order as we loop
loop_lower_taxa <- list()                    # store the genera (lower_taxa) order as we loop
for ( choosen_taxa in list_uppper_taxa){
  temp1 <- dsPlot[ dsPlot$Class == choosen_taxa , ]
  # subset and take only top 20 genera present in the class
  subset_taxa_names <- intersect(ordered_Features, unique(temp1$Feature))[1:20]
  subset_taxa_names <- subset_taxa_names[ !is.na(subset_taxa_names) ]
  loop_lower_taxa <- append( loop_lower_taxa, subset_taxa_names )
  temp1 <- temp1[  temp1$Feature %in% subset_taxa_names , ]
  sub_dsPlot <- rbind( sub_dsPlot , temp1 )
  # Create a reference table that set the order of (unique) F_Type, and then we 
  # can take the (unique) Colors as list ordered to match the F-Type
  temp2 <- temp1 %>% distinct( .[[upper_rank]], .keep_all=TRUE)
  reference_Factors <- rbind( reference_Factors, temp2 )
}
# Factorize the "Feature" and "F_Type" to the loooping order (list_uppper_taxa and loop_lower_taxa)
sub_dsPlot$Feature <- factor( sub_dsPlot$Feature, levels=rev(unlist(loop_lower_taxa)))                      # rev(ordered_Features)
sub_dsPlot[upper_rank] <- factor( sub_dsPlot[[upper_rank]], levels=list_uppper_taxa )

# Filter out clsses that have a median = 0
temp <- sub_dsPlot %>% group_by(Feature) %>% summarise(median = median(CoV), IQR = IQR(CoV) )
mask <- temp[ temp$median != 0 , "Feature"]
sub_dsPlot <- sub_dsPlot[ sub_dsPlot$Feature %in% unlist(mask) , ]


# Set Plotting theme and colors
my_palette_3 <- reference_Factors$Colour
figFontSize <- 11

scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = my_palette_3)
}
scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = my_palette_3)
}


# Plot sub-figure C - layer 1
sub_dsPlot_1 <- sub_dsPlot
hnd_3 <- ggplot( sub_dsPlot_1, aes(x=Feature, y=CoV) ) +  
      geom_boxplot( aes(fill=.data[[upper_rank]]), color="#666666", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
      geom_jitter(  color="#666666", alpha=0.5, size=0.5, width=0.2) +
      coord_flip() + 
      ylim( -0.1, 2.1) +
      ggtitle("CoV - by Genus") +
      labs(fill="Class:") +
      xlab("") + ylab("CoV") +
      theme_BoxPlot
hnd_3

# --- Prepare QC of METABOLITES 
# Now lets include the repeated measurements control (CoV_analysis_Control.R)
# Load ASV table(s) and Metadata
ASV  <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/Reg_All/comb_ASV_Freq.txt",
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) )
TAXO <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/Reg_All/comb_TAXO.txt",
                  skip_empty_rows = TRUE ,show_col_types = FALSE,
                  guess_max  = min( 10 , 300 ) ) 
MetaData <- read_tsv( "/Users/sm/Documents/R/16SMT_ctrl/DataTables/comb_Metadata.txt",
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

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

# Take only features EXP2 datasets
mask <- str_detect( MetaData$External_ID , "Exp2_DNA_" )
mask <- mask + str_detect( MetaData$External_ID , "A09" )    # A01, A06, A09, A16, A20
mask <- mask + str_detect( MetaData$External_ID , "C[123]" )     # C1 to C6    
mask <- mask == 3
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]

sub_MetaData <- MetaData


# Group by individual and calculate the coefficient of variation
cov_val <- cbind(sub_MetaData,  reduc_ASVs) %>% 
  group_by(Exp_Name) %>%
  reframe(  across( colnames(reduc_ASVs), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
  as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)

dsCTRL_3 <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                      names_to  = "Feature",     # selected column for a new variable ("groupings")
                                      values_to = "CoV")         # values all into one new variables

# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsCTRL_3$Feature <- sapply( dsCTRL_3$Feature , function(x)  gsub("\\[", "", x) )
dsCTRL_3$Feature <- sapply( dsCTRL_3$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" 
dsCTRL_3$Feature <- factor( dsCTRL_3$Feature, levels=rev( sort(unique(sub_dsPlot_1$Feature), decreasing=T)) )         
dsCTRL_3 <- dsCTRL_3[ !is.na(dsCTRL_3$Feature) , ]

# Plot sub-figure C - overlay layer 2, controls
hnd_3 <- hnd_3 + #ggplot(NULL) + 
  geom_jitter(  data=dsCTRL_3, aes(x=Feature, y=CoV),
                color="#8a4949", fill="#963939", shape=23, size=1.5, width=0) +
  coord_flip() + 
  ylim( -0.1, 2.1)+
  theme_BoxPlot +
  scale_colour_manual( values = my_palette_3) +
  scale_fill_manual( values = my_palette_3)

hnd_3




# --- FIGURE - combine sub-plots -------------------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")

ggdraw() +
  draw_plot(hnd_2, x = 0.05,  y = 0.69,  width = .40, height = .31) +
  draw_plot(hnd_1, x = 0.0,  y = 0.12,  width = .45, height = .56) +
  draw_plot(hnd_3, x = .5,    y = 0,    width = .47, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 17, x = c( .08, .08, .54), y = c(1,0.68,1))

ggsave( paste0(path_saveFig, "Suppl_Figure_2.tiff"), width=10, height=12, dpi=300)
  



























