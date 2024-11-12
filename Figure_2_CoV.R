# +++++ DESCRIPTION - Figure_2_CoV.R +++++++++++++++++++++++++++++++++++++++++++
#
# Create the 2nd figure of the Longitudinal paper using this script. There are 
# two subfigures:
# - Fig 2 A - CoV of microbial Genera       --> hnd_A
# - Fig 2 B - CoV of metabolites            --> hnd_B
#
# This figures summarise the results from the analysis of the Coefficient of
# Variation (CoV) of microbiome and metabolites data. CoV analysis is done
# separately for the microbiome genera (Fig. 2A) and metabolites (Fig. 2B).
#
# Lastly, at section Combine all the "handles" of the 2 subplots (hnd_X) are
# arranged into one publish read figure. 
# NOTE: legends are not plotted.
# The taxa and metabolites are colored consistently across the paper, using two 
# file that rank and assign specific color tone and saturation, so that features 
# from the same class share the same tone and the plots appear better structured 
# and easier to understand
#
# Essentially this combines the previous scripts:
# - CoV_analysis.R
# - CoV_analysis_Control.R
# We plot the CoV of features, either taxa or metabolites, so the section that 
# calculatese CoVs grouping samples by individuals has been removed.
# (This we will put in a separate Suppl. Figure)
# 
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



# ---- Load Data tables --------------------------------------------------------------
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

# --- Pruning the two datasets
aggreg_at_rank <- "Genus"
source(paste0( path_2_data, "pruning_2Datasets.R"))



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Fig. A - METABOLITES ------------------------------------------------------
# Calculate CoV - [OPT 1] - Group by individual and calculate the coefficient of variation
cov_val <- cbind(sub_MetaData,  MB_normal) %>% 
             group_by(Exp_Name) %>%
             reframe(  across( colnames(MB_normal), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
             as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)

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
dsPlot$Feature <- factor( dsPlot$Feature, levels=rev(metabo_Ranking$Feature) ) 
dsPlot$F_Type  <- factor( dsPlot$F_Type, levels= unique(sort(dsPlot$F_Type)) )  
# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot %>% distinct( F_Type, .keep_all=TRUE)
reference_Factors <- reference_Factors[ order(reference_Factors$F_Type) , ]
 
dsPlot <- dsPlot[ !is.na(dsPlot$CoV), ] # drop NA only rows

# Plot sub-figure A - first layer
my_palette_1 <- reference_Factors$Colour
dsPlot_1 <- dsPlot
hnd_A <- ggplot( dsPlot_1, aes(x=Feature, y=CoV)) + 
          geom_boxplot( aes(fill=F_Type), color="#777777", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
          geom_jitter( color="#555555", alpha=0.5, size=0.6, width=0.2) +
          coord_flip() +
          ylim( -0.1, 2.1) +  
          labs(fill="Metabolites Type:") +
          xlab("") + ylab("CoV") +  
          theme_BoxPlot +
          scale_colour_manual( values = my_palette_1) +
          scale_fill_manual( values = my_palette_1)
hnd_A


# --- Prepare QC of METABOLITES 
# Now lets include the repeated measurements control (CoV_analysis_Control.R)
ctrl_metabo <- read_tsv( paste0(path_2_data, "Data_NMR/Ctrl_Repeated_Measurments.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

Sample_ID   <- data.frame( "Exp_Name"=lapply( ctrl_metabo[,1], function(x) substring(x,1,1) ))
ctrl_metabo <- ctrl_metabo[,-1]
cov_val <- cbind( Sample_ID, ctrl_metabo) %>% 
            group_by(Exp_Name) %>%
            reframe(  across( colnames(ctrl_metabo), function(x) sd(x, na.rm=T)/mean(x, na.rm=T) ) ) %>%
            as.data.frame()
cov_val <- replace(cov_val, is.na(cov_val), 0)
cov_val <- cov_val[ cov_val$Exp_Name=="H" , ]

dsCTRL_1 <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                    names_to  = "Feature",       # selected column for a new variable ("groupings")
                                    values_to = "CoV")           # values all into one new variables

# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(MB_normal)[ order(colSums(MB_normal), decreasing=T) ]
# Factorize the "Feature" 
dsCTRL_1$Feature <- factor( dsCTRL_1$Feature, levels=rev( sort(unique(dsPlot_1$Feature), decreasing=T)) )         
dsCTRL_1 <- dsCTRL_1[ !is.na(dsCTRL_1$Feature) , ]


# Plot sub-figure A - overlay layer 2, controls
hnd_A <- hnd_A + #ggplot(NULL) + 
          geom_jitter( data=dsCTRL_1, aes(x=Feature, y=CoV), 
                       color="#AA4949", fill="#963939", shape=23, size=1.5, width=0) +
          coord_flip() + 
          ylim( -0.1, 2.1)+
          theme_BoxPlot +
          scale_colour_manual( values = my_palette_1) +
          scale_fill_manual( values = my_palette_1)

hnd_A



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Fig. B - MICROBIOME -----------------------------------------------
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
dsPlot$Feature <- factor( dsPlot$Feature, levels=rev(ordered_Features))      
dsPlot[upper_rank] <- factor( dsPlot[[upper_rank]], levels=ref_ranks[[upper_rank]] )
# Now subset to take only the "top" taxa
dsPlot <- dsPlot[dsPlot$Feature %in% ordered_Features[1:51], ]

# Remove the features that have CoV at zeros
dsPlot_X <- dsPlot
dsPlot <- dsPlot[ dsPlot$Feature != "Megamonas", ]  

# Create a reference table that set the order of (unique) F_Type, and then we 
# can take the (unique) Colors as list ordered to match the F-Type
reference_Factors <- dsPlot %>% distinct( .[[upper_rank]], .keep_all=TRUE)
reference_Factors[upper_rank] <- factor( reference_Factors[[upper_rank]], levels=ref_ranks[[upper_rank]] ) 
reference_Factors <- reference_Factors[ order(reference_Factors[[upper_rank]]), ]

# Plot sub-figure A - first layer, CoV distributions
my_palette_2 <- reference_Factors$Colour
dsPlot_2 <- dsPlot
hnd_B <- ggplot( dsPlot_2, aes(x=Feature, y=CoV) ) +  
        geom_boxplot( aes(fill=.data[[upper_rank]]), color="#777777", linewidth=0.35, alpha=0.6, outlier.shape=NA) +  
        geom_jitter(  color="#555555", alpha=0.5, size=0.6, width=0.2) +
        # stat_summary(fun.y=mean, geom="point", shape=18, size=2, color="grey20", fill="grey20") +
        coord_flip() + 
        ylim( -0.1, 2.1) + 
        labs(fill="Class:") +
        xlab("") + ylab("CoV") +        # ggtitle("CoV - Microbiome") +
        theme_BoxPlot +
        scale_colour_manual( values = my_palette_2) +
        scale_fill_manual( values = my_palette_2)
hnd_B

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create Legend for Figure 2 and 3
# Remove the features that have CoV at zeros
my_palette_2 <- reference_Factors$Colour
dsPlot_X2 <- dsPlot_X[ dsPlot_X$Exp_Name == "A10" | dsPlot_X$Exp_Name == "A12", ]
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

ggsave( paste0(path_saveFig, "Legend_uBiome_Class2Genus.tiff"), width=8, height=15, dpi=300)
 # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Prepare QC of Microbiome 
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
mask <- mask + str_detect( MetaData$External_ID , "A09" )
mask <- mask + str_detect( MetaData$External_ID , "C[123]" )    
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

dsCTRL <- cov_val %>% pivot_longer( cols      = -Exp_Name,     # exclude specific column
                                    names_to  = "Feature",     # selected column for a new variable ("groupings")
                                    values_to = "CoV")         # values all into one new variables

# Order the Features to be plotted by the relative frequency in which they appear
# Also, adjust the Features string, removing unwanted characters
ordered_Features <- colnames(reduc_ASVs)[ order(colSums(reduc_ASVs), decreasing=T) ]
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\[", "", x) )
ordered_Features <- sapply( ordered_Features , function(x)  gsub("\\]", "", x) )
dsCTRL$Feature <- sapply( dsCTRL$Feature , function(x)  gsub("\\[", "", x) )
dsCTRL$Feature <- sapply( dsCTRL$Feature , function(x)  gsub("\\]", "", x) )

# Factorize the "Feature" 
dsCTRL$Feature <- factor( dsCTRL$Feature, levels=rev( sort(unique(dsPlot_2$Feature), decreasing=T)) )         
dsCTRL <- dsCTRL[ !is.na(dsCTRL$Feature) , ]


# Plot sub-figure B - overlay layer 2, controls
hnd_B <- hnd_B + #ggplot(NULL) + 
          geom_jitter(  data=dsCTRL, aes(x=Feature, y=CoV),
                        color="#8a4949", fill="#963939", shape=23, size=1.5, width=0) +
          coord_flip() + 
          ylim( -0.1, 2.1)+
          theme_BoxPlot +
          scale_colour_manual( values = my_palette_2) +
          scale_fill_manual( values = my_palette_2)
hnd_B



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- COMBINE sub-Figures ------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")
ggdraw()

ggdraw() +
  draw_plot(hnd_B, x = 0, y = .28, width = .55, height = .71) +
  draw_plot(hnd_A, x = 0.55, y = 0, width = .45, height = .99) +
  draw_plot_label(label = c("A", "B"), size = 17, x = c(0.07, 0.58), y = c(1, 1))

ggsave( paste0(path_saveFig, "Figure_2.tiff"), width=10, height=12.5, dpi=300)




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Obtain numeric values of CoV mean, median amd IQR that were charted in box-plot
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# The CoV is calculated grouping by time, one for each individual (and feature)
# Now calculate the average CoV of each feature, for the cohort (group_by ID)
CoV_Table_AllFeat <- dsPlot %>% group_by(Feature, Colour) %>% summarise(Average = mean(CoV), Median = median(CoV), IQR = IQR(CoV) )
CoV_Table_AllFeat <- CoV_Table_AllFeat %>% select(Feature,Average,Median,IQR) %>% 
                     distinct( ., Feature, .keep_all = T)
# Sort by feature names
microbiome_Ranking <- read_tsv( paste0(path_2_data, "/Res_Tables/uBiome_RelAbu_Genus_byCohort.txt"),
                                col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                                guess_max  = min( 10 , 300 ) )
sort_microbeRank  <- microbiome_Ranking$Genus
sort_microbeRank  <- gsub("\\[", "", sort_microbeRank) 
sort_microbeRank  <- gsub("\\]", "", sort_microbeRank) 
CoV_Table_AllFeat <- CoV_Table_AllFeat[ match(sort_microbeRank, CoV_Table_AllFeat$Feature ) , ] 


# --- TABLE of CoV mean values and IQR -----------------------------------------
# RUN the code below to calculate the average CoV values for each plotted category in the boxplot
CoV_vals_metabo <- dsPlot_1 %>% group_by(Feature, F_Type, Colour) %>% summarise(Average = mean(CoV), Median = median(CoV), IQR = IQR(CoV) )
CoV_vals_taxa   <- dsPlot_2 %>% group_by(Feature, Class, Colour)  %>% summarise(Average = mean(CoV), Median = median(CoV), IQR = IQR(CoV) )

# Since we "flip" the coordinates plotting, the features are in reverse order in the dataframe
# Flip the dataframes to have them sorted as they appear in the plot.
CoV_vals_metabo <- CoV_vals_metabo[seq(dim(CoV_vals_metabo)[1],1) , ]
CoV_vals_taxa   <- CoV_vals_taxa[seq(dim(CoV_vals_taxa)[1],1) , ]

# Save the data in .txt files; it is already sorted by abundacy
write_tsv( CoV_vals_metabo, paste0( path_2_data, "Res_Tables/Table_CoV_Metabolites.txt") )
write_tsv( CoV_vals_taxa, paste0( path_2_data, "Res_Tables/Table_CoV_uB_",aggreg_at_rank,".txt") )


# Quick test to find statistical differences in in distribution of means or median values of features
t.test( CoV_vals_metabo$Average, CoV_vals_taxa$Average )
t.test( CoV_vals_metabo$Median, CoV_vals_taxa$Median )


