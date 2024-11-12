# +++++ DESCRIPTION - Figure_Suppl_1 +++++++++++++++++++++++++++++
#
# Create Supplementary Figure 2, which are the Alpha (Shannon index) and 
# Beta diversity (Jaccard Similarity) of microbial taxa grouped by individuals. 
# There are two subfigures:
# - SF 1 A - Shannon index          --> hnd_3
# - SF 1 B - Jaccard Similarity     --> hnd_4
#
# Essentially, this scripts is based on:
# - Metrics_microbiome_AlphaBeta.R
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2023
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# --- Load Libraries
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)

cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures/")

# Import custom functions, specifically contains those to calculated microbiome diversity metrics
source( paste0(cDir_Rscript, "/SangerMicrobiome.R") )

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
        legend.background = element_blank() ,
        legend.title = element_text(size = figFontSize), 
        legend.text  = element_text(size = figFontSize-2),
        axis.line.x = element_line(size = 1, colour = "#4e4e4e"),
        axis.text.x = element_text(angle = 45, hjust=1) ,
        axis.text.y  = element_text( color="#444444", size=figFontSize-2 ) 
)


# ---- Load DATA tables
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


# --- Pruning the data.frame ---------------------------------------------------
# Sample_ID is the sequencing ID. First, order (ascent) the two tables by
# their shared column (Sample_ID). Then, we can copy paste the Exp_Name and 
# guarantees respect match row order.
MetaData <- MetaData[ order( MetaData$Sample_ID ), ]
ASV <- ASV[ order( ASV$Sample_ID ), ]
ASV <- mutate( ASV, Exp_Name = MetaData[ order(ASV$Sample_ID), "Exp_Name"] [[1]], .after=Sample_ID ) 

# --- Set the DF to use for plotting
# Create Datasets based on TimePoint
temp_Meta <- filter(MetaData, Exp_TimePoint == "0")
tab_00 <- ASV[ ASV$Sample_ID %in% temp_Meta$Sample_ID , ]
tab_00 <- tab_00[ order( tab_00$Exp_Name ), ]
temp_Meta <- filter(MetaData, Exp_TimePoint == "7")
tab_01 <- ASV[ ASV$Sample_ID %in% temp_Meta$Sample_ID , ]
tab_01 <- tab_01[ order( tab_01$Exp_Name ), ]
temp_Meta <- filter(MetaData, Exp_TimePoint == "14")
tab_02 <- ASV[ ASV$Sample_ID %in% temp_Meta$Sample_ID , ]
tab_02 <- tab_02[ order( tab_02$Exp_Name ), ]
temp_Meta <- filter(MetaData, Exp_TimePoint == "28")
tab_03 <- ASV[ ASV$Sample_ID %in% temp_Meta$Sample_ID , ]
tab_03 <- tab_03[ order( tab_03$Exp_Name ), ]

temp_names <- tab_00$Exp_Name        # Rename and remove not numeric data
tab_00 <- tab_00[ , -c(1,2)]                  
rownames(tab_00) <- temp_names
temp_names <- tab_01$Exp_Name        # Rename and remove
tab_01 <- tab_01[ , -c(1,2)]
rownames(tab_01) <- temp_names
temp_names <- tab_02$Exp_Name        # Rename and remove
tab_02 <- tab_02[ , -c(1,2)] 
rownames(tab_02) <- temp_names
temp_names <- tab_03$Exp_Name        # Rename and remove
tab_03 <- tab_03[ , -c(1,2) ] 
rownames(tab_03) <- temp_names

# (Optional) Check
rownames(tab_00) == rownames(tab_01) 
rownames(tab_02) == rownames(tab_03)



# ---- A - SHANNON diversity -------------------------------------------------------
# Calculate the values to use for plotting and Transform into long format DF
Shan <- NULL
Shan <- data.frame( Shan_00 = apply( tab_00, 1, function(x) met_Shannon( x ) ) ,
                    Shan_07 = apply( tab_01, 1, function(x) met_Shannon( x ) ) ,
                    Shan_14 = apply( tab_02, 1, function(x) met_Shannon( x ) ) ,
                    Shan_28 = apply( tab_03, 1, function(x) met_Shannon( x ) ) )

Shan <- as.data.frame( t(Shan) )
colnames(Shan) <- lapply( sprintf('%0.2d', 1:14), function(x) paste0("A",x))
Shan$RowID <- c("d00", "d07", "d14", "d28") 
dfPlot <- Shan %>% pivot_longer( cols = starts_with("A"),            # Aggregate all columns starting with "d"
                                 names_to  = "Shannon",
                                 values_to = "values" )
dfPlot$Shannon <- factor(dfPlot$Shannon, levels=colnames(Shan) )     # Factorize x-axis to display in temporal order

# Create a boxplot, which also connects with lines matching observations
hnd_1 <- ggplot( dfPlot, aes(x=Shannon, y=values, fill=Shannon)) + 
          geom_boxplot( fill="#2593BD", color="#666666", alpha=0.7, outlier.shape=NA) + 
          geom_jitter(  aes(fill=Shannon), color="#666666", alpha=0.5, size=1.2, width=0.2) +
          scale_y_continuous( expand = c(0, NA), limits = c(2.5, 3.75)) +
          theme(legend.position = "none") + 
          xlab("")  +  ylab("Shannon index") +
          theme_BoxPlot 
hnd_1 


# ---- B - JACCARD diversity -------------------------------------------------------
# Calculate metric index for both tables, comparing the timepoints. 
Jacc <- data.frame( Jacc_01 = met_Jaccard_Similarity_2DF( tab_00, tab_01 ) ,
                    Jacc_02 = met_Jaccard_Similarity_2DF( tab_00, tab_02 ) ,
                    Jacc_03 = met_Jaccard_Similarity_2DF( tab_00, tab_03 ) )
Jacc <- as.data.frame( t(Jacc) )
# REMOVE the problematic timepoint of day 28 for individual A12; reads <1000 for this sample
Jacc[3,"A12"] <- NaN
colnames(Jacc) <- lapply( sprintf('%0.2d', 1:14), function(x) paste0("A",x))
Jacc$RowID <- c( "d.0 - d.7", "d.0 - d.14", "d.0 - d.28")  
dfPlot <- Jacc %>% pivot_longer( cols = starts_with("A"),
                                 names_to  = "Jaccard",
                                 values_to = "values" )

hnd_2 <- ggplot( dfPlot, aes(x=Jaccard, y=values, fill=Jaccard)) + 
  geom_boxplot( fill="#F8BD39", color="#666666", alpha=0.7, outlier.shape=NA) +  
  geom_jitter(  aes(fill=Jaccard), color="#666666", alpha=0.5, size=1.2, width=0.2 ) +
  scale_y_continuous( expand = c(0, NA), limits = c(0, 0.8)) +
  theme(legend.position = "none") + 
  xlab("")  +  ylab("Jaccard Similarity") +
  theme_BoxPlot
hnd_2 


# --- FIGURE - Combine sub-plots -----------------------------------------------
# Place sub-plots together into one figure
library("cowplot")

ggdraw() +
  draw_plot(hnd_1, x = .20, y = .55, width = .60, height = .4) +
  draw_plot(hnd_2, x = .20, y = .07, width = .60, height = .4) +
  draw_plot_label(label = c("A", "B"), size = 17, x = c(.18,.18), y = c(.98, 0.50))

ggsave( paste0(path_saveFig, "Suppl_Figure_1.tiff"), width=10, height=10, dpi=300)



