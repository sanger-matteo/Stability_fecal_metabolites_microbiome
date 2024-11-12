# +++++ DESCRIPTION - Figure_Suppl_5.R ++++++++++++++++++++++++++++++++++++++
#
# Create Supplementary Figure 5, which plots the features (of either microbiome
# or metabolites) in a scatter plot, as IQR versus Median of the calculated CoV
# values. 
# There are two subfigures:
# - SF 5 A - IQR-vs-Median  of microbial (genera)      --> hnd_A
# - SF 5 B - IQR-vs-Median of metabolites              --> hnd_B
#
# Essentially, this scripts is based on:
# - Scatter_IQR-vs-CoV.R
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
  panel.background = element_rect(size = 2, colour = "#9e9e9e", linetype = "solid", fill=NA) ,
  plot.title   = element_text( color="#222222", size=figFontSize+4, face = "bold", hjust = .5) ,
  axis.title.x = element_text( color="#444444", size=figFontSize, vjust = -.6, hjust = .5) ,
  axis.title.y = element_text( color="#444444", size=figFontSize, vjust =  .6, hjust = .5) ,
  axis.text.x  = element_text( color="#444444", size=figFontSize-3,   vjust = -1), 
  axis.text.y  = element_text( color="#444444", size=figFontSize-3 ),
  axis.ticks  = element_blank(),
  panel.grid.major.x = element_line( size = 0.5, colour = "#dddddd", linetype="dashed") ,
  panel.grid.major.y = element_line( size = 0.5, colour = "#dddddd", linetype="dashed") ,
  panel.grid.minor = element_blank() ,
  legend.position = "none",
  legend.key   = element_blank(),
  legend.title = element_text(size = figFontSize-2), 
  legend.text  = element_text(size = figFontSize-2)
)




# ---- Load DATA tables --------------------------------------------------------------
cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures/")

aggreg_at_rank <- "Genus"

# Load ASV table(s) and Metadata
df_uBiome <- read_tsv( paste0( path_2_data, "Res_Tables/Table_CoV_uB_",aggreg_at_rank,".txt"),
                       skip_empty_rows = TRUE ,show_col_types = FALSE,
                       guess_max  = min( 10 , 300 ) ) 
df_metabo <- read_tsv( paste0( path_2_data, "Res_Tables/Table_CoV_Metabolites.txt"),
                       skip_empty_rows = TRUE ,show_col_types = FALSE,
                       guess_max  = min( 10 , 300 ) )


uBiome_RelAbu <- read_tsv( paste0( path_2_data, "Res_Tables/uBiome_RelAbu_Genus_byCohort.txt"),
                           skip_empty_rows = TRUE ,show_col_types = FALSE,
                           guess_max  = min( 10 , 300 ) )
Table_Avg_Metabolites  <- read_tsv( paste0( path_2_data, "Res_Tables/Table_Avg_Metabolites.txt"),
                                    skip_empty_rows = TRUE ,show_col_types = FALSE,
                                    guess_max  = min( 10 , 300 ) ) 

# Give unique columns names that are manageable in R
colnames(uBiome_RelAbu) <- c( "Genus", "Frequency", "Average_RA", "StDev_RA", "Median_RA", "IQR_RA",
                              "Phyla", "Class", "Order", "Family" )
colnames(Table_Avg_Metabolites) <- c( "Metabolite", "Average_RA", "StDev_RA", "Median_RA", "IQR_RA", "Class")

colnames(df_uBiome) <- c( "Feature", "Class",  "Colour", "Average_CoV", "Median_CoV", "IQR_CoV" )
colnames(df_metabo) <- c( "Feature", "F_Type", "Colour", "Average_CoV", "Median_CoV", "IQR_CoV" )

# --- PRUNE the data
# to avoid confusion, drop the column "Average" because it refers to CoVs. 
# Later,we mearge with the other tables that also contains a column "Average", 
# which refers to the Rel. Abundance and it is the one we need to plot.
df_uBiome <- subset( df_uBiome, select=-c(Average_CoV))
df_metabo <- subset( df_metabo, select=-c(Average_CoV))
uBiome_RelAbu <- subset( uBiome_RelAbu, select=-c(Class))
Table_Avg_Metabolites <- subset( Table_Avg_Metabolites, select=-c(Class))

# Adjust the microbiome features names
uBiome_RelAbu$Genus <- sapply( uBiome_RelAbu$Genus , function(x)  gsub("\\[", "", x) )
uBiome_RelAbu$Genus <- sapply( uBiome_RelAbu$Genus , function(x)  gsub("\\]", "", x) )

# Change the column to ensure that the linkage column in both tables have same 
# name and then merge the two tables
colnames(df_uBiome)[1] <- "Genus"
colnames(df_metabo)[1] <- "Metabolite"
df_uBiome <- merge( df_uBiome, uBiome_RelAbu, by="Genus")
df_metabo <- merge( df_metabo, Table_Avg_Metabolites, by="Metabolite")


# --- TAXA IQR vs Median CoVs --------------------------------------------------
# Set Plotting theme and colors
df_uBiome$Class <- factor( df_uBiome$Class, levels=unique(df_uBiome$Class))
my_palette <- df_uBiome$Colour
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = my_palette)
}
scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = my_palette)
}

hnd_A <- ggplot( df_uBiome, aes(x=Median_CoV, y=Average_RA, color=Colour )) + 
        geom_point( alpha = 0.75, size = 2.5  ) +
        scale_color_identity()  + 
        ggtitle("Microbiome Genera") +
        xlab("median CoV") + ylab("Mean Rel. Abundance (%)") +
        theme_BoxPlot 
hnd_A

# --- METABO IQR vs Median CoVs --------------------------------------------------
# Set Plotting theme and colors
df_metabo$F_Type <- factor( df_metabo$F_Type, levels=unique(df_metabo$F_Type))
my_palette <- df_metabo$Colour
scale_colour_discrete <- function(...) {
  scale_colour_manual(..., values = my_palette)
}
scale_fill_discrete <- function(...) {
  scale_fill_manual(..., values = my_palette)
}

hnd_B <- ggplot( df_metabo, aes(x=Median_CoV, y=Average_RA, color=Colour )) + 
          geom_point( alpha = 0.75, size = 2.5  ) +
          scale_color_identity()  + 
          ggtitle("Metabolites") +
          xlab("median CoV") + ylab("Mean Rel. Abundance (%)") +
          theme_BoxPlot 
hnd_B


# --- COMBINE sub-Figures ------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")
ggdraw()

ggdraw() +
  draw_plot(hnd_A, x = 0, y = 0, width = .48, height =.95) +
  draw_plot(hnd_B, x = .51, y = 0, width = .48, height = .95) +
  draw_plot_label(label = c("A", "B"), size = 17, x = c(0, 0.5), y = c(.98, .98))

ggsave( paste0(path_saveFig, "Suppl_Figure_5_AB.tiff"), width=10, height=5, dpi=300)







