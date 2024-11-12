# +++ Description +++ Figure_3_PartialSpearman +++++++++++++++++++++++++++++++++++
# 
# Create part of the 3nd figure of the Longitudinal paper using this script. 
# The script plots the subfigures:
# - SubF 3 A - Partial Spearman Correlation, which is a correlation that can be 
#              adjusted by specific grouping factors. Then it can be plotted as
#              a heat map. The correlations is between metabolites and microbe 
#              at a specified taxonomic rank, where non-significant correlations
#              (those with nominal p<0.05) are not showed.
#
# Moreover, the script perform p-value adjustment (Bonferroni or BH), and 
# creates a plot that mirrors the correlations plot is between metabolites and 
# microbe. Instead of the correlation values, it simply places a star "*" where 
# the correlation remain significant (p<0.05) after adjustment. The two images 
# can then be overlaid outside R to generate one heat plot that shows correlation
# with two levels of p-values (nominal and adjusted)
#
# NOTE 1: why this convoluted approach? Because there is no other way to make 
#         corrplot filter by p-value and show the symbol using two different 
#         p-value tables (nominal and adjusted)
# NOTE 2: There is no easy way to overlay two images either. When re-importing,  
#         the transparency causes difficulties. It is simpler to combine them 
#         manually in A.Illustrator or PowerPoint.
#
# Essentially, this is the script:
# - Corr_Spearman.R
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2024
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




# --- Load Libraries -----------------------------------------------------------
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(PResiduals)
library(tiff)

# General plotting setting
figFontSize <- 10
theme_HeatPlot <- theme(
        plot.margin = margin(t=0.3, r=0.5, b=0.5, l=0.3, "cm"),
        panel.background = element_blank() ,
        panel.grid.minor = element_blank() ,
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_blank() ,
        plot.title  = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
        axis.title  = element_text( color="#444444", size=figFontSize+1 ) ,
        axis.text   = element_text( color="#444444", size=figFontSize   ) ,
        axis.ticks  = element_blank(),
        legend.position = "none",
        legend.background = element_blank() ,
        legend.title = element_text(size = figFontSize), 
        legend.text  = element_text(size = figFontSize-2),
        axis.text.x  = element_text(color="#444444", size=figFontSize, angle = 90, hjust=0) ,
        axis.text.y  = element_text( color="#444444", size=figFontSize ) 
)


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



# --- Pruning MICROBIOME data --------------------------------------------------
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



# --- Pruning METABOLITES data -------------------------------------------------
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


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Perform Partial Spearman Correlation ---------------------------------------------
# Perform the correlation and the significance (nominal p-value)

# Bring in the ordered names of Genus and Metabolites, which we use to set
# consistent order of row-cols in the heatplot
metabo_Ranking <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Ranking.txt"),
                            col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                            guess_max  = min( 10 , 300 ) )
microbiome_Ranking <- read_tsv( paste0(path_2_data, "/Res_Tables/uBiome_RelAbu_",aggreg_at_rank,"_byCohort.txt"),
                                col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                                guess_max  = min( 10 , 300 ) )

metabo_Ranking  <- as.data.frame(metabo_Ranking)
sort_metaboRank <- metabo_Ranking$Feature_Name
# Subset the data to limit the number of microbiome features to cross-correlate
if(aggreg_at_rank == "Genus") {
  sort_microbeRank <- microbiome_Ranking[[aggreg_at_rank]][ 1:51 ]
}else if (aggreg_at_rank == "Family"){
  sort_microbeRank <- microbiome_Ranking[[aggreg_at_rank]][ 1:25 ]
}

# Order the microbial features according to the microbiome_Ranking
sub_ASVs   <- reduc_ASVs[ , match(sort_microbeRank, colnames(reduc_ASVs)) ]   # reduc_ASVs[  , 1:51 ]    
sub_metabo <- MB_normal

# Keep a list of the original, ordered feature names
ori_ASVs_names   <- colnames(sub_ASVs)
ori_Metabo_names <- colnames(sub_metabo)
# In order to use parametrize the fit.formula in LME or correlation, we must 
# adjust the names of the features (columns). Some names start with numbers, and
# some have "-" or "," symbol: those are invalid as variable names in the formula
colnames(sub_ASVs) <- sapply( colnames(sub_ASVs) , function(x)  gsub("\\[", "", x) )
colnames(sub_ASVs) <- sapply( colnames(sub_ASVs) , function(x)  gsub("\\]", "", x) )
colnames(sub_ASVs) <- sapply( colnames(sub_ASVs) , function(x)  gsub(",", "", x) )
colnames(sub_ASVs) <- sapply( colnames(sub_ASVs) , function(x)  gsub(" ", "_", x) )
colnames(sub_ASVs) <- sapply( colnames(sub_ASVs) , function(x)  gsub("-", "", x) )
colnames(sub_metabo) <- sapply( colnames(sub_metabo) , function(x)  gsub("\\[", "", x) )
colnames(sub_metabo) <- sapply( colnames(sub_metabo) , function(x)  gsub("\\]", "", x) )
colnames(sub_metabo) <- sapply( colnames(sub_metabo) , function(x)  gsub(",", "", x) )
colnames(sub_metabo) <- sapply( colnames(sub_metabo) , function(x)  gsub(" ", "_", x) )
colnames(sub_metabo) <- sapply( colnames(sub_metabo) , function(x)  gsub("-", "", x) )
colnames(sub_metabo) <- paste0("m_", colnames(sub_metabo)) 

# Combine in one data.frame all the features: microbes, metabolites, response var (as.numeric)
combDF   <- cbind( sub_ASVs, sub_metabo, sub_MetaData )   

# Create an empty dataframe to store correlation results
corr_long <- data.frame( matrix(ncol=5, nrow = 0))
corr_colnames <- c( "Feat_metabo", "Feat_uBiome", "p_value", "Var", "Corr_Est"  )
# Loop through all x (metabolites) and y (microbe) features, one pair at a time, 
# calculate the Partial Spearman Correlation, adjust for grouping variables,
# such as Age, Sex and (random effect) individual ID.
for( xx in  seq(1,dim(sub_metabo)[2]) ){
  for( yy in  seq(1,dim(sub_ASVs)[2]) ){
    mod_name_metabo <- colnames(sub_metabo)[xx]
    mod_name_asv    <- colnames(sub_ASVs)[yy]
    # Create the formula for fixed effects
    cor.formula <- as.formula(paste0(mod_name_asv," | ",mod_name_metabo," ~ Age_Rank + Sex + Exp_Name"))

    # Perform the partial Spearman correlation
    pcor.fit <- PResiduals::partial_Spearman( cor.formula, data=combDF, # fit.x="poisson", fit.y="orm",
                                              fisher=TRUE,  conf.int=0.95)

    # Store results from the correlation model; long format: one row for each xx-yy pair of features
    corr_long <- rbind( corr_long, c( ori_Metabo_names[xx], ori_ASVs_names[yy], 
                                  pcor.fit$TS$TB$pval, pcor.fit$TS$TB$var, pcor.fit$TS$TB$ts))
  }
}
# Adjust the columns types and names 
corr_long[, 3:5] <- sapply(corr_long[, 3:5], function(x) as.numeric(x))
colnames(corr_long) <- corr_colnames

# Factorize the X and Y feature's names, so that the heatplot rows and cols are fixed consistently with other scripts
corr_long$Feat_uBiome <- factor( corr_long$Feat_uBiome, levels=rev(sort_microbeRank) )
corr_long$Feat_metabo <- factor( corr_long$Feat_metabo, levels=sort_metaboRank )



# --- Adjust p_value -----------------------------------------------------------
# - by Bonferroni adjustment
# corr_long$adjusted_p_value <- corr_long$p_value * (dim(sub_ASVs)[2]*dim(sub_metabo)[2])

# - by Benjamin-Hochber adjustment
# corr_long$adjusted_p_value_ORI <- p.adjust(corr_long$p_value, method = "BH")

# - my manual Benjamin-Hochber adjustment
m <- length(corr_long$p_value)          # Total number of tests
ranked_p <- rank(corr_long$p_value)     # Ranks of the p-values
sorted_p <- sort(corr_long$p_value)     # Sorted p-values
adjusted_p <- rep(NA, m)                # Initialize to store adjusted p-values
# BH adjustment
for (jj in m:1) {
  adjusted_p[jj]   <- min( 1 , sorted_p[jj] * m/jj)
  if (jj < m) { 
    adjusted_p[jj] <- min( adjusted_p[jj] , adjusted_p[jj+1])
  }
}
corr_long$adjusted_p_value <- adjusted_p[ranked_p]



# --- Prepare wide corr and p-val tables ---------------------------------------

# Transform the results into a wide format, needed to plot results.
df_wide <- corr_long %>% 
            select( "Feat_metabo", "Feat_uBiome", "Corr_Est" ) %>%
            pivot_wider( names_from  = Feat_metabo,
                         values_from = Corr_Est  )
# Adjust the data to a matrix format, rownames and columns order
df_wide$Feat_uBiome <- gsub("\\[", "", df_wide$Feat_uBiome)
df_wide$Feat_uBiome <- gsub("\\]", "", df_wide$Feat_uBiome) 
corr_wide <- as.data.frame(df_wide)
rownames(corr_wide) <- corr_wide$Feat_uBiome
corr_wide <- corr_wide[ ,-1]
corr_wide <- corr_wide[ , match(sort_metaboRank, colnames(corr_wide)) ] 


df_wide <- corr_long %>% 
            select( "Feat_metabo", "Feat_uBiome", "p_value" ) %>%
            pivot_wider( names_from  = Feat_metabo,
                         values_from = p_value  )
# Adjust the data to a matrix format, rownames and columns order
wide_pNom <- as.data.frame(df_wide)
rownames(wide_pNom) <- wide_pNom$Feat_uBiome
wide_pNom <- wide_pNom[ ,-1]
wide_pNom <- wide_pNom[ , match(sort_metaboRank, colnames(wide_pNom)) ] 

df_wide <- corr_long %>% 
            select( "Feat_metabo", "Feat_uBiome", "adjusted_p_value" ) %>%
            pivot_wider( names_from  = Feat_metabo,
                         values_from = adjusted_p_value  )
# Adjust the data to a matrix format, rownames and columns order
wide_pAdj <- as.data.frame(df_wide)
rownames(wide_pAdj) <- wide_pAdj$Feat_uBiome
wide_pAdj <- wide_pAdj[ ,-1]
wide_pAdj <- wide_pAdj[ , match(sort_metaboRank, colnames(wide_pAdj)) ] 

most_significant_features <- as.data.frame( sort((rowSums (wide_pAdj<0.05 ) / ncol(wide_pAdj)) *100) )
most_significant_metabo   <- as.data.frame( sort((colSums (wide_pAdj<0.05 ) / nrow(wide_pAdj)) *100) )

# --- Plot the heatmap  -----------------------------------------------------------

# Create and save the correlation and adjusted p-value as two separate images, 
# The correlation plot is created to show all correlation with nominal p-value 
# < 0.05, whereas the symbolic plot is mostly transparent and contains only the 
# symbol "*" to indicate the correlation whose p-value is significant after 
# adjustment (via Benjamin-Hochberg)
# We will then use correlation plot as first layer, and the symbolic plot is the
# second layer (mostly transparent except the symbols) that can be overlaid on 
# top of the first
# +++++  see NOTE 1 and NOTE 2 in Description at the beginning of script  ++++++

# First, plot and save the Correlation plot
if(aggreg_at_rank == "Genus") {
  tiff_corr <- paste0(path_saveFig, "Figure_3_A_Layer1.png")
}else if (aggreg_at_rank == "Family"){
  tiff_corr <- paste0(path_saveFig, "Suppl_Figure_6_A_Layer1.png")
}
png( tiff_corr, res=200, pointsize=7, width=2800, height=2100)
corrplot::corrplot( as.matrix(corr_wide),
                    p.mat=as.matrix(wide_pNom),
                    sig.level=0.05,
                    is.corr=F,  # to plot a matrix with corrplot func(), but no operations
                    addgrid.col='#BBBBBB',
                    insig='blank',  method='color',
                    cl.pos='n',
                    tl.cex = 1.75,             # size and color of the text labels
                    tl.col='#333333')
dev.off()

# Second, plot and save the adjusted p-value sign
# Make it all transparent, but keep the text labels, to maintain proportion 
# that match and overlap exactly with the correlation plot
if(aggreg_at_rank == "Genus") {
  tiff_pval <- paste0(path_saveFig, "Figure_3_A_Layer2.png")
}else if (aggreg_at_rank == "Family"){
  tiff_pval <- paste0(path_saveFig, "Suppl_Figure_6_A_Layer2.png")
}
png( tiff_pval, res=200, pointsize=7, width=2800, height=2100, bg="transparent")
corrplot::corrplot( as.matrix(corr_wide),   
                    p.mat=as.matrix( wide_pAdj ), 
                    sig.level=0.05,
                    is.corr=F,    # to plot a matrix with corrplot func(), but no operations
                    insig='label_sig' ,  method='square', 
                    col=NA, bg=NA, addgrid.col=NA,       # NA color, i.e., they are transparent
                    cl.pos='n', 
                    pch.cex=3,  pch.col="#EEEEEE",    # options for mark "*"  "#BDD800"
                    addshade = "negative", 
                    tl.cex = 1.75,             # size and color of the text labels
                    tl.col='#333333' )
dev.off()


# Combined the two layers and maintain the transparencies
library(magick)
image_corr <- magick::image_read( tiff_corr )
image_pval <- magick::image_read( tiff_pval )
combined_layers <- magick::image_composite( image_corr, image_pval, operator="over")
plot(combined_layers)
image_write(combined_layers, path=paste0(path_saveFig, "Figure_3_A_Combined.png") )
# Save location
if(aggreg_at_rank == "Genus") {
  image_write(combined_layers, path=paste0(path_saveFig, "Figure_3_A_PartialSpearman_Genus.tiff") )
}else if (aggreg_at_rank == "Family"){
  image_write(combined_layers, path=paste0(path_saveFig, "Suppl_Figure_6_PartialSpearman_Family.tiff") )
}



