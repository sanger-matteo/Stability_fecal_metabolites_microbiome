# +++++ DESCRIPTION - Figure_1_Composition.R +++++++++++++++++++++++++++++++++++
#
# This script creates the results for the first figure of the Longitudinal paper
# There are four sub-figures:
# - Fig 1 A - composition stacked bar-plot of microbiome features     --> hnd_A
# - Fig 1 B - composition stacked bar-plot of metabolites             --> hnd_B
# - Fig 1 C - t-SNE clustering using microbiome data                  --> hnd_C
# - Fig 1 D - t-SNE clustering using measured metabolites             --> hnd_D
#
# Lastly, at section Combine all the "handles" of the 4 subplots (hnd_X) are
# arranged into one publication-ready figure. 
# NOTE: legends are not plotted.
#
# The taxa and metabolites are colored consistently across the paper, using two 
# file that rank and assign specific color tone and saturation, so that features 
# from the same class share the same tone and the plots appear better structured 
# and easier to understand
#
# Essentially it combines 4 previous scripts into one:
# - microbiome_Compositional_Plot.R
# - metabo_Compositional_Plot.R
# - microbiome_tSNE.R
# - metabo_tSNE.R
# It is not streamlined or efficient, because each section re-loads data and 
# repeats many operations, but this ensure consistency with original scripts 
# and it is easier to implement and change each section without affecting the 
# others.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Date:     2023
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




# --- Load ALL Libraries -----------------------------------------------------------
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(Rtsne)

cDir_Rscript <- setwd("/Users/sm/Documents/R/Longitudinal_Analysis")
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures_Microbiome/" )

figFontSize <- 12

# Import custom functions, specifically contains those to calculated microbiome diversity metrics
source( paste0(cDir_Rscript, "/SangerMicrobiome.R") )



# --- Fig 1 A ------------------------------------------------------------------
# Load ASV table(s) and Metadata
ASV  <- read_tsv( paste0( path_2_data,"Data_Seq/ASV_longit.tsv"), 
                  skip_empty_rows=TRUE, show_col_types=FALSE,
                  guess_max=min( 10 , 300 ) )
TAXO <- read_tsv( paste0( path_2_data,"TData_Seq/AXO_longit.tsv"), 
                  skip_empty_rows=TRUE, show_col_types=FALSE,
                  guess_max=min( 10 , 300 ) ) 
MetaData <- read_tsv( paste0(path_2_data, "Data_Seq/META_longit.txt"), 
                      skip_empty_rows=TRUE, show_col_types=FALSE,
                      guess_max=min( 10 , 300 ) ) 

# We must remove samples with zero reads
MetaData <- filter(MetaData, Sample_TotReads >= 500)
ASV  <- ASV[ ASV$Sample_ID %in% MetaData$Sample_ID , ]

# --- Pruning MICROBIOME data
# Convert from count feature table to proportion (normalized)
ASV <- as.data.frame( dplyr::select( ASV, -c(Sample_ID, Sample_TotReads) ) )
ASV <- ASV/rowSums(ASV)
# Filtered taxonomic features with a Rel. abundance less than (0.03%) in greater than 10% of all samples
thres_abundance  <- 0.0003
thres_freqSample <- 0.010
mask <- unlist( colSums(ASV > thres_abundance) / dim(ASV)[1] >= thres_freqSample, use.names=FALSE )
ASV  <- ASV[,mask]
TAXO <- TAXO[,c(T,mask)]

# --- Prepare data 
# Transpose both ASV and TAXO, placing features as rows. Then they can be joined
# (cbind) as one DF and input to function groupBy_Ranks. Also, strip taxonomic 
# names of prefix (D__) and remove unnecessary columns (Sample_ID, RowName)
listDF <- combineReady_ASV_TAXO( ASV, TAXO, MetaData$Sample_ID ) 
ASV    <- listDF[[1]]
TAXO   <- listDF[[2]]

# -- Group data
# Group by 2 hierarchical taxonomic ranks
UpperRank_1 <- "Class"
LowerRank_1 <- "Genus"
thres_UP <- 0.01
thres_LW <- 0.03
listDF <- groupBy_2_Ranks( ASV, TAXO , UpperRank_1, LowerRank_1, thres_UP, thres_LW )
g_ASV  <- listDF[[1]]
g_TAXO <- listDF[[2]]
dsPlot <- listDF[[3]]

# [OPTION] Add samples' information, that can help organize data along x-axis 
# and display in the plot as groups of samples sharing a common aspect
dsPlot <- dsPlot %>% mutate( ., Exp_Name="", .after=Sample_ID)
LoopList <- unique(dsPlot$Sample_ID)
for (kk in seq(1,length(LoopList))){
  mask <- dsPlot$Sample_ID==LoopList[kk]
  dsPlot[ mask, "Exp_Name"] <- MetaData[ MetaData$Sample_ID==LoopList[kk] , "Exp_Name" ][[1]]
}
# Replace two ID tags, to simplify and make the name sequential series 
dsPlot$Exp_Name <- sub( "A15", "A13", dsPlot$Exp_Name)
dsPlot$Exp_Name <- sub( "A20", "A14", dsPlot$Exp_Name)

dsPlot <- dsPlot %>% mutate( ., Exp_TimePoint="", .after=Sample_ID)
LoopList <- unique(dsPlot$Sample_ID)
for (kk in seq(1,length(LoopList))){
  mask <- dsPlot$Sample_ID==LoopList[kk]
  dsPlot[ mask, "Exp_TimePoint"] <-  paste0( "Day ",
                                             MetaData[ MetaData$Sample_ID==LoopList[kk] , 
                                             "Exp_TimePoint" ][[1]] )
}

# ----- Plot compositional data
# Organize and order the x-axis samples         x=factor( Sample_ID, levels=xAxOrder$Sample_ID )
xAxOrder <- MetaData %>% dplyr::select( c(Sample_ID, Exp_Name, Exp_TimePoint))
xAxOrder <- xAxOrder[ order(xAxOrder$Exp_Name, xAxOrder$Exp_TimePoint), ]

# Define color scale
#---[ OPTION 2 ] 
gradColors_1 <- ColourMultiGradient( dsPlot, UpperRank_1, LowerRank_1, LightenAllOthers=TRUE )

dsPlot_1 <- dsPlot
# Stacked bar-plot
hnd_A <- ggplot( dsPlot_1, aes( x=Sample_ID,
                            y=Counts, 
                            fill=factor( .data[[LowerRank_1]], levels=unique(dsPlot_1[ ,LowerRank_1]))) ) + 
        geom_col( position="fill", colour="grey", size=0.1) +         
        scale_fill_manual( "", values=gradColors_1) +
        scale_y_continuous( labels = scales::percent_format(),
                            expand = c(0, NA),limits = c(0, NA)) +    # Scale y axis into % 
        
        # [OPTION] Group x axis by a variable and plot in different grids blocks
        facet_grid( .~Exp_Name, switch="y", scale="free", space="free_x" ) +
        scale_x_discrete( labels= rep( c("d-0", "d-7", "d-14", "d-28"), 13) ) +
        
        xlab("") + ylab("") +  
        guides(fill=guide_legend(ncol=1)) +   # display legend keys as N columns
        theme(
          plot.margin = margin(t=0.2, r=0.2, b=0.2, l=0.2, "cm"),
          panel.background = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.grid.major = element_blank() , 
          plot.title   = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
          axis.title.x = element_text( color="#444444", size=figFontSize  , vjust = -.6, hjust = .5) ,
          axis.title.y = element_text( color="#444444", size=figFontSize  , vjust =  .6, hjust = .5) ,
          axis.text.x  = element_text( color="#444444", size=figFontSize-4, hjust = 1, angle=60  ) , 
          axis.text.y  = element_text( color="#444444", size=figFontSize ) ,
          axis.ticks   = element_blank() ,
          axis.line.x.top    = element_line(size = .5, colour = "#4e4e4e", linetype = "solid"),
          axis.line.x.bottom = element_line(size = .5, colour = "#4e4e4e", linetype = "solid"),
          legend.position = "none",
          legend.title = element_text(size = figFontSize-2), 
          legend.text  = element_text(size = figFontSize-4),
          legend.key.size = unit(0.25, 'cm'),
          strip.text.x     = element_text( color="#444444", size=figFontSize ),
          strip.background = element_blank() 
        )
hnd_A



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Fig 1 B ------------------------------------------------------------------
# ---- Load DATA tables 
path_2_data  <- paste0( cDir_Rscript, "/" )
path_saveFig <- paste0( cDir_Rscript, "/Figures_Metabolites/" )

# Load Feature table and Metadata
MetaData <- read_tsv( paste0(path_2_data, "Data_Seq/META_longit.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

metabols <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Concentration.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

# --- Pruning METABOLITES data 
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
MB_normal <- MB_normal *100

# --- Prepare dataset for plotting
# Create a long Format table for the dataset
temp   <- MetaData %>% dplyr::select( c(Sample_ID, External_ID, Exp_Name, Exp_TimePoint) )
dsPlot <- cbind( temp, MB_normal )
dsPlot <- dsPlot %>% pivot_longer( cols= -c(Sample_ID, External_ID, Exp_Name, Exp_TimePoint) , 
                                   values_to="Counts",
                                   names_to ="Feature_Name")
# Import metabolites general classification and set order to match color palette
metabo_Ranking <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Ranking.txt"),
                            col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                            guess_max  = min( 10 , 300 ) )
metabo_Ranking <- as.data.frame(metabo_Ranking)
metabo_Ranking <- metabo_Ranking %>% dplyr::select( c("Feature_Name","F_Type","Colour") )
dsPlot <- merge( dsPlot, metabo_Ranking, by="Feature_Name", all=TRUE)

# Define color scale
UpperRank_2 <- "F_Type"
LowerRank_2 <- "Feature_Name"
list_Return <- ColourMultiGradient_Metabolites( dsPlot, UpperRank_2, LowerRank_2, 
                                                LightenAllOthers=TRUE, dynamicToneRange=TRUE)
summ_UP     <- list_Return[[2]]
summ_LW     <- list_Return[[3]]
TaxaColors  <- list_Return[[1]]
gradColors_2  <- array(unlist(TaxaColors))

# Rename ID tags for the x-axis to be simple and sequential series
dsPlot$Exp_TimePoint <- sub(  "0",  "d.00", dsPlot$Exp_TimePoint)
dsPlot$Exp_TimePoint <- sub(  "7",  "d.07", dsPlot$Exp_TimePoint)
dsPlot$Exp_TimePoint <- sub(  "14", "d.14", dsPlot$Exp_TimePoint)
dsPlot$Exp_TimePoint <- sub(  "28", "d.28", dsPlot$Exp_TimePoint)

dsPlot$Exp_Name <- sub( "A15", "A13", dsPlot$Exp_Name)
dsPlot$Exp_Name <- sub( "A20", "A14", dsPlot$Exp_Name)

dsPlot_2 <- dsPlot

# ----- Plot compositional data 
hnd_B <- ggplot( dsPlot_2, aes( x=Exp_TimePoint ,
                            y=Counts, 
                                  fill=factor( .data[[LowerRank_2]], levels=metabo_Ranking$Feature_Name ) )) + 
        geom_col( position="fill", colour="grey", size=0.025) +         
        scale_fill_manual( "", values=gradColors_2) +
        scale_y_continuous( labels = scales::percent_format(),
                            expand = c(0, NA),limits = c(0, NA)) +    # Scale y axis into % 
        # [OPTION] Group x axis by a variable and plot in different grids blocks
        facet_grid( .~Exp_Name, switch="y", scale="free", space="free_x" ) +
        scale_x_discrete( labels= rep( c("d-0", "d-7", "d-14", "d-28"), 13) ) + 
        
        xlab("") + ylab("") +  
        guides(fill=guide_legend(ncol=1)) +   # display legend keys as N columns
        theme(
          plot.margin = margin(t=0.2, r=0.2, b=0.2, l=0.2, "cm"),
          panel.background = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.grid.major = element_blank() , 
          plot.title   = element_text( color="#222222", size=figFontSize+2, face = "bold" ) ,
          axis.title.x = element_text( color="#444444", size=figFontSize  , vjust = -.6, hjust = .5) ,
          axis.title.y = element_text( color="#444444", size=figFontSize  , vjust =  .6, hjust = .5) ,
          axis.text.x  = element_text( color="#444444", size=figFontSize-4, hjust = 1, angle=60  ) , 
          axis.text.y  = element_text( color="#444444", size=figFontSize ) ,
          axis.ticks   = element_blank() ,
          axis.line.x.top    = element_line(size = .5, colour = "#4e4e4e", linetype = "solid"),
          axis.line.x.bottom = element_line(size = .5, colour = "#4e4e4e", linetype = "solid"),
          legend.position = "none",
          legend.title = element_text(size = figFontSize-2), 
          legend.text  = element_text(size = figFontSize-4),
          legend.key.size = unit(0.25, 'cm'),
          
          strip.text.x     = element_text( color="#444444", size=figFontSize ),
          strip.background = element_blank() 
        )
hnd_B



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Fig 1 C ------------------------------------------------------------------
# Set the palette of colors to use for the tSNE plots
my_palette <- c( "#299CD2", "#F6D21F", "#ff5b4f", "#11d421", "#b071e3", 
                 "#3c54de", "#f89239", "#c43d61", "#7cd9b8", "#e3b0ff", 
                 "#85d4ff", "#B79D28", "#FFB1B4", "#20911c", "#C6C9D5" )

# ---- Load DATA tables 
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

# --- Pruning MICROBIOME data 
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
# CALC GeomMean >> geom_mean <- apply( ASV, 1, function(x) exp(mean(log(x[x>0]))) )
ASV <- ASV/rowSums(ASV)
# Filtered taxonomic features with a Rel. abundance less than (0.03%) in greater than 10% of all samples
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

# Remove point A12 at d_28, because it has only about 10 reads in total !
mask <- (MetaData$Exp_Name=="A12") + (MetaData$Exp_TimePoint=="28") != 2
MetaData <- MetaData[ mask, ]
reduc_ASVs <- reduc_ASVs[ mask, ]


# --- t-SNE Plot 
# Perform t-SNE on the dataset, but then extract the data and plot using ggplot

plotMetaData <- MetaData[ , c("Exp_TimePoint", "Exp_Name", "BMI_cat", "Age_Rank", "Sex" )]     
# Revalue the tag to be ordered from A01 -> A14
plotMetaData$Exp_Name <- replace( plotMetaData$Exp_Name, plotMetaData$Exp_Name=="A15", "A13")
plotMetaData$Exp_Name <- replace( plotMetaData$Exp_Name, plotMetaData$Exp_Name=="A20", "A14")

dsPlot   <- reduc_ASVs                             
dsPlot   <- dsPlot[ , !colSums(dsPlot)==0]  # take only non-zero sum columns
colnames(dsPlot) <- seq(1, dim(dsPlot)[2])

perplex   <- 18                    # appr ~ floor((nrow(dsPlot)-1)/3)
ScaleData <- scale(dsPlot)         # Scale data
tSNE_fit  <- Rtsne( X = data.matrix( ScaleData ),  perplexity = perplex,
                    pca_center = TRUE, pca_scale  = TRUE, 
                    dims = 2 , max_iter = 1000, check_duplicates = FALSE )
tSNE_df   <- tSNE_fit$Y %>% as.data.frame()        # %>% rename(tSNE1_score="V1", tSNE2_score="V2")

# Add grouping tag accordingly: "Exp_Name", "Exp_TimePoint" 
tSNE_df$Type <- factor(plotMetaData$Exp_Name)               # add region tag to color points accordingly
tSNE_df$TimePoint <- factor(plotMetaData$Exp_TimePoint)     # add region tag to color points accordingly
fig_type <- "Exp_Name"

maxLim <- max( c( abs(floor( min(tSNE_df$V1))), abs(ceiling( max(tSNE_df$V2))), 
                  abs(floor( min(tSNE_df$V1))), abs(ceiling( max(tSNE_df$V1))) ) ) +0.5
ax_lims <- maxLim

hnd_C <- ggplot( tSNE_df, aes(x=V1, y=V2, color=Type )) +
          geom_hline(aes(yintercept = 0), size = 0.5, colour = "#cecece", linetype = "dotted") +
          geom_vline(aes(xintercept = 0), size = 0.5, colour = "#cecece", linetype = "dotted") +
          geom_point( alpha = 0.8, size=2 ) +       
          scale_colour_manual( values=my_palette )  +
          xlim( -ax_lims, ax_lims) +
          ylim( -ax_lims, ax_lims) +
          coord_fixed(ratio = 1) +
        xlab("tSNE_1") + 
        ylab("tSNE_2 ") +
        labs( color= paste0(fig_type, ":")) +
        theme(
          # panel.background = element_blank() 
          panel.background = element_rect(size = 2, colour = "#9e9e9e", linetype = "solid", fill=NA) ,
          panel.grid.minor = element_blank() ,
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_blank() , 
          plot.title   = element_text( color="#222222", size=figFontSize+4, face = "bold", hjust = .5) ,
          axis.title.x = element_text( color="#444444", size=figFontSize, vjust = -.6, hjust = .5) ,
          axis.title.y = element_text( color="#444444", size=figFontSize, vjust =  .6, hjust = .5) ,
          axis.text.x  = element_text( color="#444444", size=figFontSize-4,   vjust = -1), 
          axis.text.y  = element_text( color="#444444", size=figFontSize-4 ),
          # axis.ticks   = element_blank(),
          legend.position = "none",
          legend.key   = element_blank(),
          legend.title = element_text(size = figFontSize-2), 
          legend.text  = element_text(size = figFontSize-2)
        )
hnd_C
floor(   min(tSNE_df$V1) )
ceiling( max(tSNE_df$V1) )
floor(   min(tSNE_df$V2) )
ceiling( max(tSNE_df$V2) )



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Fig 1 D ------------------------------------------------------------------
# ---- Load DATA tables 
path_2_data  <- paste0( cDir_Rscript, "/" )

# Load Feature table and Metadata
MetaData <- read_tsv( paste0(path_2_data, "Data_Seq/META_longit.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

metabols <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Concentration.txt"),
                      skip_empty_rows = TRUE ,show_col_types = FALSE,
                      guess_max  = min( 10 , 300 ) ) 

# --- Pruning METABOLITES data 
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


# --- t-SNE Plot 
# Perform t-SNE on the dataset, but then extract the data and plot using ggplot
# The response variable is parametric (i.e. fig_type) and can be changed
plotMetaData <- MetaData[ , c("Exp_TimePoint", "Exp_Name", "BMI_cat", "Age_Rank", "Sex" )]             
# Revalue the tag to be ordered from A01 -> A14
plotMetaData$Exp_Name <- replace( plotMetaData$Exp_Name, plotMetaData$Exp_Name=="A15", "A13")
plotMetaData$Exp_Name <- replace( plotMetaData$Exp_Name, plotMetaData$Exp_Name=="A20", "A14")

dsPlot   <- MB_normal                             
dsPlot   <- dsPlot[ , !colSums(dsPlot)==0]     # take only non-zero sum columns
colnames(dsPlot) <- seq(1, dim(dsPlot)[2])

perplex   <- 14                     # appr ~ floor((nrow(dsPlot)-1)/3)
ScaleData <- scale(dsPlot)          # Scale data
tSNE_fit  <- Rtsne( X = data.matrix( ScaleData ),  perplexity = perplex,
                    pca_center = TRUE, pca_scale  = TRUE, 
                    dims = 2 , max_iter = 1000, check_duplicates = FALSE )
tSNE_df   <- tSNE_fit$Y %>% as.data.frame()     # %>% rename(tSNE1_score="V1", tSNE2_score="V2")

# Add grouping tag accordingly: "Exp_Name", "Exp_TimePoint"
tSNE_df$Type <- factor(plotMetaData$Exp_Name)               # add region tag to color points accordingly
tSNE_df$TimePoint <- factor(plotMetaData$Exp_TimePoint)     # add region tag to color points accordingly
fig_type <- "Exp_Name"

maxLim <- max( c( abs(floor( min(tSNE_df$V1))), abs(ceiling( max(tSNE_df$V2))), 
                  abs(floor( min(tSNE_df$V1))), abs(ceiling( max(tSNE_df$V1))) ) ) +0.5
ax_lims <- maxLim

hnd_D <- ggplot( tSNE_df, aes(x=V1, y=V2, color=Type )) +
        geom_hline(aes(yintercept = 0), size = 0.5, colour = "#cecece", linetype = "dotted") +
        geom_vline(aes(xintercept = 0), size = 0.5, colour = "#cecece", linetype = "dotted") +
        geom_point( alpha = 0.8, size=2 ) +
        scale_colour_manual( values=my_palette )  +
        xlim( -ax_lims, ax_lims) +
        ylim( -ax_lims, ax_lims) +
        coord_fixed(ratio = 1) +
        
        xlab("tSNE_1") + 
        ylab("tSNE_2 ") +
        labs( color= paste0(fig_type, ":")) +
        theme(
          # panel.background = element_blank() 
          panel.background = element_rect(size = 2, colour = "#9e9e9e", linetype = "solid", fill=NA) ,
          panel.grid.minor = element_blank() ,
          panel.grid.major.x = element_blank() ,
          panel.grid.major.y = element_blank() , 
          plot.title   = element_text( color="#222222", size=figFontSize+4, face = "bold", hjust = .5) ,
          axis.title.x = element_text( color="#444444", size=figFontSize, vjust = -.6, hjust = .5) ,
          axis.title.y = element_text( color="#444444", size=figFontSize, vjust =  .6, hjust = .5) ,
          axis.text.x  = element_text( color="#444444", size=figFontSize-4,   vjust = -1), 
          axis.text.y  = element_text( color="#444444", size=figFontSize-4 ),
          # axis.ticks   = element_blank(),
          legend.position = "none",
          legend.key   = element_blank(),
          legend.title = element_text(size = figFontSize-2), 
          legend.text  = element_text(size = figFontSize-2)
        )
hnd_D
floor(   min(tSNE_df$V1) )
ceiling( max(tSNE_df$V1) )
floor(   min(tSNE_df$V2) )
ceiling( max(tSNE_df$V2) )



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- COMBINE sub-Figures ------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")
ggdraw()

ggdraw() +
  draw_plot(hnd_A, x = 0.0,  y = 0.66,  width = .8, height = .33) +
  draw_plot(hnd_B, x = 0.0,  y = 0.32,  width = .8, height = .33) +
  draw_plot(hnd_C, x = 0.06, y = 0.0,   width = .33, height = .33) +
  draw_plot(hnd_D, x = 0.42, y = 0.0,   width = .33, height = .33) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 17, x = c(0, 0, 0.06, 0.42 ), y = c(1, .66, .30, .30))

ggsave( paste0(path_saveFig, "Figure_1.tiff"), width=10, height=14, dpi=300)








