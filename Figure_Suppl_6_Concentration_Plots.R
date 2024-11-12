# +++++ DESCRIPTION - Figure_Suppl_6_Concentration_Plots.R +++++++++++++++++++++++++++
# 
# Create Supplementary Figure number 6, specifically the concentration plots 
# subfigures; line plots that follow the evolution of the relative abundance of 
# selected metabolite(s) and microbe(s) over time to observe how the two changes
# in respect to each others.
# - SF 6 B - One Class of metabolites vs one microbe       --> hnd_E
# - SF 6 C - One metabolite vs one microbe                 --> hnd_C
# - SF 6 D - One Class of metabolites vs one microbe       --> hnd_D
# - SF 6 E - One metabolite vs a group of microbes         --> hnd_E
# 
# NOTE: sub-figure A can be created with script "Figure_3_Concentration_Plots.R"
#       by simply running it for the taxonomic rank "Family".
#       So, we avoid a ful copy-paste of that code here.
# 
# Essentially, this scripts is devide from script:
# - Concentration_Plot.R
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- Load Libraries -----------------------------------------------------------
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(dplyr)
library(ggrepel)
library(ggtext)
library(showtext)


# Define function to perform division by zero and return zero (not NaN)
"%/0%" <- function(x,y) { res <- x/y;   res[is.na(res)] <-0;  return(res) }

# General plotting setting
figFontSize <- 15
theme_Plot <- theme(
  plot.margin = margin(t=0.3, r=0.5, b=0.5, l=0.3, "cm"),
  panel.background = element_blank(),
  panel.grid  = element_blank(),
  plot.title  = element_text( color="#222222", size=figFontSize-1, face = "bold" ) ,
  axis.text   = element_blank(),
  axis.text.x = element_text( color="#444444", size=figFontSize-2),
  axis.text.y = element_text( color="#444444", size=figFontSize-2),
  axis.ticks  = element_blank(),
  legend.position = "none"
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
sele_cols <- c("Exp_TimePoint", "Exp_Name", "Sex", "BMI_cat", "Age_Rank")
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

# Upload refence table of ranking for metabolites and microbes
metabo_Ranking <- read_tsv( paste0(path_2_data, "Data_NMR/Metabolites_Ranking.txt"),
                            col_names = TRUE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                            guess_max  = min( 10 , 300 ) )
metabo_Ranking <- as.data.frame(metabo_Ranking)
metabo_Ranking <- metabo_Ranking %>% dplyr::select( c("Feature_Name","F_Type","Sub_Color") )
colnames(metabo_Ranking) <- c("Feature","F_Type","Colour")



# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# --- B - One Class of metabolites vs one microbe taxa --------------------------------------------
sele_taxa   <- "Ruminococcaceae"
sele_F_Type <- "SCFA"                 # unique(metabo_Ranking$F_Type)
sele_metabo <- metabo_Ranking[ metabo_Ranking$F_Type==sele_F_Type , "Feature" ]
sele_metabo <- sele_metabo[ sele_metabo %in% colnames(MB_normal) ]
basic_palette <- metabo_Ranking[ metabo_Ranking$Feature %in% sele_metabo , "Colour" ]
basic_palette <- c( "#2593BD",basic_palette )
sub_MetaData <- MetaData %>% dplyr::select( all_of(  c("Exp_TimePoint", "Exp_Name") ) )

scale_B <- 0.5
comb_df <- data.frame( sub_MetaData, "Microbe"=reduc_ASVs[[sele_taxa]]*100 /scale_B, MB_normal[sele_metabo]*100 )
colnames(comb_df)[colnames(comb_df)=="Microbe"] <- sele_taxa
comb_df <- comb_df %>% pivot_longer( cols     = -c(Exp_Name, Exp_TimePoint) ,  # exclude specific column
                                     names_to  = "Name_Feature",               # selected column for a new variable ("groupings")
                                     values_to = "Values")                     # values all into one new variables
comb_df$Name_Feature <- factor( comb_df$Name_Feature, levels=c(sele_taxa , sele_metabo))

dsPlot <- comb_df[comb_df$Exp_Name=="A04",] 
ymax <- 55 # ceiling( max(dsPlot$Values) )   

hnd_B <- ggplot( dsPlot, aes(x=Exp_TimePoint, y=Values, group=Name_Feature)) +
        ## Geometric annotations that act as grid lines
        geom_segment(
          data = tibble(x=c(0,7,14,28), y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey91", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(
          data = tibble(y=seq(0, ymax, by=5), x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
          color = "grey91", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # X-axis line
          data = tibble(y=0, x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
          color = "grey10", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # left Y-axis line
          data = tibble(x=0, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey10", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # right Y-axis line
          data = tibble(x=41, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey0", size = 1 , inherit.aes = FALSE
        ) +
        ## Add the labels first
        geom_text_repel(
          data = dsPlot %>% filter( Exp_TimePoint == "28"),
          aes( label = Name_Feature), color = "grey40",
          size = 3.5, direction = "y", xlim = c(32, 40),
          hjust = 0, segment.size=.4, segment.alpha=.6, segment.linetype = "dotted",
          box.padding = .4, segment.curvature = 0.5, segment.ncp = 0.1, segment.angle = 20, 
        ) +
        ## Lines of the Features   
        geom_line(  aes(color=Name_Feature), linewidth=1) +
        geom_point( aes(color=Name_Feature), size=2 ) +
        scale_x_continuous( expand = c(0,0), limits = c(-0.5, 41), breaks = c(0,7,14,28) ) +
        scale_y_continuous(
              name = paste0("Metabolite adundance (%)") ,   
              sec.axis = sec_axis( ~.*scale_B, name= paste0(sele_taxa," adundance (%)")) ) +
        ggtitle( paste0("SCFAs"," vs ",sele_taxa) ) + ylab("") +
        xlab("Time (days)") +
        scale_color_manual( values=basic_palette ) + 
        theme_Plot +
        theme(
          axis.title.y.left  = element_text(color = "#B244BE"),
          axis.title.y.right = element_text(color ="#2593BD")
        )
hnd_B



# # --- C - One metabolite vs one microbe --------------------------------------------
sele_taxa   <- "Ruminococcaceae"
sele_F_Type <- "Carboxylic acids"                 # unique(metabo_Ranking$F_Type)
sele_metabo <- metabo_Ranking[ metabo_Ranking$F_Type==sele_F_Type , "Feature" ]
sele_metabo <- sele_metabo[ sele_metabo %in% colnames(MB_normal) ]
basic_palette <- metabo_Ranking[ metabo_Ranking$Feature %in% sele_metabo , "Colour" ]
basic_palette <- c( "#2593BD",basic_palette )
sub_MetaData <- MetaData %>% dplyr::select( all_of(  c("Exp_TimePoint", "Exp_Name") ) )

scale_D <- 33
comb_df <- data.frame( sub_MetaData, "Microbe"=reduc_ASVs[[sele_taxa]]*100 /scale_D, MB_normal[sele_metabo]*100 )
colnames(comb_df)[colnames(comb_df)=="Microbe"] <- sele_taxa
comb_df <- comb_df %>% pivot_longer( cols     = -c(Exp_Name, Exp_TimePoint) ,  # exclude specific column
                                     names_to  = "Name_Feature",               # selected column for a new variable ("groupings")
                                     values_to = "Values")                     # values all into one new variables
comb_df$Name_Feature <- factor( comb_df$Name_Feature, levels=c(sele_taxa , sele_metabo))

dsPlot <- comb_df[comb_df$Exp_Name=="A05",] 
dsPlot <- dsPlot[ !is.na(dsPlot$Name_Feature), ]
ymax <- 1.1 #ceiling( max(dsPlot$Values) )               # ceiling( max(dsPlot$Values) )

hnd_C <- ggplot( dsPlot, aes(x=Exp_TimePoint, y=Values, group=Name_Feature)) +
          ## Geometric annotations that act as grid lines
          geom_segment(
            data = tibble(x=c(0,7,14,28), y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
            color = "grey91", size = .5 , inherit.aes = FALSE
          ) +
          geom_segment(
            data = tibble(y=seq(0, ymax, by=0.1), x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
            color = "grey91", size = .5 , inherit.aes = FALSE
          ) +
          geom_segment(     # X-axis line
            data = tibble(y=0, x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
            color = "grey10", size = .5 , inherit.aes = FALSE
          ) +
          geom_segment(     # left Y-axis line
            data = tibble(x=0, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
            color = "grey10", size = .5 , inherit.aes = FALSE
          ) +
          geom_segment(     # right Y-axis line
            data = tibble(x=41, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
            color = "grey0", size = 1 , inherit.aes = FALSE
          ) +
          ## Add the labels first
          geom_text_repel(
            data = dsPlot %>% filter( Exp_TimePoint == "28"),
            aes( label = Name_Feature), color = "grey40",
            size = 3.5, direction = "y", xlim = c(32, 40),
            hjust = 0, segment.size=.4, segment.alpha=.6, segment.linetype = "dotted",
            box.padding = .4, segment.curvature = 0.5, segment.ncp = 0.1, segment.angle = 20, 
          ) +
          ## Lines of the Features   
          geom_line(  aes(color=Name_Feature), linewidth=1) +
          geom_point( aes(color=Name_Feature), size=2 ) +
          scale_x_continuous( expand = c(0,0), limits = c(-0.5, 41), breaks = c(0,7,14,28) ) +
          scale_y_continuous(
            name = paste0("Metabolite adundance (%)") ,   
            sec.axis = sec_axis( ~.*scale_D, name= paste0(sele_taxa," adundance (%)")) ) +
          ggtitle( paste0("Carboxylic acids"," vs ",sele_taxa) ) + ylab("") +
          xlab("Time (days)") +
          scale_color_manual( values=basic_palette ) + 
          theme_Plot +
          theme(
            axis.title.y.left  = element_text(color = "#800028"),
            axis.title.y.right = element_text(color ="#2593BD")
          )
hnd_C



# --- D - One Class of metabolites vs one microbe --------------------------------------------
sele_taxa   <- "Bacteroidaceae"
sele_F_Type <- "Amino acids (essential)"                 # unique(metabo_Ranking$F_Type)
sele_metabo <- metabo_Ranking[ metabo_Ranking$F_Type==sele_F_Type , "Feature" ]
sele_metabo <- sele_metabo[ sele_metabo %in% colnames(MB_normal) ]
basic_palette <- metabo_Ranking[ metabo_Ranking$Feature %in% sele_metabo , "Colour" ]
basic_palette <- c( "#2593BD",basic_palette )
sub_MetaData <- MetaData %>% dplyr::select( all_of(  c("Exp_TimePoint", "Exp_Name") ) )

scale_D <- 3.5
comb_df <- data.frame( sub_MetaData, "Microbe"=reduc_ASVs[[sele_taxa]]*100 /scale_D, MB_normal[sele_metabo]*100 )
colnames(comb_df)[colnames(comb_df)=="Microbe"] <- sele_taxa
comb_df <- comb_df %>% pivot_longer( cols     = -c(Exp_Name, Exp_TimePoint) , # exclude specific column
                                     names_to  = "Name_Feature",              # selected column for a new variable ("groupings")
                                     values_to = "Values")                    # values all into one new variables
comb_df$Name_Feature <- factor( comb_df$Name_Feature, levels=c(sele_taxa , sele_metabo))

dsPlot <- comb_df[comb_df$Exp_Name=="A04",] 
ymax <- 3 #ceiling( max(dsPlot$Values) )               # ceiling( max(dsPlot$Values) )
  
hnd_D <- ggplot( dsPlot, aes(x=Exp_TimePoint, y=Values, group=Name_Feature)) +
        ## Geometric annotations that act as grid lines
        geom_segment(
          data = tibble(x=c(0,7,14,28), y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey91", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(
          data = tibble(y=seq(0, ymax, by=1), x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
          color = "grey91", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # X-axis line
          data = tibble(y=0, x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
          color = "grey10", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # left Y-axis line
          data = tibble(x=0, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey10", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # right Y-axis line
          data = tibble(x=41, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey0", size = 1 , inherit.aes = FALSE
        ) +
        ## Add the labels first
        geom_text_repel(
          data = dsPlot %>% filter( Exp_TimePoint == "28"),
          aes( label = Name_Feature), color = "grey40",
          size = 3.5, direction = "y", xlim = c(32, 40),
          hjust = 0, segment.size=.4, segment.alpha=.6, segment.linetype = "dotted",
          box.padding = .4, segment.curvature = 0.5, segment.ncp = 0.1, segment.angle = 20, 
        ) +
        ## Lines of the Features   
        geom_line(  aes(color=Name_Feature), linewidth=1) +
        geom_point( aes(color=Name_Feature), size=2 ) +
        scale_x_continuous( expand = c(0,0), limits = c(-0.5, 41), breaks = c(0,7,14,28) ) +
        scale_y_continuous(
              name = paste0("Metabolite adundance (%)") ,   
              sec.axis = sec_axis( ~.*scale_D, name= paste0(sele_taxa," adundance (%)")) ) +
        ggtitle( paste0("Amino acids"," vs ",sele_taxa) ) + ylab("") +
        xlab("Time (days)") +
        scale_color_manual( values=basic_palette ) + 
        theme_Plot +
        theme(
          axis.title.y.left  = element_text(color = "#71CB17"),
          axis.title.y.right = element_text(color = "#2593BD")
        )
hnd_D
 


# - E - One metabolite vs a group of microbes ----------------------------------
# sele_taxa   <- colnames(reduc_ASVs)[1:10] 
sele_taxa   <- c("Lachnospiraceae", "Ruminococcaceae","Bacteroidaceae","Bifidobacteriaceae","Christensenellaceae","Rikenellaceae","Eggerthellaceae")
sele_metabo <- "Butyrate"
sub_MetaData <- MetaData %>% dplyr::select( all_of(  c("Exp_TimePoint", "Exp_Name") ) )

scale_E <- 0.5
comb_df <- data.frame( sub_MetaData, reduc_ASVs[sele_taxa]*100, MB_normal[sele_metabo]*100 /scale_E )
colnames(comb_df)[colnames(comb_df)=="Microbe"] <- sele_taxa
comb_df <- comb_df %>% pivot_longer( cols     = -c(Exp_Name, Exp_TimePoint) ,  # exclude specific column
                                     names_to  = "Name_Feature",               # selected column for a new variable ("groupings")
                                     values_to = "Values")                     # values all into one new variables

# Import color LUT for microbiome ranks and set order to match color palette
upper_rank <- "Phyla"
ref_ranks  <- read_tsv(paste0(path_2_data, "Data_Seq/RefList_",upper_rank,".txt"),
                       col_names = FALSE, skip_empty_rows = TRUE ,show_col_types = FALSE,
                       guess_max  = min( 10 , 300 ) )

# The columns names from aggreg_at_rank and upper_rank, those that appear in comb_df dataset
mask      <- reduc_TAXO[[aggreg_at_rank]] %in% sele_taxa
sub_TAXO  <- reduc_TAXO[ mask , c(upper_rank, aggreg_at_rank)]
ref_ranks <- ref_ranks[ ,-c(2)]
colnames(ref_ranks) <- c(upper_rank, "Colour")
# First, merge reduc_TAXO with ref_rank to add the "Color" column (matching upper_rank ID) to reduc_TAXO
sub_TAXO <- merge( sub_TAXO, ref_ranks, by=upper_rank)
# Second, merge with the comb_df, to add Colors and upper_rank to final plotting dataset
temp <- unique( sub_TAXO[, c(aggreg_at_rank,"Colour",upper_rank)])
colnames(temp)[1] <- "Name_Feature"
comb_df <- merge( temp, comb_df, by="Name_Feature", all=TRUE)
comb_df[ comb_df$Name_Feature == sele_metabo, "Colour" ] <- "#323232"

comb_df <- comb_df[ !is.na(comb_df$Name_Feature), ]

comb_df$Name_Feature <- factor( comb_df$Name_Feature, levels=c(sele_taxa , sele_metabo))
# Finally assign the colors to the palette
temp <- distinct(comb_df, Name_Feature, .keep_all=T)
basic_palette <- temp[ order(temp$Name_Feature), "Colour" ]

dsPlot <- comb_df[comb_df$Exp_Name=="A10",]
ymax <- 40                     # ceiling( max(dsPlot$Values, na.rm=T) )   

hnd_E <- ggplot( dsPlot, aes(x=Exp_TimePoint, y=Values, group=Name_Feature)) +
        ## Geometric annotations that act as grid lines
        geom_segment(
          data = tibble(x=c(0,7,14,28), y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey91", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(
          data = tibble(y=seq(0, ymax, by=5), x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
          color = "grey91", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # X-axis line
          data = tibble(y=0, x1=0, x2=28),  aes(x=x1, xend=x2, y=y, yend=y),
          color = "grey10", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # left Y-axis line
          data = tibble(x=0, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey10", size = .5 , inherit.aes = FALSE
        ) +
        geom_segment(     # right Y-axis line
          data = tibble(x=46, y1=0, y2=ymax),  aes(x=x, xend=x, y=y1, yend=y2),
          color = "grey0", size = 1 , inherit.aes = FALSE
        ) +
        ## Add the labels first
        geom_text_repel(
          data = dsPlot %>% filter( Exp_TimePoint == "28"),
          aes( label = Name_Feature), color = "grey40", max.overlaps=14,
          size = 2.5, direction = "y", xlim = c(32, 46),
          hjust = 0, segment.size=.4, segment.alpha=.6, segment.linetype = "dotted",
          box.padding = .4, segment.curvature = 0.5, segment.ncp = 0.1, segment.angle = 20,
        ) +
        ## Lines of the Features
        geom_line(  aes(color=Name_Feature), linewidth=1) +
        geom_point( aes(color=Name_Feature), size=2 ) +
        scale_x_continuous( expand = c(0,0), limits = c(-0.5, 46), breaks = c(0,7,14,28) ) +
        scale_y_continuous(
          name =  paste0("microbe Family adundance (%)") ,   
          sec.axis = sec_axis( ~.*scale_E, name=paste0("Metabolite adundance (%)")) ) +
        ggtitle( paste0(sele_metabo," vs microbe families") ) + ylab("") +
        xlab("Time (days)") + 
        theme_Plot +
        scale_color_manual( values=basic_palette )  +
        theme(
          axis.title.y.left  = element_text(color = "#005E76"),
          axis.title.y.right = element_text(color = "#323232")
        )
hnd_E


# --- COMBINE sub-Figures ------------------------------------------------------
# Place sub-plots together into one figure
library("cowplot")

ggdraw()
ggdraw() +
  draw_plot(hnd_D, x = .01, y = .52, width =.44, height =.44) +
  draw_plot(hnd_C, x = .53, y = .52, width =.44, height =.44) +
  draw_plot(hnd_B, x = .01, y = .01, width =.44, height =.44) +
  draw_plot(hnd_E, x = .53, y = .01, width =.44, height =.44) +
  draw_plot_label(label = c("B", "C", "D", "E"), size = 17, x = c(.01, .53, .01, .53), y = c(1,1,.49,.49) )

ggsave( paste0(path_saveFig, "Suppl_Figure_6_BCDE.tiff"), width=10, height=7.5, dpi=300)   # height= 4 or 7.5


















