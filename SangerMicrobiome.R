# +++++ SangerMicrobiome.R +++++++++++++++++++++++++++++++++++++++++++++++++++++
# 
# Set of custom function created to calculate metrics specific to microbiome 
# data, as well as manipulate and simplify dataframes with compositional data.
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Author:   Matteo Sangermani
# e-mail:   matteo.sangermani@ntnu.no
# Release:  1.0
# Release   date: 2021-2022
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



# ----- ALPHA diversity --------------------------------------------------------

met_Shannon <- function( vec ){
  # SHANNON diversity quantifies both abundance and diversity of an ecosystem.
  #     It quantifies the uncertainty in predicting the species identity of  
  #     an observation taken at random from a dataset.
  # Input:
  # - vec: a numeric vector (1-by-n) where each element is a feature
  # Output: 
  # - Shan_val: a numeric vector with calculated Shannon index
  propV    <- vec / sum(vec)
  Shan_val <-  - sum( propV*logb(propV, exp(1) ), na.rm = TRUE)
  return( Shan_val )
}


met_Simpson <- function( vec ){
  # SIMPSON is a dominance index, that assesses the probability that two 
  #     individuals randomly selected from a set will belong to the same taxa
  # Input:
  # - vec: a numeric vector (1-by-n) where each element is a feature
  # Output: 
  # - Simp_val: a numeric vector with calculated Simpson index
  propV    <- vec / sum(vec)
  Simp_val <-  1 -sum( propV^2, na.rm = TRUE)
  return( Simp_val )
}


met_Chao <- function( vec ){
  # CHAO_1 estimates the richness of a sample. It simply define singleton (F1) 
  #    and doubletons (S2) as the number of taxa with exactly 1 or 2 count.
  # Input:
  # - vec: a numeric vector (1-by-n) where each element is a feature
  # Output: 
  # - Chao_val: a numeric vector with calculated Chao value
  
  S_obs <- length( vec[vec >  0])    # number of taxa
  F_1   <- length( vec[vec == 1])    # number of taxa with exactly 1 count
  F_2   <- length( vec[vec == 2])    # number of taxa with exactly 2 count
  # Calculate Chao1 value
  if ( (F_1-F_2)^2 == (F_1+F_2)^2 ) {
    Chao_val <- S_obs +  (F_1*(F_1-1)) /(2*(F_2+1))
  }else {
    Chao_val <- S_obs +  F_1**2 /(2*F_2)
  }
  return( Chao_val )
}

# GUIDELINEs - alpha diversity metrics with threshold
# Here is how to use the q (or qq) value relates to different diversity indices:
# q = 0:    When q is 0, Hill numbers represent species richness. 
#           The resulting diversity index is equivalent to species richness 
#           (i.e., the total number of species present in the community).
# 0<q<1:    As q increases from 0 to 1, Hill numbers give increasing weight 
#           to rare species, and the resulting diversity index tends to 
#           emphasize rare species diversity more than abundant species  
#           diversity. This range is often used to capture the contribution of 
#           rare taxa  to overall diversity.
# q = 1:    When q is 1, Hill numbers represent exponential of Shannon entropy, 
#           which reflects both species richness and evenness in the community. 
#           The resulting diversity index is equivalent to Shannon diversity ,
#           index.
# q > 1:    As q increases beyond 1, Hill numbers give more weight to abundant 
#           species, and the resulting diversity index tends to emphasize 
#           evenness more than richness. Higher values of q are often used to 
#           focus on the contributions of abundant species to overall diversity.
#
# The choice of q value depends on the specific objectives of the analysis and 
# the characteristics of the microbial community being studied. In practice, 
# q values for Hill numbers are often chosen in the range of 0 to 3, with common 
# values at: 0 (richness), 1 (Shannon diversity), and 2 (Simpson diversity). 
# However, researchers may also explore other values within this range or even 
# beyond, depending on the specific research questions and the diversity  
# patterns observed in the data.
# For ACE (Abundance-based Coverage Estimator), the q value represents the
# sequencing coverage depth, which is the minimum number of reads required for 
# a species to be considered present. Typically, q values for ACE range from 1 
# to the total number of reads in the dataset, with the choice of q depending on 
# the desired level of sensitivity and accuracy in estimating species richness.

met_ACE <- function(vec, qq) {
  # ACE is and extension of Chao, estimating species richness with a variable 
  #   threshold (qq) that delimit the number of counts that define taxa as 
  #   abundant (feature_counts > qq) and rare feature_counts <= qq).
  # Input:
  # - vec: a numeric vector (1-by-n) where each element is a feature
  # Output: 
  # - ACE_val: a numeric vector with calculated ACE value
  
  vec <- vec[ vec!=0 ]                  # select non-zero feature (i.e. taxa)
  N_rare <-    sum( vec[ vec <= qq])    # total counts from all rare taxa (for qq=1 + qq=2 + ...)
  S_abun <- length( vec[ vec >  qq])    # number of abundant taxa
  S_rare <- length( vec[ vec <= qq])    # number of rare taxa
  F_1    <- length( vec[ vec == 1])     # number of taxa with exactly 1 count
  
  # Sum number of counts for rare ALL rare taxa that appears cc times
  F_kk   <- 0
  for (cc in seq(1,qq)) {
    F_kk <- F_kk + (cc* length( vec[vec==cc]))     # ALT:  (cc*(cc-1)* length( vec[vec==cc]))
  }
  # C_ace: sample coverage, estimate proportion of all individuals in rare species that are NOT singletons
  C_ace  <- 1 - (F_1)/(N_rare)
  # gamma^2 : estimate the coefficient of variation of all A_kk rare species
  gamma2 <- max( 0, (S_rare/C_ace) * (F_kk /(N_rare*(N_rare-1))) -1 )
  
  # ACE estimator of species richness:
  ACE_val <- S_abun + (S_rare/C_ace) + (F_1/C_ace) * gamma2
  if ( is.nan(ACE_val)==TRUE | ACE_val==Inf ) 
    ACE_val <- met_Chao( as.data.frame( vec ) )
  
  return( ACE_val )
}


met_HillNumber <- function( vec, qq ){
  # Hill numbers of order qq (or effective number of species) is a diversity 
  #     order. It shows the effective number of species of a community.
  #     It is unitless, but unlike Shan or Simp, has meaningful magnitude.
  # Input:
  # - vec: a numeric vector (1-by-n) where each element is a feature count
  # - qq: the order of the diversity, it defines the sensitivity of the 
  #       true diversity to rare vs. abundant species.
  #       I.E: increasing qq value increases weight given to abundant species
  # Output: 
  # - Hill_val: index calculated the vector (i.e. observation) 
  #
  # ALTernative: (good for double checks)  Hill = exp( Shannon ) 
  # but the order of diversity (qq) would be fixed at 1.0 
  
  propV    <- vec / sum(vec)
  Hill_val <- (sum(propV**qq))**(1/(1-qq))
  return( Hill_val )
}


met_RenyiEntropy <- function( vec, qq ) {
  # Renyi entropy, is simply the natural log of Hill number
  # Input:
  # - vec: a numeric vector (1-by-n) where each element is a feature count
  # - qq: the order of the diversity, weight ratio of rare vs. abundant species.
  # Output: 
  # - Renyi_val: index calculated the vector (i.e. observation)
  #
  # ALTernative: (good for double checks) RenEntr = logb( Hill_Number, exp(1))  
  
  propV   <- vec / sum(vec)
  Reny_val <- 1/(1-qq) * logb( sum(propV**qq, na.rm=TRUE) ,exp(1)) 
  return( Reny_val )
}




# ----- BETA diversity ---------------------------------------------------------

met_Jaccard <- function( vec_A, vec_B ) {
  # Jaccard index, which measure the similarity of 2 sample vectors
  # Input:
  # - vec_X: a numeric vector (1-by-n) where each element is a feature count
  # Output: 
  # - J_val: index calculated the vector (i.e. observation)
  # NOTE: it is assumed (and necessary!) that the two vector provided list a
  #       shared group of features in the exact same order.
  vec_A <- (vec_A>0) *1                 # Find non-zero features in A
  vec_B <- (vec_B>0) *1                 # Find non-zero features in B
  setAB <- ((vec_A + vec_B) >=2) *1     # Find where features appear in both DF
  J_val <- rowSums( ((vec_A + vec_B) >=2) *1 ) /
           rowSums( ((vec_A + vec_B) >=1) *1 )
  return( J_val )
}


# NOTE: essentially the ssame as met_Jaccard, but with different algorithm
met_Jaccard_Similarity_2DF <- function( df_A, df_B ) {
  # Provide DFs as either count or proportions tables.  Rows are observations 
  # and cols are features
  # Output: 
  # J_val - a one column DF, that contains the calculated Jaccard index
  #         Jaccard = intersection( A, B) / union( A, B )
  # The two DFs must have same number of rows and cols, and rows (observations)
  # must be organized in same order. The function uses the i-th row of both 
  # A and B to calculate the index for the i-th observation.
  # NOTE: Jacc_Simil-> 0 shows no similarity; Jacc_Simil-> 1 indicate similarity
  # NOTE: Jaccard_Distance = 1 - Jaccard_Similarity
  
  # Checks that the two DFs have same dimension and names rows and cols, 
  # albeit in a different order 
  if (dim(df_A)[1] != dim(df_B)[1] ) 
    stop_quietly("Error met_Jaccard_2DF(): The two tables have different row length") 
  if (dim(df_A)[2] != dim(df_B)[2] )  
    stop_quietly("Error met_Jaccard_2DF(): The two tables have different columns length") 
  
  ID_tags <- rownames(df_A)
  df_A <- df_A[,-1]
  df_B <- df_B[,-1]
  J_val <- c()
  for (ii in seq(1, length(ID_tags)) ) {
    AA <- colnames(df_A)[df_A[ii,]>=1]
    BB <- colnames(df_B)[df_B[ii,]>=1]
    intersect_AB <- length(intersect( AA, BB ))
    union_AB <- length(AA) + length(BB) - intersect_AB 
    J_val <- c( J_val, intersect_AB/union_AB) 
  }
  J_val <- data.frame( Jaccard = J_val)
  rownames(J_val) <- ID_tags
  return( J_val )
}


met_Jaccard_2DF <- function( df_A, df_B ) {
  # Provide DFs as either count or proportions tables.  Rows are observations 
  # and cols are features
  # Output: 
  # J_val - a one column DF, that contains the calculated Jaccard index
  #         Jaccard = intersection( A, B) / union( A, B )
  # The two DFs must have same number of rows and cols, and rows (observations)
  # must be organized in same order. The function uses the i-th row of both 
  # A and B to calculate the index for the i-th observation.
  # NOTE: Jacc_Simil-> 0 shows no similarity; Jacc_Simil-> 1 indicate similarity
  # NOTE: Jaccard_Distance = 1 - Jaccard_Similarity
  
  # Checks that the two DFs have same dimension and names rows and cols, 
  # albeit in a different order 
  if (dim(df_A)[1] != dim(df_B)[1] ) 
    stop_quietly("Error met_Jaccard_2DF(): The two tables have different row length") 
  if (dim(df_A)[2] != dim(df_B)[2] )  
    stop_quietly("Error met_Jaccard_2DF(): The two tables have different columns length") 
  df_A  <- (df_A[,-1]>0) *1                 # Find non-zero features in A
  df_B  <- (df_B[,-1]>0) *1                 # Find non-zero features in B
  setAB <- ((df_A + df_B) >=2) *1      # Find where features appear in both DF
  J_val <- rowSums( ((df_A + df_B) >=2) *1 ) /
    rowSums( ((df_A + df_B) >=1) *1 )
  return( data.frame( Jaccard = J_val) )
}


met_BrayCurtis <- function( vec_A, vec_B ) {
  # Bray-Curtis index, which measure the similarity of 2 sample vectors 
  # Input:
  # - vec_X: a numeric vector (1-by-n) where each element is a feature count
  # Output: 
  # - Bray-Curtis = Lesser_Values( A, B) / Sum_AllFeatures( A, B )
  #        - Lesser_Values:   sum of the "lesser values" between the set A and B
  #        - Sum_AllFeatures: sum all feature' counts present in both A and B 
  # NOTE: it is assumed (and necessary!) that the two vector provided list a
  #       shared group of features in the exact same order.
  sample_AB <- rbind(aa,bb)
  Lesser_AB <- apply( sample_AB[ c(ii,ii+1), ], 2, min )
  Lesser_AB <- sum(Lesser_AB)
  SumTot    <- sum(sample_AB)
  BC_val    <- 1 - ((2*Lesser_AB)/ SumTot) 
  return( BC_val )
}


met_BrayCurtis_2DF <- function( df_A, df_B ) {
  # Provide DF as count tables. Rows are observations and cols are features
  # Output: 
  # BC - a one column DF, that contains the calculated Jaccard index
  #      Bray-Curtis = Lesser_Values( A, B) / Sum_AllFeatures( A, B )
  #        - Lesser_Values: sum of the "lesser values" between the set A and B
  #        - Sum_AllFeatures: sum all feature' counts present in both A and B 
  # The two DFs must have same number of rows and cols, and be organized in 
  # same order. Following row order, the function uses the i-th row of both 
  # A and B to calculate the index for the i-th observation.
  
  # Checks that the two DFs have same dimension and names rows and cols, 
  # albeit in a different order 
  if (dim(df_A)[1] != dim(df_B)[1] ) 
    stop_quietly("Error met_BrayCurtis_2DF(): The two tables have different row length") 
  if (dim(df_A)[2] != dim(df_B)[2] )  
    stop_quietly("Error met_BrayCurtis_2DF(): The two tables have different columns length") 
  # if (! all( colnames(df_A) %in% colnames(df_B)) )         # Intersect cols
  #   stop_quietly("Error met_BrayCurtis_2DF(): Columns in slaveDF do not match masterDF") 
  # if (! all( rownames(df_A) %in% rownames(df_B)) )         # Intersect rows
  #   stop_quietly("Error met_BrayCurtis_2DF(): Rows in slaveDF do not match masterDF") 
  
  # The two DF have same rownames. We have to create unique row_ID by adding
  # a tag; then we Row-bind and sort by rownames
  rownames(df_A) <- paste0(rownames(df_A), "__1")
  rownames(df_B) <- paste0(rownames(df_B), "__2")
  sample_AB <- rbind( df_A, df_B)
  sample_AB <- sample_AB[ order(row.names(sample_AB)), ]
  # We can now go 2-by-2 and be sure to take rows that belong to same Sample/Patient
  BC_val <- NULL
  for (ii in seq(1,dim(sample_AB)[1], 2) ){
    Lesser_AB   <- apply( sample_AB[ c(ii,ii+1), ], 2, min )
    Lesser_AB   <- sum(Lesser_AB)
    SumTot_Feat <- sum( sample_AB[ii,]) + sum( sample_AB[ii+1,] )
    BC_val <- c( BC_val,   1 - ((2*Lesser_AB)/ SumTot_Feat) )
    # print( paste0( (ii+1)/2, " Couple: ", rownames(sample_AB[ii, ]), " - ", 
    #                rownames(sample_AB[ii+1, ]) , " = ", Lesser_AB / SumTot_Feat ) )
  }
  res <- data.frame( BrayCurtis = BC_val )
  # We need to split rownames and reassign correct ones to the result BC df
  rownames(res) <- sapply( strsplit(rownames(df_A), '__'), "[", 1 )
  return( res )
}





# ----- WRANGLING functions --------------------------------------------------------

stop_quietly <- function( print_err ) {
  # Turn off all error messages just before calling stop(), and prevent debug. 
  # Alternative, set in R-Studio menu "Debug" -> "On Error" to "Message only". 
  # (however this reset across restarts)
  opt <- options(error = NULL)    # Stop debug intervention
  on.exit(options(opt))
  stop( print_err )
}


Match_2DF_RowCols <- function( masterDF, slaveDF ){
  # The function will reorder the "slave" DF to have the same ordering of rows
  # and cols as that of "master" DF.
  # Input:
  #   - 2 DF with same dimensions (Row, Col), and matching names for R and C 
  #     (union == intersect)
  # R and C (union == intersect). However, rows and cols may be arranged 
  # in different order.
  
  mask_cols <- colnames(slaveDF) %in% colnames(masterDF)      # Intersect cols
  mask_rows <- rownames(slaveDF) %in% rownames(masterDF)      # Intersect rows
  
  # Checks that the two DFs have same dimension and names rows and cols, 
  # albeit in a different order 
  if (dim(masterDF)[1] != dim(slaveDF)[1] ) 
    stop_quietly("Error Match_2DF_RowCols(): The two tables have different row length") 
  if (dim(masterDF)[2] != dim(slaveDF)[2] )  
    stop_quietly("Error Match_2DF_RowCols(): The two tables have different columns length") 
  if (! all(mask_cols) )
    stop_quietly("Error Match_2DF_RowCols(): Columns in slaveDF do not match masterDF") 
  if (! all(mask_rows) )   
    stop_quietly("Error Match_2DF_RowCols(): Rows in slaveDF do not match masterDF") 
  
  # Create index for slaveDF columns that matches the masterDF columns order
  idx_cols  <- match( as.vector(colnames(masterDF)), as.vector(colnames(slaveDF)) )
  slaveDF   <- slaveDF[ , idx_cols ]
  # Create index for slaveDF rows that matches the masterDF rows order
  idx_rows  <- match( as.vector(rownames(masterDF)), as.vector(rownames(slaveDF)) )
  slaveDF   <- slaveDF[ idx_rows , ]
  
  return( slaveDF )
}


MatchOrder_2DF <- function( smallDF, largeDF, Daxs) {
  # The function takes the names of rows (Daxs==1) or cols (Daxs==2) in smallDF
  # and use those to create a subset of largeDF. The result (subsetDF) should 
  # have same order or rows (or cols) as smallDF.
  # Input:
  #  - smallDF and largeDF: two dataframes. It is assumed smallDF is already 
  #                         a subset of largeDF.
  #  - Daxs: specifies whether to work on rows (=1) or columns (=2)
  # Output:
  #  - subsetDF: a subset of largeDF matching row/cols order of smallDF
  #
  # NOTE: Discrepancies in dimensions (empty rows or cols in subsetDF) occur if 
  # content in smallDF does not exist in largeDF.
  
  if (! (Daxs == 1 || Daxs == 2) )   
    stop_quietly("Error MatchOrder_2DF(): Daxs must be either 1 (rows) or 2 (columns)") 
  
  if (Daxs == 1) {
    mask_rows <- rownames(largeDF) %in% rownames(smallDF)      # Intersect rows
    largeDF   <- largeDF[ mask_rows,  ]
    # Create index for largeDF rows that matches the smallDF rows order
    idx_rows  <- match( as.vector(rownames(smallDF)), as.vector(rownames(largeDF)) )
    subsetDF  <- largeDF[ idx_rows , ]
    
  }else if (Daxs == 2) {
    mask_cols <- colnames(largeDF) %in% colnames(smallDF)      # Intersect cols
    largeDF   <- largeDF[ , mask_cols ]
    # Create index for largeDF columns that matches the smallDF columns order
    idx_cols  <- match( as.vector(colnames(smallDF)), as.vector(colnames(largeDF)) )
    subsetDF  <- largeDF[ , idx_cols ]
  }
  return( subsetDF )
}


Match_SampleID_and_Groups <- function( df_A, df_B, TimeP_1_META, TimeP_2_META) {
  # The function takes two feature tables and their Metadata table divided in
  # two groups (such as two time points: 1 and 2).
  #  1 - Remove duplicate samples, from META (1 and 2)
  #  2 - Replace ID with patient ID. Then manipulate to ensure that time points 
  #      related tables (X_1 and X_2) contain same patients in same order
  #
  # Output: is returned as a list of Data.Frames: DF_t1_A, DF_t2_A, DF_t1_B, DF_t2_B
  
  # - Part 1 - 
  # We need to remove duplicate samples, those we sequenced twice to improve yield
  # First create df returning ALL duplicated value:   ListAll_IDs %in% IDs_duplicated
  
  # Requires: library(tibble) --> tibble::add_column
  library(tibble, include.only = 'add_column')
  
  # ---> Group 1 
  # For each couple, select and exclude those which has the LOWEST total reads
  dupl_Samples <- TimeP_1_META[ duplicated(TimeP_1_META$PatientID_num) | 
                                  duplicated(TimeP_1_META$PatientID_num, fromLast=TRUE) , ]
  uniq_ID <- unique(dupl_Samples$PatientID_num)
  excl_SampleID <- NULL
  for (ii in seq(1, length(uniq_ID)) ){
    couple <- dupl_Samples[ dupl_Samples$PatientID_num == uniq_ID[ii] , ]
    excl_SampleID <- c( excl_SampleID, couple[ which.min(couple$Sample_TotReads), "Sample_ID"] [[1]])
  }
  TimeP_1_META <- TimeP_1_META[  !(TimeP_1_META$Sample_ID %in% excl_SampleID) , ]
  
  # ---> Group 2 
  # For each couple, select and exclude those which has the LOWEST total reads
  dupl_Samples <- TimeP_2_META[ duplicated(TimeP_2_META$PatientID_num) | 
                                  duplicated(TimeP_2_META$PatientID_num, fromLast=TRUE) , ]
  uniq_ID <- unique(dupl_Samples$PatientID_num)
  excl_SampleID <- NULL
  for (ii in seq(1, length(uniq_ID)) ){
    couple <- dupl_Samples[ dupl_Samples$PatientID_num == uniq_ID[ii] , ]
    excl_SampleID <- c( excl_SampleID, couple[ which.min(couple$Sample_TotReads), "Sample_ID"] [[1]])
  }
  TimeP_2_META <- TimeP_2_META[  !(TimeP_2_META$Sample_ID %in% excl_SampleID) , ]
  
  # Now divide the "timepoints" DF for each dataset A and B 
  DF_t1_A <- df_A[ rownames(df_A) %in% TimeP_1_META$Sample_ID , ]
  DF_t2_A <- df_A[ rownames(df_A) %in% TimeP_2_META$Sample_ID , ]
  DF_t1_B <- df_B[ rownames(df_B) %in% TimeP_1_META$Sample_ID , ]
  DF_t2_B <- df_B[ rownames(df_B) %in% TimeP_2_META$Sample_ID , ]
  
  # - Part 2 -
  # Move the rownames to a new column  column "Sample_ID". Replace rownames with
  # Patient ID, thus we can pair Timepoints of different datasets. 
  rownames(DF_t1_A) <- TimeP_1_META$PatientID_num
  rownames(DF_t1_B) <- TimeP_1_META$PatientID_num
  rownames(DF_t2_A) <- TimeP_2_META$PatientID_num
  rownames(DF_t2_B) <- TimeP_2_META$PatientID_num
  DF_t1_A <- DF_t1_A %>% add_column( Sample_ID=TimeP_1_META$Sample_ID, .before = 1) 
  DF_t1_B <- DF_t1_B %>% add_column( Sample_ID=TimeP_1_META$Sample_ID, .before = 1) 
  DF_t2_A <- DF_t2_A %>% add_column( Sample_ID=TimeP_2_META$Sample_ID, .before = 1) 
  DF_t2_B <- DF_t2_B %>% add_column( Sample_ID=TimeP_2_META$Sample_ID, .before = 1) 
  
  # Now we use intersect to select only patients that have two time points, in both A and B
  common_A <- intersect( rownames(DF_t1_A), rownames(DF_t2_A) )
  DF_t1_A  <- DF_t1_A[ rownames(DF_t1_A) %in% common_A, ]
  DF_t2_A  <- DF_t2_A[ rownames(DF_t2_A) %in% common_A, ]
  common_B <- intersect( rownames(DF_t1_B), rownames(DF_t2_B) )
  DF_t1_B  <- DF_t1_B[ rownames(DF_t1_B) %in% common_B, ]
  DF_t2_B  <- DF_t2_B[ rownames(DF_t2_B) %in% common_B, ]
  
  # Lastly we ensure that the two groups of each table (A and B) have the same order of rows and cols 
  DF_t1_A <- MatchOrder_2DF( DF_t2_A, DF_t1_A, 1)
  DF_t1_A <- MatchOrder_2DF( DF_t2_A, DF_t1_A, 2)
  DF_t1_B <- MatchOrder_2DF( DF_t2_B, DF_t1_B, 1)
  DF_t1_B <- MatchOrder_2DF( DF_t2_B, DF_t1_B, 2)
  
  res <- list( DF_t1_A, DF_t2_A, DF_t1_B, DF_t2_B )
  
  return( res )
}




# ----- PHYLOSEQ related functions ---------------------------------------------

sangerTables_to_Phyloseq <- function( ASV_df, Taxo_df, Meta_df ){
  # Transform three tables that represent a feature table, a taxonomic nomenclature
  # table and the sample's metadata into formats that are compatible to be 
  # imported to generate a phyloseq object
  #
  # Input: 
  #   - feature tables (ASV_df)
  #   - taxonomic table (Taxo_df)
  #   - metadata (Meta_df) 
  # These are passed in our standard format (sanger scripts). The variable are
  # then manipulated and transformed to generate...
  # Output:
  #   - ASV (=OTU) table (ps_OTU)  as otu_table()
  #   - Taxonomic table  (ps_TAXO) as tax_table()
  #   - Metadata table   (ps_SAMP) as sample_data()
  # These the necessary variables that are ready to be combined into a phyloseq 
  # object:   >> phyloseq( ps_OTU, ps_TAXO, ps_SAMP )
  
  # --- OTU table ---> as data.matrix()
  t_Otu     <- ASV_df
  keepNames <- ASV_df$Sample_ID
  t_Otu     <- t_Otu %>% dplyr::select( -c("Sample_ID") )
  rownames( t_Otu ) <- keepNames
  t_Otu  <- data.matrix( t(t_Otu) )
  ps_OTU   <- otu_table( t_Otu, taxa_are_rows = TRUE )
  
  # --- TAXOnomic table ---> as data.matrix()
  t_Tax     <- Taxo_df
  keepNames <- colnames(t_Tax)[-1]
  t_Tax     <- t_Tax[-c(1,2,3),-c(1)] %>% t()
  # Split all Taxo names to remove the prefix "D_Nn__" from the names
  for ( ii in seq(1,dim(t_Tax)[2]) ){
    temp <- as.data.frame( sapply( strsplit( t_Tax[,ii], "__" ), head, 2) )
    t_Tax[,ii] <- as.character(as.vector( temp[2,] ))
  }
  t_Tax  <- data.matrix(t_Tax)
  ps_TAXO <- tax_table( t_Tax )
  # NOTE: Function tax_table() reindex rows and cols. Reassign here to ensure correct names
  taxa_names(ps_TAXO) <- keepNames
  colnames(ps_TAXO)   <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  
  # --- SAMPle table ---> as data.frame() 
  t_Sam <- Meta_df
  t_Sam <- as.data.frame(t_Sam)
  rownames(t_Sam) <- ASV_df$Sample_ID
  ps_SAMP  <- sample_data( data.frame( t_Sam, row.names=ASV_df$Sample_ID, stringsAsFactors=FALSE ))
  
  # --- Return variables ready to join into a PHYLOSEQ object
  return( list( ps_OTU, ps_TAXO, ps_SAMP ) )
  
}




# ----- Grouping by RANK for compositional plots -------------------------------

groupBy_1_Rank <- function( ori_ASV, ori_TAXO, sele_Rank, abund_threshold ) { 
  # The aim of the function is to simplify the dataset, by reducing the number 
  # of features by combining those that have same taxonomy at a specific rank.
  #
  # First, classify unique taxa in a specific Rank (sele_Rank) and aggregate
  # all the features belonging to the same taxonomic "name".
  # Second, those that have abundance below a target level (abund_threshold)
  # will be combined into a new group ("Others"), simplifying the dataset even 
  # further
  #
  # Input:
  # - ori_ASV:   count table; rows are observations, and cols are features
  # - ori_TAXO:  taxonomy table; rows are observ., and cols are taxonomic ranks
  # - sele_Rank:  Rank to use for grouping features
  # - abund_threshold:  minimum abundance (proportion). All taxa below this are 
  #                     collected and summed together as "others". Those above 
  #                     are kept as separate features 
  # Output:
  # - grouped_ASV:   grouped ori_ASV
  # - grouped_TAXO:  grouped ori_TAXO
  # - pivotLong:  ASV data in long format, plus taxonomic ranks as separate 
  #               columns. Ready to use for bar-plots.
  
  # Combine (and pivot ASV) to have features as rows
  original_DF <- cbind( ori_ASV, ori_TAXO )
  # original_DF <- original_DF[ order(original_DF[,sele_Rank]) , ]
  
  grouped_DF  <- NULL      # Store reorganized dataset
  
  # COOL PIPE: classify taxa in (column) sele_Rank. Build summary table of unique 
  # taxa in (column) sele_Rank, sum number of counts and order from most to lowest abundant.
  reviewRank <- original_DF %>% 
                group_by( .data[[sele_Rank]] ) %>% 
                summarise( across(where(is.numeric), list( ColTot = sum))) %>%
                mutate( "Total" = rowSums( across( 2:dim(.)[2] )), .keep = "unused")
  reviewRank       <- reviewRank[ order( - reviewRank$Total), ]
  reviewRank$Total <- round( reviewRank$Total /sum(reviewRank$Total), 3 )
  
  # ---- Classify by Upp/Low Ranks ----------------------------------------------------
  # Split sele_Rank taxas in two groups:
  mask <- reviewRank$Total >= abund_threshold
  Taxa_distinct <- as.list( reviewRank[ mask , sele_Rank ] )[[1]]
  Taxa_combine  <- as.list( reviewRank[ !mask, sele_Rank ] )[[1]]
  
  # A1 - Taxa_distinct - those with at least abundancy above abund_threshold (i.e. 5%)
  #      These will be further divide into taxas sub-groups based on sele_Rank.
  if ( length(Taxa_distinct) > 0 ){
    # Group all ii-th taxa in Taxa_distinct as separated features
    for (ww in seq(1, length(Taxa_distinct))) {
      DF_1Rank <- original_DF[ original_DF[[sele_Rank]] == Taxa_distinct[ww] , ]
      temp  <- DF_1Rank[1,]
      temp[ , 1:(dim(DF_1Rank)[2]-7) ] <- colSums( DF_1Rank[ , 1:(dim(DF_1Rank)[2]-7) ] ) 
      grouped_DF <- rbind( grouped_DF, temp) 
    }
  } else if ( length(Taxa_distinct) == 0 ){
    # None in sele_Rank with abundancy > abund_threshold; then combine them all into one
    DF_1Rank <- original_DF[ original_DF[[sele_Rank]] %in% Taxa_combine , ]
    temp  <- DF_1Rank[1,]
    temp[ , 1:(dim(DF_1Rank)[2]-7) ] <- colSums( DF_1Rank[ , 1:(dim(DF_1Rank)[2]-7) ] ) 
    temp[ , sele_Rank] <- paste0( Taxa_distinct[uu], " (All)")       # Rename the sele_Rank
    grouped_DF <- rbind( grouped_DF, temp)
  }
  
  # A2 - Taxa_combine - those with low abundacy (< abund_threshold)
  #      These will be combined into one group, named  "(others)"
  if ( length(Taxa_combine) != 0 && length(Taxa_distinct) > 0 ){
    # Combined into one group all sele_Rank with low abundancy 
    reviewUpR <- original_DF[ original_DF[[sele_Rank]] %in% Taxa_combine , ]
    temp  <- DF_1Rank[1,]
    temp[ , 1:(dim(DF_1Rank)[2]-7) ] <- colSums( DF_1Rank[ , 1:(dim(DF_1Rank)[2]-7) ] ) 
    temp[ , sele_Rank] <- paste0( sele_Rank, " (others)")     # Rename the UppRank
    grouped_DF <- rbind( grouped_DF, temp)
  }
  
  # ---- Pivot to Long Format ----------------------------------------------------
  # We have to split the grouped_DF into its original ASV and TAXO, ... 
  grouped_ASV  <- as.data.frame(t( grouped_DF[ , 1 : (dim(grouped_DF)[2]-7) ] ))
  grouped_TAXO <- grouped_DF[ , (dim(grouped_DF)[2]-6) : dim(grouped_DF)[2] ] 
  # ... thus be able to convert to long format 
  grouped_ASV$Sample_ID <- rownames(grouped_ASV)              # NOTE: add as last columns
  pivotLong <- gather( grouped_ASV, FeatureName, Counts, 
                       1:dim(grouped_ASV)[2]-1, factor_key=TRUE)
  # Initialize empty columns
  pivotLong$Species <- pivotLong$Genus <- pivotLong$Family <- 
    pivotLong$Order <- pivotLong$Class <- pivotLong$Phylum <- ""
  # Add rank names to the long format
  for ( ii in seq(1,dim(grouped_TAXO)[1]) ){
    mask <- pivotLong$FeatureName == rownames(grouped_TAXO)[ii]
    pivotLong$Phylum[ mask ] <- grouped_TAXO[ii, "Phylum"]
    pivotLong$Class[  mask ] <- grouped_TAXO[ii, "Class"]
    pivotLong$Order[  mask ] <- grouped_TAXO[ii, "Order"]
    pivotLong$Family[ mask ] <- grouped_TAXO[ii, "Family"]
    pivotLong$Genus[  mask ] <- grouped_TAXO[ii, "Genus"]
    pivotLong$Species[ mask] <- grouped_TAXO[ii, "Species"]
  }
  
  return( list( grouped_ASV, grouped_TAXO, pivotLong ) )
} # END groupBy_1_Ranks




groupBy_2_Ranks <- function( ori_ASV, ori_TAXO, UppRank, LowRank, 
                             threshold_UP, threshold_LW ) { 
  # The aim of the function is to simplify the dataset (as in groupBy_1_Ranks)
  # This is done by grouping at two levels:
  # A - It will try to collect all the features belonging to unique taxa at
  #     a specific upper rank (e.g. Class: Clostridia, Bacilli, ...) 
  # B - Then, features at a specific i-th upper rank (say, Clostridia) will be
  #     further classified into smaller group using a lower rank 
  #     (e.g. Families: Ruminococcaceae, Lachnospiraceae,...) 
  # C - At each level, features with low abundance will be collected and 
  #     combined properly (for example, as either "Class (all others)" or 
  #     "Clostridia (others)")
  #
  # Input:
  # - ori_ASV:   count table; rows are observations, and cols are features
  # - ori_TAXO:  taxonomy table; rows are observ., and cols are taxonomic ranks
  # - UppRank:   upper Rank to use for first level grouping 
  # - LowRank:   lower Rank to use for second level grouping 
  # - threshold_UP / _DW:  minimum abundance (proportion), at either upper or 
  #              lower rank. All taxa below this are collected and summed 
  #              together as "others". Those above are kept as separate features 
  # Output:
  # - grouped_ASV:   grouped ori_ASV
  # - grouped_TAXO:  grouped ori_TAXO
  # - pivotLong:  ASV data in long format, plus taxonomic ranks as separate 
  #               columns. Ready to use for bar-plots.
  
  # Combine (and pivot ASV) to have features as rows
  original_DF <- cbind( ori_ASV, ori_TAXO )
  # original_DF <- original_DF[ order( original_DF[,UppRank], original_DF[,LowRank] ), ]
  grouped_DF  <- NULL      # Store reorganized dataset
  
  # COOL PIPE: classify taxa in (column) UppRank. Build summary table of unique 
  # taxa in (column) UppRank, sum number of counts and order from most to lowest abundant.
  reviewUpR <- original_DF %>% 
                group_by( .data[[UppRank]] ) %>% 
                summarise( across(where(is.numeric), list( ColTot = sum))) %>%
                mutate( "Total" = rowSums( across( 2:dim(.)[2] )), .keep = "unused")
  reviewUpR <- reviewUpR[ order( - reviewUpR$Total), ]
  reviewUpR$Total <- round( reviewUpR$Total /sum(reviewUpR$Total), 3 )
  
  
  # ---- Classify by Upp/Low Ranks ----------------------------------------------------
  # Split UppRank taxas in two groups:
  mask <- reviewUpR$Total >= threshold_UP
  Upp_distinct <- as.list( reviewUpR[ mask , UppRank ] )[[1]]
  Upp_combine  <- as.list( reviewUpR[ !mask, UppRank ] )[[1]]
  
  # A1 - Upp_distinct - those with at least abundancy above threshold_UP (i.e. 5%)
  #      These will be further divide into taxas sub-groups based on LowRank.
  for (uu in seq(1, length(Upp_distinct))) {
    
    # Subset of original_DF containing only features that belong to the distinct
    # ii-th UppRank
    DF_1UpRank <- original_DF[ original_DF[[UppRank]] == Upp_distinct[uu], ]
    # RowSum all observations (columns of original_DF) 
    DF_1UpRank <- cbind( data.frame( Total = rowSums( DF_1UpRank[ , 1:(dim(DF_1UpRank)[2]-7) ] ) ) ,
                         DF_1UpRank[ , (dim(DF_1UpRank)[2]-6):dim(DF_1UpRank)[2] ]
    )  # or try mean() (excluding zeros) !!!
    
    if ( sum(DF_1UpRank$Total)!= 0){
      # COOL PIPE: Classify taxas in LowRank that belongs to the ii-th upper rank
      reviewLwR <- DF_1UpRank %>% 
                    group_by(.data[[LowRank]]) %>% 
                    summarise( across(where(is.numeric), list( sumRows = sum))) %>%
                    mutate(  "Total" = rowSums( across( 2:dim(.)[2] )), .keep = "unused")
      reviewLwR <- reviewLwR[ order( - reviewLwR$Total), ]
      reviewLwR$Total <- reviewLwR$Total /sum(reviewLwR$Total)
      
      # Split LowRank taxas in two groups:
      mask <- reviewLwR$Total >= threshold_LW
      Low_distinct <- as.list(reviewLwR[ mask,  LowRank ] )[[1]]
      Low_combine  <- as.list(reviewLwR[ !mask, LowRank ] )[[1]]
      
      # Move the unassigned to the combined
      if (any(Low_distinct %in% "Unassigned")){
        Low_distinct <- Low_distinct[ !Low_distinct %in% "Unassigned" ]
        Low_combine  <- c(Low_combine, Low_distinct[ Low_distinct %in% "Unassigned" ])
      }
      if (any(Low_distinct %in% "unknown")){
        Low_distinct <- Low_distinct[ !Low_distinct %in% "unknown" ]
        Low_combine  <- c(Low_combine, Low_distinct[ Low_distinct %in% "unknown" ])
      }
      
      # B1 - Low_distinct - Group features that have high abundance as distinct 
      #      feature of their own in grouped_DF
      if ( length(Low_distinct) > 0 ){
        # Group all ii-th taxa in Low_distinct as separated features
        for (ww in seq(1, length(Low_distinct))) {
          DF_1LwRank <- original_DF[ original_DF[[LowRank]] == Low_distinct[ww] , ]
          temp  <- DF_1LwRank[1,]
          temp[ , 1:(dim(DF_1LwRank)[2]-7) ] <- colSums( DF_1LwRank[ , 1:(dim(DF_1LwRank)[2]-7) ] ) 
          grouped_DF <- rbind( grouped_DF, temp) 
        }
      } else if ( length(Low_distinct) == 0 ){
        # No LowRank with abundancy > threshold_LW; then combine them all into one
        # Select using UpRank subset, then use rownames to obtain all lower ranks rows from original_DF
        subgroupLW <- DF_1UpRank[ DF_1UpRank[[LowRank]] %in% Low_combine , ]  
        DF_1LwRank <- original_DF[ rownames(subgroupLW), ] 
        temp  <- DF_1LwRank[1,]
        temp[ , 1:(dim(DF_1LwRank)[2]-7) ] <- colSums( DF_1LwRank[ , 1:(dim(DF_1LwRank)[2]-7) ] ) 
        temp[ , LowRank] <- paste0( Upp_distinct[uu], " (All)")       # Rename the LowRank
        grouped_DF <- rbind( grouped_DF, temp)
      }
      
      # B2 - Low_combine - Combined into one group all LowRank with low 
      #      abundancy (< threshold_LW)  
      if ( length(Low_combine) != 0 && length(Low_distinct) > 0 ){
        # Combined into one group all LowRank with low abundancy (< threshold_LW)  
        # Select using UpRank subset, then use rownames to obtain all lower ranks rows from original_DF
        subgroupLW <- DF_1UpRank[ DF_1UpRank[[LowRank]] %in% Low_combine , ]  
        DF_1LwRank <- original_DF[ rownames(subgroupLW), ] 
        temp  <- DF_1LwRank[1,]
        temp[ , 1:(dim(DF_1LwRank)[2]-7) ] <- colSums( DF_1LwRank[ , 1:(dim(DF_1LwRank)[2]-7) ] ) 
        temp[ , LowRank] <- paste0( Upp_distinct[uu], " (others)")     # Rename the LowRank
        grouped_DF <- rbind( grouped_DF, temp)
      }
    }
  }
  
  # A2 - Upp_combine  - those abundancy below threshold_UP
  #      These will be combined into one group, named "UppRank, (All others)"
  if ( length(Upp_combine) != 0 && length(Upp_distinct) > 0 ){
    # Combined into one group all UppRank with low abundancy (< threshold_UP) 
    reviewUpR <- original_DF[ original_DF[[UppRank]] %in% Upp_combine , ]
    temp  <- reviewUpR[1,]
    temp[ , 1:(dim(reviewUpR)[2]-7) ] <- colSums( reviewUpR[ , 1:(dim(reviewUpR)[2]-7) ] ) 
    # Rename the UppRank, as well as LowRank, since grouping for barplot will be done at the LowRank
    temp[ , UppRank] <- paste0( UppRank, " (All others)")
    temp[ , LowRank] <- paste0( UppRank, " (All others)") 
    grouped_DF <- rbind( grouped_DF, temp)
  }
  
  
  # ---- Pivot to Long Format ----------------------------------------------------
  # We have to split the grouped_DF into its original ASV and TAXO, ... 
  grouped_ASV  <- as.data.frame(t( grouped_DF[ , 1 : (dim(grouped_DF)[2]-7) ] ))
  grouped_TAXO <- grouped_DF[ , (dim(grouped_DF)[2]-6) : dim(grouped_DF)[2] ] 
  # ... thus be able to convert to long format 
  grouped_ASV$Sample_ID <- as.character( rownames(grouped_ASV) )             # NOTE: add as last columns
  pivotLong <- gather( grouped_ASV, FeatureName, Counts, 
                       1:dim(grouped_ASV)[2]-1, factor_key=TRUE)
  # Initialize empty columns
  pivotLong$Species <- pivotLong$Genus <- pivotLong$Family <- 
    pivotLong$Order <- pivotLong$Class <- pivotLong$Phylum <- ""
  # Add rank names to the long format
  for ( ii in seq(1,dim(grouped_TAXO)[1]) ){
    mask <- pivotLong$FeatureName == rownames(grouped_TAXO)[ii]
    pivotLong$Phylum[ mask ] <- grouped_TAXO[ii, "Phylum"]
    pivotLong$Class[  mask ] <- grouped_TAXO[ii, "Class"]
    pivotLong$Order[  mask ] <- grouped_TAXO[ii, "Order"]
    pivotLong$Family[ mask ] <- grouped_TAXO[ii, "Family"]
    pivotLong$Genus[  mask ] <- grouped_TAXO[ii, "Genus"]
    pivotLong$Species[ mask] <- grouped_TAXO[ii, "Species"]
  }
  
  return( list( grouped_ASV, grouped_TAXO, pivotLong ) )   # RES_N <- listDF[[N]]
  
} # END groupBy_2_Ranks



combineReady_ASV_TAXO <- function( DF_asv, DF_taxo, NameSamples ) {
  # Take two tables (ASV and Taxo), transpose and adjust to ensure that they 
  # both have the same number of features as rows.
  # DF_asv  - ASV count table (as N-obs by X-feat) 
  # DF_taxo - Taxonomic table (as 10 by X-feat)
  # Return same data transpose and with all taxonomic names stripped of prefix
  # (they are received in form: D__X__name)
  DF_asv <- as.data.frame( t(DF_asv) )
  colnames(DF_asv) <- as.character( NameSamples )
  
  NameRows <- DF_taxo$RowName
  DF_taxo  <- DF_taxo %>% dplyr::select( -c("RowName") )
  DF_taxo  <- as.data.frame( t(DF_taxo) )
  colnames(DF_taxo) <- NameRows
  
  # Split all Taxo names to remove the prefix "D_RR__" from the names
  for ( ii in seq(1,dim(DF_taxo)[2]) ){
    temp <- as.data.frame( sapply( strsplit( DF_taxo[,ii], "__" ), head, 2) )
    DF_taxo[,ii] <- as.character(as.vector( temp[2,] ))
  }
  DF_taxo <- DF_taxo %>% dplyr::select( -c("Region", "Confidence", "ConsensusLineage") )
  return( list(DF_asv, DF_taxo) )
}



# ----- Color by RANK for compositional plots -------------------------------

ColourMultiGradient <- function( dataset, group, subgroup, CustPalette=0 ,
                                 LightenOthers=FALSE, LightenAllOthers=FALSE ) {
  # The function creates a custom color palette, returned as a list of
  # hexadecimal values.
  # The function finds how many base-colors to create (=unique(groups)). Then, 
  # it take a range of the graded tones equal to the number of "subgroup", using 
  # one column of CustPalette_X, which contains the tonalities for a specific 
  # base color
  #
  # Input:
  # - dataset = the output[3] returned by function groupBy_2_Ranks.
  #             long format compositional dataset 
  #      rows : all information about a specific feature for one sample
  #      cols : are SampleID, FeatureID, total counts for a feature in specific 
  #             sample; lastly 6 columns with the taxonomic names at each rank
  # - group    = rank to used for first-level grouping, 
  #              determine base color 
  # - subgroup = rank to used for second-level grouping, 
  #              define number of tones to create (one for each subgroup)
  # - CustPalette   = Provide custom palette [OPTION not implemented yet]
  # - LightenOthers = choice to assign lightest tone color for subgroups that 
  #                   have tag "other"
  # - LightenAllOthers = choice to assign Gray color for subgroups that 
  #                      have tag "All others"
  # Output:
  # - colorScale = a vector with colors, which have same length as
  #                unique( dataset[, subgroup]).
  #                Colors are arranged to match the order in columns subgroup
  #                of dataset.
  
  # Define defaults color palattes
  if (CustPalette==0) {   
    # A - set defaults 12x15 palette (15 base-colors in 12 graded tones)
    CustPalette_A  <- data.frame(
      "CyanBlue"   =c("#1781B2","#1B88BB","#1F8FC3","#2495CB","#299CD2","#31A7DC","#39B2E4",
                      "#43BDEC","#4DC7F3","#59D1F9","#66D9FD","#77E2FF","#88E9FF") ,
      "YellowLight"=c("#A18E59","#B29C5E","#C3AA62","#D5B865","#E8C668","#FBD56A","#FFDA71",
                      "#FFDF7A","#FFE385","#FFE692","#FFE9A2","#FFEDB4","#FFF1C8") ,
      "RedLight"   =c("#A3615A","#B4675E","#C66C62","#D97165","#EC7568","#FF796A","#FF7C71",
                      "#FF8179","#FF8984","#FF9292","#FFA1A3","#FFB3B7","#FFC7CC") ,
      "GreenDark"  =c("#567450","#5C8055","#628C5A","#68985F","#6EA463","#73B167","#81BE73",
                      "#8FCA7F","#9ED58D","#ADDF9C","#BDE7AC","#CCEFBD","#DCF5D0") ,
      "Pinky"      =c("#925E7F","#A1638B","#B16996","#C06EA2","#D172AE","#E176BA","#EB7EC6",
                      "#F387D1","#FA93DC","#FF9FE6","#FFAEEE","#FFBEF5","#FFCFFA") ,
      "BlueDark"   =c("#565FA2","#5A65B4","#5D6AC6","#606ED8","#6272EB","#6476FF","#6B81FF",
                      "#748EFF","#7F9CFF","#8DABFF","#9DBAFF","#B0CBFF","#C5DBFF") ,
      "Orange"     =c("#A3775A","#B4805E","#C68A62","#D99365","#EC9C68","#FFA56A","#FFA971",
                      "#FFAE79","#FFB384","#FFBA92","#FFC2A1","#FFCCB3","#FFD8C7") ,
      "RedDark"    =c("#785159","#845660","#905A67","#9D5F6D","#AA6373","#B76779","#C47287",
                      "#CF7F95","#DA8CA3","#E39BB2","#EBABC1","#F2BDD0","#F7CFDF") ,
      "GreenLight" =c("#55995E","#59AA63","#5CBB68","#5FCC6D","#62DE72","#64F076","#6CFB7B",
                      "#76FF81","#82FF89","#90FF94","#A0FFA1","#B4FFB2","#C9FFC7") ,
      "Purple"     =c("#7A5F96","#8464A5","#8F69B5","#996EC5","#A472D6","#AE76E7","#B47EF1",
                      "#BA87F9","#C092FF","#C79EFF","#CFADFF","#D7BDFF","#E1CFFF") ,
      "BlueLight"  =c("#5A7EA3","#5F89B4","#6394C6","#669FD9","#69A9EC","#6BB4FF","#71C0FF",
                      "#7ACBFF","#85D6FF","#92E1FF","#A2EAFF","#B4F2FF","#C8F8FF") ,
      "YellowOlive"=c("#8A8A55","#98985A","#A7A75F","#B6B663","#C5C567","#D5D56A","#E0DE73",
                      "#EAE57E","#F2EB8A","#F9EF98","#FDF3A8","#FFF6B9","#FFF9CC") ,
      "CyanGreen"  =c("#579B85","#5BAC92","#5FBD9F","#62CEAC","#65E1B9","#67F3C6","#6EFDCC",
                      "#78FFD1","#84FFD6","#91FFDB","#A1FFE0","#B3FFE5","#C8FFEB") ,
      "OrangeDark" =c("#896D5B","#977660","#A57F65","#B4886A","#C3906F","#D29973","#DDA07C",
                      "#E6A886","#EEB092","#F5BAA0","#FAC4AE","#FECFBF","#FFDBD0") ,
      "Grey"       =c("#887979","#908181","#988989","#A09191","#A8999A","#AFA1A2","#B7A9AA",
                      "#BFB1B3","#C6BABB","#CDC2C3","#D5CACC","#DCD3D4","#E3DCDD")
    )
    # B - set defaults 21x14 palette (14 base-colors in 21 graded tones)
    CustPalette_B  <- data.frame(
      "CyanBlue"   =c("#09638D","#0C6B97","#0F72A0","#1379A9","#1781B2","#1B88BB","#1F8FC3",
                      "#2495CB","#299CD2","#31A7DC","#39B2E4","#43BDEC","#4DC7F3","#59D1F9",
                      "#66D9FD","#77E2FF","#88E9FF","#99EFFF","#AAF5FF","#BBF9FF","#CCFDFF") ,
      "YellowLight"=c("#996F00","#A67B02","#B38705","#BF9408","#CCA00C","#D9AD10","#E5B915",
                      "#EEC61A","#F6D21F","#FFD826","#FFDC2D","#FFE036","#FFE444","#FFE755",
                      "#FFE966","#FFEC77","#FFEE88","#FFF099","#FFF2AA","#FFF4BB","#FFF6CC") ,
      "RedLight"   =c("#993252","#A63855","#B33F57","#BF4559","#CC4C5A","#D9535C","#E65A5E",
                      "#F26562","#FF736A","#FF776E","#FF7C73","#FF8279","#FF8880","#FF9087",
                      "#FF978F","#FFA098","#FFA9A2","#FFB3AD","#FFBEB9","#FFCAC5","#FFD6D2") ,
      "GreenDark"  =c("#297127","#2E792B","#328130","#378935","#3C913A","#41983F","#47A044",
                      "#4CA749","#52AE4F","#5BB757","#65C060","#70C869","#7ACF72","#86D67D",
                      "#91DD87","#9DE392","#A9E89E","#B5EDAA","#C1F1B7","#CDF5C5","#DAF8D2") ,
      "Pinky"      =c("#664B99","#7353A6","#805BB3","#8E63BF","#9C6BCC","#AA74D9","#B87CE6",
                      "#C785F0","#D58EFB","#D691FF","#D794FF","#D798FF","#D89DFF","#D8A2FF",
                      "#D9A7FF","#DAAEFF","#DBB4FF","#DDBCFF","#E0C4FF","#E3CCFF","#E7D5FF") ,
      "BluesDark"  =c("#354599","#3B4CA6","#4153B3","#485BBF","#4F62CC","#566AD9","#5D72E4",
                      "#657AEE","#6D82F8","#7187FE","#768CFF","#7C91FF","#8398FF","#8A9EFF",
                      "#93A6FF","#9CADFF","#A5B6FF","#B0BEFF","#BBC8FF","#C7D2FF","#D4DCFF") ,
      "Orange"     =c("#994F11","#A65915","#B36419","#BF6F1D","#CC7A22","#D98527","#E6912D",
                      "#EF9C33","#F8A839","#FFAC3F","#FFB146","#FFB54E","#FFB957","#FFBD61",
                      "#FFC16C","#FFC678","#FFCB88","#FFD099","#FFD6AA","#FFDDBB","#FFE4CC") ,
      "RedDark"    =c("#6D2450","#792956","#862F5C","#923561","#9D3C66","#A9426A","#B4496E",
                      "#BE5172","#C95876","#D15E7E","#D86587","#DF6C8F","#E57498","#EA7CA1",
                      "#EF85AA","#F38FB3","#F799BC","#FAA3C5","#FCAFCE","#FEBBD7","#FFC7E0") ,
      "CyanGreen"  =c("#3D9164","#449B6E","#4AA678","#51B083","#58BA8D","#5FC498","#67CEA3",
                      "#6ED8AE","#76E1B9","#7BE7BD","#81EDC1","#87F2C5","#8EF6C9","#95FACD",
                      "#9DFDD1","#A5FFD5","#AEFFD9","#B8FFDE","#C3FFE2","#CDFFE7","#D9FFED") ,
      "Purple"     =c("#724C99","#7D54A6","#885BB2","#9463BD","#A06CC9","#AB74D3","#B77CDE",
                      "#C385E9","#CF8EF3","#D292F8","#D596FC","#D79BFF","#DAA0FF","#DCA6FF",
                      "#DFACFF","#E2B3FF","#E4BBFF","#E7C3FF","#EBCBFF","#EED5FF","#F2DEFF") , 
      "BlueLight"  =c("#617299","#6A7BA6","#7385B3","#7D90BF","#869ACC","#90A4D9","#9AAEE6",
                      "#A4B9F2","#AEC3FF","#B0C6FF","#B3C8FF","#B6CBFF","#BACEFF","#BED2FF",
                      "#C2D5FF","#C7D9FF","#CDDDFF","#D3E1FF","#D9E6FF","#E0EAFF","#E7EFFF") ,
      "YellowDark" =c("#7B540B","#835D0E","#8B6611","#936F14","#9B7817","#A2811B","#AA8B1F",
                      "#B09424","#B79D28","#C1A531","#CBAD3A","#D4B445","#DCBB50","#E3C15C",
                      "#E9C769","#EFCD77","#F4D388","#F8D999","#FBDFAA","#FDE5BB","#FFEBCC") ,
      "PalePink"   =c("#99637A","#A66C82","#B37689","#BF7F90","#CC8998","#D9939F","#E69DA6",
                      "#F2A7AD","#FFB1B4","#FFB3B4","#FFB7B6","#FFBCB9","#FFC2BC","#FFC8C0",
                      "#FFCEC5","#FFD4C9","#FFDACF","#FFE0D4","#FFE6DA","#FFEBE1","#FFF1E8") ,
      "GreenLight" =c("#55993E","#59A445","#5DAF4C","#60BA53","#63C55A","#67CF61","#6AD969",
                      "#71E375","#79ED82","#7DF383","#81F884","#87FC86","#90FF8C","#99FF92",
                      "#A3FF99","#ADFFA1","#B7FFA9","#C1FFB1","#CBFFBB","#D4FFC5","#DEFFCF") , 
      "Gray"       =c("#757C81","#7F868C","#898F96","#9399A1","#9EA3AB","#A8ACB6","#B2B6C0",
                      "#BCBFCB","#C6C9D5","#C9CCD8","#CCD0DB","#CFD4DE","#D2D7E1","#D6DBE4",
                      "#D9DEE7","#DDE2EA","#E0E5ED","#E4E9EF","#E7ECF2","#EBF0F4","#EFF3F6")
    )
  }
  # Create summary of unique lower ranks and total read counts. In effect, this 
  # is a list of all elements to plot in a stacked bar-plot
  subsetDF <- dataset[ , c( group, subgroup, "Counts")]
  summ_LW  <- data.frame()
  uniq_Lw  <- unique( subsetDF[,subgroup] )
  for ( ii in seq(1, length(uniq_Lw)) ){
    mask    <- subsetDF[,subgroup] == uniq_Lw[ii]
    summ_LW <- rbind( summ_LW, c( subsetDF[ which(subsetDF[,subgroup]==uniq_Lw[ii])[1], group ], 
                                  uniq_Lw[ii],   sum(subsetDF[ mask,"Counts"]) ) )
  }
  colnames(summ_LW)   <- c( group, subgroup, "Counts" )
  summ_LW[,"Counts" ] <- as.numeric( summ_LW[,"Counts" ] )
  # Create the same for the upper ranks, and counting the number of lower ranks
  # present in each upper rank (e.g. how many Families there are in each Class)
  summ_UP <- data.frame()
  uniq_UP <- unique( summ_LW[,group] )
  for ( uu in seq(1, length(uniq_UP)) ){
    mask    <- summ_LW[,group] == uniq_UP[uu]
    summ_UP <- rbind( summ_UP, c( uniq_UP[uu] , sum( mask*1 )) )
  }
  colnames(summ_UP) <- c( group, "N_subgroups" )
  summ_UP[ ,2] <- as.numeric( summ_UP[,2] )
  
  # Set tones dynamically -  the first tone color is chosen such to centered the 
  # scale around the middle of the available range, so we don't have too many dark colors
  colorScale <- NULL
  for (cc in seq(1, nrow( summ_UP )) ) {
    # According to the number of lower ranks in the cc-th upper ranks, ...
    if      (summ_UP[ cc, 2] <= 10)       CustPalette <- CustPalette_A  # ...use short range palette
    else if (summ_UP[ cc, 2]  > 10)       CustPalette <- CustPalette_B  # ...or use long range palette
    
    firstClr  <- round(nrow(CustPalette)/2) - round(summ_UP[ cc,2] / 2 )
    if (firstClr <= 0)         firstClr <- 1
    lastClr   <- summ_UP[ cc, 2 ]
    toneRange <- CustPalette[ firstClr:(firstClr-1+lastClr) , cc]
    
    # Give the lightest tone color to the last, if lower rank is tagged "other"  
    temp <- summ_LW[ summ_LW[,1]==summ_UP[cc,1], ]
    if (LightenOthers){
      if ( grepl("other", temp[ dim(temp)[1],2 ]) )     
        toneRange[length(toneRange)] <- CustPalette[ dim(CustPalette)[1], cc]
    }
    colorScale <- c(colorScale, toneRange ) 
  }
  
  # Assign last palette color as Grey to the upper rank grouping, if tagged "All others" 
  if (LightenAllOthers){
    if ( grepl("All other", summ_UP[ dim(summ_UP)[1], 1]) )  
      colorScale[length(colorScale)] <- CustPalette_B$Gray[ round( dim(CustPalette_B)[1]/2 ) ]
  }
  
  return(colorScale)
  
}




# ----- Color by RANK for compositional plots -------------------------------
ColourMultiGradient_Metabolites <- function( dataset, group, subgroup, CustPalette=0 ,
                                 LightenOthers=TRUE, LightenAllOthers=TRUE, 
                                 dynamicToneRange=FALSE) {
  # The function creates a custom color palette, returned as a list of
  # hexadecimal values.
  # The function finds how many base-colors to create (=unique(groups)). Then, 
  # it take a range of the graded tones equal to the number of "subgroup", using 
  # one column of CustPalette_X, which contains the tonalities for a specific 
  # base color
  #
  # Input:
  # - dataset = the output[3] returned by function groupBy_2_Ranks.
  #             long format compositional dataset 
  #      rows : all information about a specific feature for one sample
  #      cols : are SampleID, FeatureID, total counts for a feature in specific 
  #             sample; lastly 6 columns with the taxonomic names at each rank
  # - group    = rank to used for first-level grouping, 
  #              determine base color 
  # - subgroup = rank to used for second-level grouping, 
  #              define number of tones to create (one for each subgroup)
  # - CustPalette   = Provide custom palette [OPTION not implemented yet]
  # - LightenOthers = choice to assign lightest tone color for subgroups that 
  #                   have tag "other"
  # - LightenAllOthers = choice to assign Gray color for subgroups that 
  #                      have tag "All others"
  # Output:
  # - colorScale = a vector with colors, which have same length as
  #                unique( dataset[, subgroup]).
  #                Colors are arranged to match the order in columns subgroup
  #                of dataset.
  
  # Define defaults color palattes
  library(readr)
  
  if (CustPalette==0) {   
    # A - set defaults 15x17 palette (15 graded tones in 17 base-colors)
    CustPalette_B <- read_tsv("/Users/mattesa/ZenBook/R/LUT_Metabolites.txt", skip_empty_rows = TRUE ,show_col_types = FALSE)
  }
  # Create summary of unique lower ranks and total read counts. In effect, this 
  # is a list of all elements to plot in a stacked bar-plot
  subsetDF <- dataset[ , c( group, subgroup, "Counts")]
  summ_LW  <- data.frame()
  uniq_Lw  <- unique( subsetDF[,subgroup] )
  for ( ii in seq(1, length(uniq_Lw)) ){
    mask    <- subsetDF[,subgroup] == uniq_Lw[ii]
    summ_LW <- rbind( summ_LW, c( subsetDF[ which(subsetDF[,subgroup]==uniq_Lw[ii])[1], group ], 
                                  uniq_Lw[ii],   sum(subsetDF[ mask,"Counts"]) ) )
  }
  colnames(summ_LW)   <- c( group, subgroup, "Counts" )
  summ_LW[,"Counts" ] <- as.numeric( summ_LW[,"Counts" ] )
  # Create the same for the upper ranks, and counting the number of lower ranks
  # present in each upper rank (e.g. how many Families there are in each Class)
  summ_UP <- data.frame()
  uniq_UP <- sort( unique( summ_LW[,group] ) )
  for ( uu in seq(1, length(uniq_UP)) ){
    mask    <- summ_LW[,group] == uniq_UP[uu]
    summ_UP <- rbind( summ_UP, c( uniq_UP[uu] , sum( mask*1 )) )
  }
  colnames(summ_UP) <- c( group, "N_subgroups" )
  summ_UP[ ,2] <- as.numeric( summ_UP[,2] )
  
  # Set tones dynamically -  the first tone color is chosen such to centered the 
  # scale around the middle of the available range, so we don't have too many dark colors
  colorScale <- NULL
  for (cc in seq(1, nrow( summ_UP )) ) {
    # According to the number of lower ranks in the cc-th upper ranks, ...
    if ( dynamicToneRange){
      if      (summ_UP[ cc, 2] <= 10)       CustPalette <- CustPalette_B  # ...use short range palette
      else if (summ_UP[ cc, 2]  > 10)       CustPalette <- CustPalette_B  # ...or use long range palette
    } else {
      CustPalette <- CustPalette_B
    }
    
    cc_group <- summ_UP[cc,1] 
    
    # Set range of tone colors. If the option of dynamic tone range is true, 
    # select first color near middle of range ...
    if ( dynamicToneRange){
      firstClr <- round(nrow(CustPalette)/2) - round(summ_UP[ cc,2] / 2 )
      if (firstClr <= 0)   {firstClr <- 1}
    } else {      # ... otherwise choose first element
      firstClr  <- 1
    }
    lastClr   <- summ_UP[ cc, 2 ]
    toneRange <- CustPalette[ firstClr:(firstClr-1+lastClr) , colnames(CustPalette)==summ_UP[ cc,1 ] ]

    
    colorScale  <- c(colorScale, toneRange ) 
  }
  
  return( list( colorScale, summ_UP, summ_LW) )
  
}
