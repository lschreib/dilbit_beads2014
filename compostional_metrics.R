library("propr")
library("EnvStats")

# The presented functions calculate compositional metrics based on OTU/ASV count tables with
# rows representing OTUs and columns representing samples

propr_modified.alr_to_phi <- function(ct,ivar){
  # This function performes a alr transformation of the count table and subsequently
  # calculates the phi dissimilarity between samples
  # The phi dissimilarities are returned in the form of a distance matrix

  # The format of ct should be samples in columns and OTUs in rows

  # The row number of the reference OTU for the alr transformation is passed to the
  # function using the "ivar" parameter 

  # Find the minimum of the complete count table and replace all zeros with that
  zeros <- ct == 0
  ct[zeros] <- min(ct[!zeros])
  
  # Calculate the alr (log-ratio) table
  lr <- apply(ct, MARGIN = 2, function(x) log(x/x[ivar]))
  
  # Calculate the phi metric using the "propr" package
  mat <- propr:::lr2phi(as.matrix(lr))
  #Fix the row and colnames
  colnames(mat) <- colnames(lr)
  rownames(mat) <- colnames(lr)
  
  mat.dist <- as.dist(mat, upper = TRUE, diag = TRUE)
  return(mat.dist)
}

propr_modified.ct_to_alr <- function(ct,ivar){
  # This function performes a alr transformation of a OTU/ASV count table

  # The format of ct should be OTUs in rows and samples in columns

  # The row number of the reference OTU for the alr transformation is passed to the
  # function using the "ivar" parameter 

  # Find the minimum of the complete count table and replace all zeros with that
  zeros <- ct == 0
  ct[zeros] <- min(ct[!zeros])
  
  #Calculate the alr (log-ratio) table
  lr <- apply(ct, MARGIN = 2, function(x) log(x/x[ivar]))
  
  return(lr)
}

propr_modified.ct_to_clr <- function(ct){
  # This function performes a clr transformation of a OTU/ASV count table

  # Find the minimum of the complete count table and replace all zeros with that
  zeros <- ct == 0
  ct[zeros] <- min(ct[!zeros])
  
  #Calculate the alr (log-ratio) table
  #Calculate the clr log-ratio table
  lr <- apply(ct, MARGIN = 2, function(x) log(x/geoMean(x)))
  
  return(lr)
}

propr_modified.clr_to_phi <- function(ct){
  # This function performes a clr transformation of a OTU/ASV count table and subsequently
  # calculates the phi dissimilarity between samples
  # The phi dissimilarities between samples are returned in the form of a distance matrix

  # The format of ct should be samples in columns and OTUs in rows

  # Find the minimum of the complete count table and replace all zeros with that
  zeros <- ct == 0
  ct[zeros] <- min(ct[!zeros])
    
  #Calculate the clr log-ratio table
  lr <- apply(ct, MARGIN = 2, function(x) log(x/geoMean(x)))
  
  #Calculate the phi metric
  mat <- propr:::lr2phi(as.matrix(lr))
  #Fix the row and colnames
  colnames(mat) <- colnames(lr)
  rownames(mat) <- colnames(lr)
  
  mat.dist <- as.dist(mat, upper = TRUE, diag = TRUE)
  return(mat.dist)
}