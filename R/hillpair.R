#' Pairwise Hill numbers-based dissimilarity computation
#' @title Pairwise Hill number dissimilarity computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV MAG diversity dissimilarity overlap similarity
#' @description Compute pairwise dissimilarity metrics based on neutral, phylogenetic and/or functional Hill numbers.
#' @param data A matrix/data.frame indicating the (relative) OTU/ASV/MAG counts of multiple samples. Columns must refer to samples and rows to OTUs/ASVs/MAGs.
#' @param q A positive integer or decimal number (>=0), or a vector of numbers, usually between 0 and 3.
#' @param metric Dissimilarity metric(s) to be computed. Default is all four: c("S","C","U","V")
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param dist A distance matrix indicating the pairwise functional distances between samples.
#' @param tau An optional maximum distance value between OTUs/ASVs/MAGs to be considered in the functional Hill numbers analyses.
#' @param out Type of output object. "dist" for a distance object of pairwise distances (default), "pair" for a table of pairwise distances in rows.
#' @import tidyverse
#' @seealso \code{\link{hilldiv}}, \code{\link{hilldiss}}, \code{\link{hillpart}}
#' @examples
#' hillpair(data=counts)
#' hillpair(data=counts,tree=tree)
#' hillpair(data=counts,dist=dist)
#' hillpair(data=counts,q=0)
#' hillpair(data=counts,q=1,tree=tree)
#' hillpair(data=counts,q=2,dist=dist)
#' @references
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.\cr\cr
# Chao et al. 2018. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89(2), e01343.\cr\cr
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export


hillpair <- function(data,q=c(0,1,2),metric=c("S","C","U","V"),tree,dist,tau,out="dist"){

  ###
  # Type of analysis detection
  ###

  #Hill numbers type definition
  if(missing(tree) && missing(dist)){hilltype="neutral"}
  if(!missing(tree) && missing(dist)){hilltype="phylogenetic"}
  if(missing(tree) && !missing(dist)){hilltype="functional"}
  if(!missing(tree) && !missing(dist)) stop("Phylogenetic and functional trait information cannot be added at once. Use either phylogenetic (tree) or functional (fun) information.")

  ###
  # Create pairwise tables
  ###

  # List of pairwise combinations
  pairs <- list()
  for (i in 1:(ncol(data)-1)) {
    for (j in (i + 1):ncol(data)) {
      pair <- data[, c(i, j), drop = FALSE]
      pairs <- append(pairs, list(pair))
    }
  }

  # List of pairwise sample names
  pairs_names <- list()
  for (i in 1:(ncol(data)-1)) {
    for (j in (i + 1):ncol(data)) {
      pair_names <- colnames(data[, c(i, j), drop = FALSE])
      pairs_names <- append(pairs_names, list(pair_names))
    }
  }

  ###
  # Estimate processing time
  ###
  pairs_length <- length(pairs_names)
  start.time <- Sys.time()
  if(hilltype == "neutral"){suppressMessages(pairs_dist <- lapply(pairs[1:10], hilldiss, q=q, metric=metric))}
  if(hilltype == "phylogenetic"){suppressMessages(pairs_dist <- lapply(pairs[1:10], hilldiss, q=q, metric=metric, tree=tree))}
  if(hilltype == "functional"){suppressMessages(pairs_dist <- lapply(pairs[1:10], hilldiss, q=q, metric=metric, dist=dist, tau=tau))}
  end.time <- Sys.time()
  time_estimate <- round((end.time - start.time)/10*pairs_length)
  options(warn=1)
  message(paste0("The estimated time for calculating the ",pairs_length," pairwise combinations is ",time_estimate+2, " seconds."))

  ###
  # Compute all pairwise distances
  ###

   if(hilltype == "neutral"){suppressMessages(pairs_dist <- lapply(pairs, hilldiss, q=q, metric=metric))}
   if(hilltype == "phylogenetic"){suppressMessages(pairs_dist <- lapply(pairs, hilldiss, q=q, metric=metric, tree=tree))}
   if(hilltype == "functional"){suppressMessages(pairs_dist <- lapply(pairs, hilldiss, q=q, metric=metric, dist=dist, tau=tau))}

  ###
  # Arrange output
  ###

  #Flatten list into matrix
  flatten_matrix <- function(matrix){
      matrix %>%
      as.data.frame() %>%
      pivot_longer(cols = everything(), names_to = "Variable") %>%
      select(value) %>%
      pull()
      }

    pairs_dist_matrix <- pairs_dist %>%
      lapply(flatten_matrix) %>%
      do.call(rbind,.) %>%
      as.data.frame()

    #List metric combinations
    if(length(q)>1){
      metric_combinations <- expand.grid(colnames(pairs_dist[[1]]),rownames(pairs_dist[[1]])) %>%
          mutate(name = paste(Var2, Var1, sep = "") ) %>%
          select(name)  %>%
          pull()
    }else{
      metric_combinations <- metric
    }

    #Assign column names to different metrics
    colnames(pairs_dist_matrix) <- metric_combinations

    #Add sample names
    pairs_dist_matrix <- pairs_dist_matrix %>%
      mutate(first=unlist(lapply(pairs_names, function(x) x[1]))) %>%
      relocate(first, .before = 1) %>%
      mutate(second=unlist(lapply(pairs_names, function(x) x[2]))) %>%
      relocate(second, .after = 1)

    #Convert table into distance object
    if (out == "dist"){
          create_distance_matrix <- function(column_name) {
            matrix_df <- pairs_dist_matrix %>%
              select(first, second, {{column_name}}) %>%
              pivot_wider(names_from = first, values_from = {{column_name}}) %>%
              as.data.frame() %>%
              column_to_rownames(., var = "second") %>%
              add_row(., .before = 1)
            rownames(matrix_df)[1] <- pairs_dist_matrix[1,1]
            suppressWarnings(matrix_df <- as.dist(matrix_df))
            return(matrix_df)
          }

      pairs_dist_matrix <- lapply(colnames(pairs_dist_matrix[3:ncol(pairs_dist_matrix)]), create_distance_matrix)
      names(pairs_dist_matrix) <- metric_combinations

      if(length(pairs_dist_matrix) == 1){
        pairs_dist_matrix <- pairs_dist_matrix[[1]]
      }
    }

  return(pairs_dist_matrix)

}
