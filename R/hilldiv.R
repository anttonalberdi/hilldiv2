#' Hill numbers computation
#' @title Hill numbers computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV MAG diversity
#' @description Compute neutral, phylogenetic and/or functional Hill numbers from a single sample (vector) or count table (matrix).
#' @param data A vector or a matrix/data.frame indicating the (relative) counts of one or multiple samples, respectively. If a matrix/data.frame is provided, must refer to samples and rows to OTUs/ASVs/MAGs.
#' @param q A positive integer or decimal number (>=0), or a vector of numbers, usually between 0 and 3.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param dist A distance matrix indicating the pairwise functional distances between samples.
#' @param tau An optional maximum distance value between OTUs/ASVs/MAGs to be considered in the functional Hill numbers analyses.
#' @import tidyverse
#' @importFrom geiger tips
#' @seealso \code{\link{hilldiss}}, \code{\link{hillpair}}, \code{\link{hillpart}}
#' @examples
#' hilldiv(data=counts)
#' hilldiv(data=counts,tree=tree)
#' hilldiv(data=counts,dist=dist)
#' hilldiv(data=counts,q=0)
#' hilldiv(data=counts,q=1,tree=tree)
#' hilldiv(data=counts,q=2,dist=dist)
#' @references
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.\cr\cr
# Chao et al. 2018. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89(2), e01343.\cr\cr
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export

hilldiv <- function(data,q=c(0,1,2),tree,dist,tau){

  ###
  # Quality-check
  ###

  #Quality-check and warnings
  if(missing(data)) stop("Data are missing")
  if(any(q < 0) == TRUE) stop("Order(s) of diversity (q-value) need(s) to be possitive (equal or higher than zero)")

  ###
  # Hill numbers functions
  ###

  #Neutral Hill numbers: Basal function for a single sample
  hilldiv.neutral <- function(vector,q=c(0,1,2)){
      pi <- tss(vector) # transform to 0-1
      results <- rep(0,length(q))
      for (i in 1:length(q)){
        qvalue <- q[i]
        if(qvalue==1){
            results[i] <- exp(-sum(pi[pi!=0] * log(pi[pi!=0]))) # special case for q=1
        }else{
            results[i] <- sum(pi[pi!=0]^qvalue)^(1/(1-qvalue)) # all other q values
            }
          }
        names(results) <- paste0("q",q)
        return(results)
      }

  #Neutral Hill numbers: Function for multiple samples and q values
  hilldiv.neutral.multi <- function(matrix,q=c(0,1,2)){
      pi.multi <- tss(matrix) # transform to 0-1
      results <- matrix(0, nrow = length(q), ncol = ncol(matrix))
      for (i in 1:length(q)){
        qvalue <- q[i]
        if(qvalue==1){
            results[i, ] <- apply(pi.multi, 2, function(pi) exp(-sum(pi[pi!=0] * log(pi[pi!=0])))) # special case for q=1
        }else{
            results[i, ] <- apply(pi.multi, 2, function(pi) sum(pi[pi!=0]^qvalue)^(1/(1-qvalue))) # all other q values
            }
      }
      rownames(results) <- paste0("q",q)
      colnames(results) <- colnames(matrix)
      return(results)
      }

  #Phylogentic Hill numbers: Basal function for a single sample
  hilldiv.phylogenetic <- function(vector,q=c(0,1,2),tree){
      Li <- tree$edge.length
      ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
      ai <- unlist(lapply(ltips, function(TipVector) sum(tss(vector)[TipVector])))
      T <- sum(Li * ai)
      results <- rep(0,length(q))
      for (i in 1:length(q)){
        qvalue <- q[i]
        if(qvalue==1){
            results[i] <- exp(-sum(Li[ai != 0] * (ai[ai != 0]/T) * log(ai[ai != 0]/T)))/T # special case for q=1
        }else{
            results[i] <- sum(Li[ai != 0] * (ai[ai != 0]/T)^qvalue)^(1/(1-qvalue))/T # all other q values
            }
        }
      names(results) <- paste0("q",q)
      return(results)
      }

  #Phylogentic Hill numbers: Function for multiple samples and q values
  #When averaging T value, it can yield imperfections in lowest diversity samples if there is a large difference with others
  hilldiv.phylogenetic.multi <- function(matrix,q=c(0,1,2),tree){
      Li <- tree$edge.length
      ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
      ai.multi <- apply(matrix, 2, function(vector) unlist(lapply(ltips, function(TipVector) sum(tss(vector)[TipVector]))))

      #Calculate reference T for an even distribution of all present MAGs
      pool <- matrix %>%
          rowwise() %>%
          mutate(pool = if_else(any(c_across(everything()) != 0), 1, 0)) %>%
          select(pool) %>%
          pull() #generate a pool of equally abundant MAGs
      names(pool) <- rownames(matrix)
      ai.all <- unlist(lapply(ltips, function(TipVector) sum(tss(pool)[TipVector])))
      T <- sum(Li * ai.all)

      #Calculate present OTUs/ASVs/MAGs
      present <- colSums(matrix != 0)

      #Calculate Hill numbers
      results <- matrix(0, nrow = length(q), ncol = ncol(matrix))
      for (i in 1:length(q)){
        qvalue <- q[i]
        if(qvalue==1){
            phylo.q1.func <- function(ai, p) {
              if(p == 0){
                return(0)
              }else if(p == 1){
                return(1)
              }else{
                return(exp(-sum(Li[ai!=0] * (ai[ai!=0]/T) * log(ai[ai!=0]/T)))/T)
              }
            }
            results[i, ] <- sapply(1:ncol(ai.multi), function(i, p) phylo.q1.func(ai.multi[, i], present[i])) # special case for q=1
        }else{
            phylo.all.func <- function(ai, p) {
              if(p == 0){
                return(0)
              }else if(p == 1){
                return(1)
              }else{
               return(sum(Li[ai!=0] * (ai[ai!=0]/T)^qvalue)^(1/(1-qvalue))/T)
              }
            }
            results[i, ] <- sapply(1:ncol(ai.multi), function(i, p) phylo.all.func(ai.multi[, i], present[i])) # all other q values
            }
        }

      #Add column and row names to results table
      rownames(results) <- paste0("q",q)
      colnames(results) <- colnames(matrix)
      return(results)
      }

  #Phylogentic Hill numbers: Function for multiple samples and q values
  #Function to calculate Hill numbers with a different T per sample (problematic for diversity partitioning)
  hilldiv.phylogenetic.multi_independentT <- function(matrix,q=c(0,1,2),tree){
      Li <- tree$edge.length
      ltips <- sapply(tree$edge[, 2], function(node) tips(tree, node))
      ai.multi <- apply(matrix, 2, function(vector) unlist(lapply(ltips, function(TipVector) sum(tss(vector)[TipVector]))))
      T.multi <- apply(ai.multi, 2, function(ai) sum(Li * ai))
      results <- matrix(0, nrow = length(q), ncol = ncol(matrix))
      for (i in 1:length(q)){
        qvalue <- q[i]
        if(qvalue==1){
            phylo.q1.func <- function(ai,T) {exp(-sum(Li[ai!=0] * (ai[ai!=0]/T) * log(ai[ai!=0]/T)))/T}
            results[i, ] <- sapply(1:ncol(ai.multi), function(i) phylo.q1.func(ai.multi[, i], T.multi[i])) # special case for q=1
        }else{
            phylo.all.func <- function(ai,T) {sum(Li[ai!=0] * (ai[ai!=0]/T)^qvalue)^(1/(1-qvalue))/T}
            results[i, ] <- sapply(1:ncol(ai.multi), function(i) phylo.all.func(ai.multi[, i], T.multi[i])) # all other q values
            }
        }
      rownames(results) <- paste0("q",q)
      colnames(results) <- colnames(matrix)
      return(results)
      }

  #Functional Hill numbers: Basal function for a single sample
  hilldiv.functional <- function(vector, q=c(0,1,2), dist, tau){
     if(missing(tau)){tau=max(dist)} #if Tau is not declared, use maximum distance
     dij <- as.matrix(dist)
     dij[which(dij>tau,arr.ind = T)] <- tau
     vector <- tss(vector)
     a <- as.vector((1 - dij/tau) %*% vector)
     vector <- vector[a!=0]
     a <- a[a!=0]
     v <- vector/a
     results <- rep(0,length(q))
     for (i in 1:length(q)){
       qvalue <- q[i]
       if(qvalue==1){
          results[i] <- exp(sum(-v*a*log(a))) # special case for q=1
       }else{
          results[i] <- (sum(v*a^qvalue))^(1 / (1-qvalue)) # all other q values
          }
       }
       names(results) <- paste0("q",q)
       return(results)
     }
  #Functional Hill numbers: Function for multiple samples
  hilldiv.functional.multi <- function(matrix, q=c(0,1,2), dist, tau){
      if(missing(tau)){tau=max(dist)} #if Tau is not declared, use maximum distance
      dij <- as.matrix(dist)
      dij[which(dij>tau,arr.ind = T)] <- tau
      matrix <- tss(matrix)
      a.multi <- apply(matrix, 2, function(col) as.vector((1 - dij/tau) %*% col))
      a.multi <- a.multi[rowSums(a.multi) != 0,]
      matrix <- matrix[rowSums(matrix) != 0,]
      v.multi <- matrix/a.multi
      results <- matrix(0, nrow = length(q), ncol = ncol(matrix))
      for (i in 1:length(q)){
        qvalue <- q[i]
        if(qvalue==1){
           func.q1.func <- function(a,v) {exp(sum(-v*a*log(a)))}
           results[i, ] <- sapply(1:ncol(a.multi), function(i) func.q1.func(a.multi[, i], v.multi[,i])) # special case for q=1
        }else{
           func.all.func <- function(a,v) {(sum(v*a^qvalue))^(1 / (1-qvalue))}
           results[i, ] <- sapply(1:ncol(a.multi), function(i) func.all.func(a.multi[, i], v.multi[,i])) # all other q values
           }
         }
      rownames(results) <- paste0("q",q)
      colnames(results) <- colnames(matrix)
      return(results)
      }

    ###
    # Type of analysis detection
    ###

    #Hill numbers type definition
    if(missing(tree) && missing(dist)){hilltype="neutral"}
    if(!missing(tree) && missing(dist)){hilltype="phylogenetic"}
    if(missing(tree) && !missing(dist)){hilltype="functional"}
    if(!missing(tree) && !missing(dist)) stop("Phylogenetic and functional trait information cannot be added at once. Use either phylogenetic (tree) or functional (fun) information.")

    #Input data type definition
    if(is.null(dim(data)) == TRUE){
    datatype="one"
    }else{
    datatype="multi"
    }

    ###
    # Run calculations
    ###

    # Neutral Hill numbers
    if(hilltype == "neutral" & datatype == "one"){
        qvalues <- paste(paste0("q",q),collapse = ", ")
        message(paste0("Neutral Hill numbers of ", qvalues))
        return(hilldiv.neutral(vector=data,q=q))
    }
    if(hilltype == "neutral" & datatype == "multi"){
        qvalues <- paste(paste0("q",q),collapse = ", ")
        message(paste0("Neutral Hill numbers of ", qvalues))
        return(hilldiv.neutral.multi(matrix=data,q=q))
    }

    # Phylogenetic Hill numbers
    if(hilltype == "phylogenetic" & datatype == "one"){
        if(is.null(names(data)) == TRUE) stop("Error: The vector needs to contain names in order to link it to the tree tips")
        if(length(data) != length(tree$tip.label)) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
        if(all(sort(names(data)) != sort(tree$tip.label))) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
        qvalues <- paste(paste0("q",q),collapse = ", ")
        message(paste0("Phylogenetic Hill numbers of ", qvalues))
        return(hilldiv.phylogenetic(vector=data,q=q,tree=tree))
    }
    if(hilltype == "phylogenetic" & datatype == "multi"){
        if(nrow(data) != length(tree$tip.label)) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
        if(all(sort(rownames(data)) != sort(tree$tip.label))) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
        qvalues <- paste(paste0("q",q),collapse = ", ")
        message(paste0("Phylogenetic Hill numbers of ", qvalues))
        return(hilldiv.phylogenetic.multi(matrix=data,q=q,tree=tree))
    }

    # Functional Hill numbers
    if(hilltype == "functional" & datatype == "one"){
        if(length(data) != nrow(dist)) stop("The count vector and the functional distance table do not match.")
        if(all(sort(names(data)) != sort(rownames(dist)))) stop("The count vector and the functional distance table do not match.")
        dist <- dist[rownames(data),rownames(data)]
        qvalues <- paste(paste0("q",q),collapse = ", ")
        message(paste0("Functional Hill numbers of ", qvalues))
        return(hilldiv.functional(vector=data,q=q,dist=dist,tau=tau))
    }
    if(hilltype == "functional" & datatype == "multi"){
        if(nrow(data) != nrow(dist)) stop("The count table and the functional distance table do not match.")
        if(all(sort(rownames(data)) != sort(rownames(dist)))) stop("The count table and the functional distance table do not match.")
        dist <- dist[rownames(data),rownames(data)]
        qvalues <- paste(paste0("q",q),collapse = ", ")
        message(paste0("Functional Hill numbers of ", qvalues))
       return(hilldiv.functional.multi(matrix=data,q=q,dist=dist,tau=tau))
    }
}
