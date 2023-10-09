#' Hill numbers diversity partitioning
#' @title Hill numbers diversity partitioning
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV MAG diversity partitioning
#' @description Compute diversity partitioning of neutral, phylogenetic or functional Hill numbers from a set of samples.
#' @param data A vector or a matrix/data.frame indicating the (relative) counts of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @param q A positive integer or decimal number (>=0), or a vector of numbers, usually between 0 and 3.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param dist A distance matrix indicating the pairwise functional distances between samples.
#' @param tau An optional maximum distance value between OTUs/ASVs/MAGs to be considered in the functional Hill numbers analyses.
#' @import tidyverse
#' @importFrom geiger tips
#' @seealso \code{\link{hilldiv}}, \code{\link{hillpair}}, \code{\link{hilldiss}}
#' @examples
#' hillpart(data=counts)
#' hillpart(data=counts,tree=tree)
#' hillpart(data=counts,dist=dist)
#' hillpart(data=counts,q=0)
#' hillpart(data=counts,q=1,tree=tree)
#' hillpart(data=counts,q=2,dist=dist)
#' @references
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.\cr\cr
# Chao et al. 2018. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89(2), e01343.\cr\cr
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export

hillpart <- function(data,q=c(0,1,2),tree,dist,tau){

  ###
  # Quality-check
  ###

  #Quality-check and warnings
  if(missing(data)) stop("Data are missing")
  if(any(q < 0) == TRUE) stop("Order(s) of diversity (q-value) need(s) to be possitive (equal or higher than zero)")

  ###
  # Hill numbers diversity partitioning functions
  ###

  hillpart.neutral <- function(data,q=c(0,1,2),tree){
    N <- ncol(data)
    pi <- tss(data)

    results <- matrix(0, nrow = length(q), ncol = 3)
    for (r in 1:length(q)){
      qvalue <- q[r]
      if(qvalue==1){
        alpha <- 1/N*exp(-sum((pi[pi!=0]/N) * log((pi[pi!=0]/N))))
        gamma <- exp(-sum(rowSums(pi/N)[rowSums(pi/N) != 0] * log(rowSums(pi/N)[rowSums(pi/N) != 0])))
        beta <- gamma / alpha
      }else{
        alpha <- 1/N*sum((pi[pi!=0]/N)^qvalue)^(1/(1-qvalue))
        gamma <- sum(rowSums(pi/N)^qvalue)^(1/(1-qvalue))
        beta <- gamma / alpha
        }
        results[r, 1] <- alpha
        results[r, 2] <- gamma
        results[r, 3] <- beta
      }
    rownames(results) <- paste0("q",q)
    colnames(results) <- c("alpha","gamma","beta")
    return(results)

  }

  hillpart.phylogenetic <- function(data,q=c(0,1,2),tree){
    N <- ncol(data)
    data <- tss(data)
    Li <- tree$edge.length
    ltips <- sapply(tree$edge[, 2], function(node) geiger::tips(tree, node))
    aij <- matrix(unlist(lapply(ltips, function(TipVector) colSums(data[TipVector,]))), ncol = N, byrow = TRUE)
    ai <- rowSums(aij)
    T <- sum(sweep(aij, 1, Li, "*"))
    L <- matrix(rep(Li, N), ncol = N)
    Li <- Li[ai != 0]
    ai <- ai[ai != 0]
    i <-  which(aij > 0)
    alpha <- 1/N*exp(-sum(L[i] * aij[i]/T * log(aij[i]/T)))
    results <- matrix(0, nrow = length(q), ncol = 3)
    for (r in 1:length(q)){
      qvalue <- q[r]
      if(qvalue==1){
        alpha <- 1/N*exp(-sum(L[i] * aij[i]/T * log(aij[i]/T)))
        gamma <- exp(-sum(Li * (ai/T) * log(ai/T)))
        beta <- gamma / alpha
      }else{
        alpha <- 1/N*sum(L[i] * (aij[i]/T)^qvalue)^(1/(1 - qvalue))
        gamma <- (sum(Li * (ai/T)^qvalue)^(1/(1 - qvalue)))
        beta <- gamma / alpha
        }
        results[r, 1] <- alpha
        results[r, 2] <- gamma
        results[r, 3] <- beta
      }
    rownames(results) <- paste0("q",q)
    colnames(results) <- c("alpha","gamma","beta")
    return(results)

  }

  hillpart.functional <- function(data,q=c(0,1,2),dist,tau){
    if(missing(tau)){tau=max(dist)} #if Tau is not declared, use maximum distance
    N <- ncol(data)
    dij <- as.matrix(dist)
    dij[which(dij>tau,arr.ind = T)] <- tau
    aik <- apply(data, 2, function(col) as.vector((1 - dij/tau) %*% col)) # aik
    aiplus <- apply(aik, 1, sum)
    vi <- apply(data, 1, sum)/aiplus
    alpha_v <- rep(vi, N)
    nplus <- sum(data)
    aik <- as.vector(aik)
    alpha_v <- alpha_v[aik!=0]
    aik <- aik[aik!=0]
    results <- matrix(0, nrow = length(q), ncol = 3)
    for (r in 1:length(q)){
      qvalue <- q[r]
      if(qvalue==1){
        alpha <- 1/N*exp(sum(-alpha_v*aik/nplus*log(aik/nplus)))
        gamma <- exp(sum(-vi*aiplus/nplus*log(aiplus/nplus)))
        beta <- gamma / alpha
      }else{
        alpha <- 1/N*(sum(alpha_v*(aik/nplus)^qvalue))^(1 / (1-qvalue))
        gamma <- (sum(vi*(aiplus/nplus)^qvalue))^(1 / (1-qvalue))
        beta <- gamma / alpha
        }
      results[r, 1] <- alpha
      results[r, 2] <- gamma
      results[r, 3] <- beta
      }
    rownames(results) <- paste0("q",q)
    colnames(results) <- c("alpha","gamma","beta")
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

  ###
  # Run calculations
  ###

  # Neutral Hill numbers
  if(hilltype == "neutral"){
      qvalues <- paste(paste0("q",q),collapse = ", ")
      message(paste0("Partitioning of neutral Hill numbers of ", qvalues))
      return(hillpart.neutral(data=data,q=q))
  }

  # Phylogenetic Hill numbers
  if(hilltype == "phylogenetic"){
      if(is.null(names(data)) == TRUE) stop("Error: The vector needs to contain names in order to link it to the tree tips")
      if(all(sort(rownames(data)) != sort(tree$tip.label))) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
      qvalues <- paste(paste0("q",q),collapse = ", ")
      message(paste0("Partitioning of phylogenetic Hill numbers of ", qvalues))
      return(hillpart.phylogenetic(data=data,q=q,tree=tree))
  }

  # Functional Hill numbers
  if(hilltype == "functional"){
      if(all(sort(rownames(data)) != sort(rownames(dist)))) stop("The count table and the functional distance table do not match.")
      qvalues <- paste(paste0("q",q),collapse = ", ")
      message(paste0("Partitioning of functional Hill numbers of ", qvalues))
      return(hillpart.functional(data=data,q=q,dist=dist))
  }

}
