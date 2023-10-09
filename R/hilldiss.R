#' Hill numbers-based dissimilarity computation
#' @title Hill number dissimilarity computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV MAG diversity dissimilarity overlap similarity
#' @description Compute dissimilarity metrics based on neutral, phylogenetic and/or functional Hill numbers.
#' @param data A matrix/data.frame indicating the (relative) OTU/ASV/MAG counts of multiple samples. Columns must refer to samples and rows to OTUs/ASVs/MAGs.
#' @param q A positive integer or decimal number (>=0), or a vector of numbers, usually between 0 and 3.
#' @param metric Dissimilarity metric(s) to be computed. Default is all four: c("S","C","U","V")
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param dist A distance matrix indicating the pairwise functional distances between samples.
#' @param tau An optional maximum distance value between OTUs/ASVs/MAGs to be considered in the functional Hill numbers analyses.
#' @import tidyverse
#' @seealso \code{\link{hilldiv}}, \code{\link{hillpair}}, \code{\link{hillpart}}
#' @examples
#' hilldiss(data=counts)
#' hilldiss(data=counts,tree=tree)
#' hilldiss(data=counts,dist=dist)
#' hilldiss(data=counts,q=0)
#' hilldiss(data=counts,q=1,tree=tree)
#' hilldiss(data=counts,q=2,dist=dist)
#' @references
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.\cr\cr
# Chao et al. 2018. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89(2), e01343.\cr\cr
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export

hilldiss <- function(data,q=c(0,1,2),metric=c("S","C","U","V"),tree,dist,tau){

  ###
  # Quality-check
  ###

  #Quality-check and warnings
  if(missing(data)) stop("Data are missing")
  if(any(q < 0) == TRUE) stop("Order(s) of diversity (q-value) need(s) to be possitive (equal or higher than zero)")

  ###
  # Type of analysis detection
  ###

  #Hill numbers type definition
  if(missing(tree) && missing(dist)){hilltype="neutral"}
  if(!missing(tree) && missing(dist)){hilltype="phylogenetic"}
  if(missing(tree) && !missing(dist)){hilltype="functional"}
  if(!missing(tree) && !missing(dist)) stop("Phylogenetic and functional trait information cannot be added at once. Use either phylogenetic (tree) or functional (fun) information.")

  ###
  # Dissimilarity functions
  ###

  #Homogeneity-complement (Jaccard-type)
  hilldiss.S <- function(beta,N){
    1-(((1/beta) - 1/N)/(1 - 1/N)) %>%
    return()
  }

  #Local overlap-complement (Sørensen type)
  hilldiss.C <- function(beta,N,q){
    if(qvalue==1){
      1-(((1/beta)^(0.999999 - 1) - (1/N)^(0.999999 - 1))/(1 - (1/N)^(0.999999 - 1))) %>%
      return()

    }else{
      1-(((1/beta)^(qvalue - 1) - (1/N)^(qvalue - 1))/(1 - (1/N)^(qvalue - 1))) %>%
      return()
    }
  }

  #Turnover (Sørensen type)
  hilldiss.V <- function(beta,N){
    1-((N - beta)/(N - 1)) %>%
    return()
  }

  #Regional overlap-complement (Jaccard-type)
  hilldiss.U <- function(beta,N,q){
    if(qvalue==1){
      1-(((1/beta)^(1 - 0.999999) - (1/N)^(1 - 0.999999))/(1 - (1/N)^(1 - 0.999999))) %>%
      return()
    }else{
      1-(((1/beta)^(1 - qvalue) - (1/N)^(1 - qvalue))/(1 - (1/N)^(1 - qvalue))) %>%
      return()
    }
  }

  ###
  # Run calculations
  ###

    N <- ncol(data)

    # Neutral Hill numbers
    results <- matrix(0, nrow = length(q), ncol = 4)
    if(hilltype == "neutral"){
      qvalues <- paste(paste0("q",q),collapse = ", ")
      message(paste0("Dissimilarities based on neutral Hill numbers of ", qvalues))
      suppressMessages(betas <- hillpart(data=data,q=q)[,"beta"]) # calculate neutral beta value(s)
      }
    if(hilltype == "phylogenetic"){
      qvalues <- paste(paste0("q",q),collapse = ", ")
      message(paste0("Dissimilarities based on phylogenetic Hill numbers of ", qvalues))
      suppressMessages(betas <- hillpart(data=data,q=q,tree=tree)[,"beta"]) # calculate phylogenetic beta values(s)
      }
    if(hilltype == "functional"){
      qvalues <- paste(paste0("q",q),collapse = ", ")
      message(paste0("Dissimilarities based on functional Hill numbers of ", qvalues))
      suppressMessages(betas <- hillpart(data=data,q=q,dist=dist,tau=tau)[,"beta"]) # calculate functional beta values(s)
    }
      for(i in c(1:length(q))){
          beta <- as.numeric(betas[i])
          qvalue <- as.numeric(q[i])
          results[i,1] <- hilldiss.S(beta,N)
          results[i,2] <- hilldiss.C(beta,N,qvalue)
          results[i,3] <- hilldiss.V(beta,N)
          results[i,4] <- hilldiss.U(beta,N,qvalue)
      }

    rownames(results) <- paste0("q",q)
    colnames(results) <- c("S","C","V","U")
    return(results[,metric])

}
