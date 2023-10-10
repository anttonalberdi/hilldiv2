#' Hill numbers redundancy
#' @title Hill numbers redundancy computation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords phylogenetic functional redundancy
#' @description Compute redundancy of phylogenetic or functional Hill numbers for a set of samples.
#' @param data A matrix/data.frame indicating the (relative) counts of multiple samples. Columns must refer to samples and rows to OTUs/ASVs/MAGs.
#' @param q A positive integer or decimal number (>=0), or a vector of numbers, usually between 0 and 3.
#' @param tree A phylogenetic tree of class 'phylo'. The tip labels must match the names of the vector values (if one sample) or matrix rows (if multiple samples). Use the function match_data() if the OTU names do not match.
#' @param dist A distance matrix indicating the pairwise functional distances between samples.
#' @param tau An optional maximum distance value between OTUs/ASVs/MAGs to be considered in the functional Hill numbers analyses.
#' @import tidyverse
#' @importFrom geiger tips
#' @seealso \code{\link{hilldiss}}, \code{\link{hillpair}}, \code{\link{hillpart}}, \code{\link{hilldiv}}
#' @examples
#' hillredun(data=counts,tree=tree)
#' hillredun(data=counts,dist=dist)
#' hillredun(data=counts,q=1,tree=tree)
#' hillredun(data=counts,q=2,dist=dist)
#' @references
#' Jost, L. (2006). Entropy and diversity. Oikos, 113, 363-375.\cr\cr
#' Hill, M. O. (1973). Diversity and evenness: a unifying notation and its consequences. Ecology, 54, 427-432.\cr\cr
# Chao et al. 2018. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89(2), e01343.\cr\cr
#' Alberdi, A., Gilbert, M.T.P. (2019). A guide to the application of Hill numbers to DNA-based diversity analyses. Molecular Ecology Resources, 19, 804-817.
#' @export

hillredun <- function(data,q=c(0,1,2),tree,dist,tau){

    ###
    # Quality-check
    ###

    #Quality-check and warnings
    if(missing(data)) stop("Data are missing")
    if(any(q < 0) == TRUE) stop("Order(s) of diversity (q-value) need(s) to be possitive (equal or higher than zero)")

    ###
    # Definition of the neutral-functional relationship function
    ###

    # We use an increasing form of an exponential decat function where:
    # a = range of the functional Hill numbers (y value) between the lowest value in the set samples or population and the saturation value
    # b = units of neutral Hill numbers required to increase half of the range of functional Hill numbers
    # c = value of the saturation point of functional Hill numbers (y value)

    relationship_function <- function(x, a, b, c) {
      return(-a*2^(-x/b)+c)
    }

    # As the model includes is writen in terms of powers of 2 instead of powers of exp, the b parameter is interpretable
    # as the units of neutral Hill numbers required to increase half of the range of functional Hill numbers

    ###
    # Type of analysis detection
    ###

    #Hill numbers type definition
    if(missing(tree) && missing(dist)) stop("Hill numbers redundancy analysis requires either phylogenetic or functional information to be inputed.")
    if(!missing(tree) && missing(dist)){hilltype="phylogenetic"}
    if(missing(tree) && !missing(dist)){hilltype="functional"}
    if(!missing(tree) && !missing(dist)) stop("Phylogenetic and functional trait information cannot be added at once. Use either phylogenetic (tree) or functional (fun) information.")

    ###
    # Run calculations
    ###

    # Phylogenetic Hill numbers
    if(hilltype == "phylogenetic"){
        if(nrow(data) != length(tree$tip.label)) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
        if(all(sort(rownames(data)) != sort(tree$tip.label))) stop("The OTU/ASV/MAG names in the count table and tree tips do not match.")
        qvalues <- paste(paste0("q",q),collapse = ", ")

        x=c(hilldiv(data=data,q=q))
        y=c(hilldiv(data=data,q=q,tree=tree))

    }

    # Functional Hill numbers
    if(hilltype == "functional"){
        if(nrow(data) != nrow(dist)) stop("The count table and the functional distance table do not match.")
        if(all(sort(rownames(data)) != sort(rownames(dist)))) stop("The count table and the functional distance table do not match.")
        dist <- dist[rownames(data),rownames(data)]
        qvalues <- paste(paste0("q",q),collapse = ", ")

        x=c(hilldiv(data=data,q=q))
        y=c(hilldiv(data=data,q=q,dist=dist))

    }


}
