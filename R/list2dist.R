#' Convert a long table of pairwise distances to a distance matrix
#' @title List to dist
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV MAG diversity
#' @description Convert a long table of pairwise distances to a distance matrix
#' @param distlist A table containing pairwise distances between samples, with sample names in the first 2 columns and the distance in the third column,
#' @examples
#' list2dist(distlist)
#' @export

#Distance matrix function
list2dist <- function(distlist){
    distlist <- distlist[order(distlist[,1],distlist[,2]),]
    distlist.name1 <- as.character(distlist[, 1])
    distlist.name2 <- as.character(distlist[, 2])
    distlist.value <- distlist[, 3]
    names1 <- sort(unique(as.character(distlist[, 1])))
    names2 <- sort(unique(as.character(distlist[, 2])))
    total.names <- unique(c(names1, names2))
    elements <- rep(NA, length(total.names)^2)
    dim(elements) <- c(length(total.names), length(total.names))
    rownames(elements) <- total.names
    colnames(elements) <- total.names
    for (i in 1:length(total.names)) {
        for (j in 1:length(total.names)) {
            for (k in 1:length(distlist.name1)) {
                if ((total.names[i] == distlist.name1[k]) & (total.names[j] ==
                  distlist.name2[k])) {
                  elements[i, j] <- distlist.value[k]
                }
            }
        }
    }
    res <- as.dist(t(elements))
    return(res)
}
