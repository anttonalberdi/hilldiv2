#' Convert traits matrix into a distance matrix
#' @title Traits to distance
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords OTU ASV MAG diversity
#' @description Convert traits matrix into a distance matrix.
#' @param traits A table containing trait information for each OTU/ASV/MAG in columns and OTUs/ASVs/MAGs in rows.
#' @param method The distance metric to be applied: "euclidean", "manhattan" or "gower" (default: "gower").
#' @import tidyverse
#' @importFrom cluster daisy
#' @examples
#' traits2dist(trait,method="gower")
#' @export

#Distance matrix function
traits2dist <- function(traits,method="gower"){
  traits %>%
  as.data.frame() %>%
  select_if(~ !all(. == .[1])) %>%
  cluster::daisy(., metric = method, warnType = F) %>% #metric could be "euclidean", "manhattan" or "gower".
      as.matrix() %>%
      return()
}
