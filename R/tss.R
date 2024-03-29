#' Total Sum Scaling normalisation
#' @title Total Sum Scaling normalisation
#' @author Antton Alberdi, \email{anttonalberdi@gmail.com}
#' @keywords normalisation Hill
#' @description Normalise a vector or count matrix to the range of 0-1.
#' @param abund A vector or a matrix/data.frame indicating the relative abundances of one or multiple samples, respectively. If a matrix/data.frame is provided, columns must refer to samples and rows to OTUs.
#' @import tidyverse
#' @return Normalised vector or matrix.
#' @seealso \code{\link{hill_div}}, \code{\link{index_div}}
#' @examples
#' data(bat.diet.otutable)
#' tss(bat.diet.otutable)
#' bat.diet.sample <- bat.diet.otutable[,1]
#' tss(bat.diet.sample)
#' @export

tss <- function(abund){
  #If input data is a vector
  if(is.null(dim(abund)) == TRUE){ abund.norm <- ifelse(is.nan(abund/sum(abund)), 0, abund/sum(abund)) }
  #If input data is an OTU table
  if(is.null(dim(abund)) == FALSE){
        abund.norm <- as.matrix(abund) %>% #convert to matrix to avoid sweep issues
                      sweep(., 2, colSums(.), FUN="/") %>%
                      replace(., is.nan(.), 0)} #convert NaNs derived from 0/0 to 0
  return(abund.norm)
}
