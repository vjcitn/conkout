#' MAF for TCGA GBM samples
#' @docType data
#' @format maftools MAF object
#' @source poisonAlien github TCGAMutations package
#' @examples
#' conkout::tcga_gbm
"tcga_gbm"
#' affy hgu133 expression data for TCGA GBM
#' @docType data
#' @format 12042 x 528 matrix
#' @source curatedTCGAData extract with full sample id as colname and hugo symbol as rowname
#' @examples
#' conkout::gbmu133
"gbmu133"
#' MSigDB-based glioblastoma gene sets, 47 in number
#' @docType data
#' @format 12042 x 528 matrix
#' @source imported using getGmt from an MSigDb download
#' @examples
#' conkout::glioSets47
#' "glioSets47"
