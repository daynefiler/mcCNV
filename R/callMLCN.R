#' @useDynLib mcCNV, .registration = TRUE, .fixes = "C_"
#' @export

callMLCN <- function(N, mu, sf, phi, prior) {
  
  .Call(C_calcMLCN, N, mu, sf, phi, prior)
  
}
