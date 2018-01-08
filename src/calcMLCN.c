#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

SEXP calcMLCN(SEXP N, SEXP mu, SEXP sf, SEXP phi, SEXP prior) {
  
  /*
  SEXP result = PROTECT(allocVector(REALSXP, 4));
  double c_sf  = asReal(sf); 
  double c_phi = asReal(phi);
  double c_mu  = asReal(mu);
  double prob = 1/(c_mu*c_phi + 1);
  double size = c_sf/c_phi;
  REAL(result)[0] = asInteger(N);
  REAL(result)[1] = size;
  REAL(result)[2] = prob;
  REAL(result)[3] = pnbinom(asInteger(N), size, prob, 1, 0);
  UNPROTECT(1);
  return result;
  */
  
  double c_sf  = asReal(sf); 
  double c_phi = asReal(phi);
  SEXP result = PROTECT(allocVector(REALSXP, 4));
  
  REAL(result)[0] = c_sf*c_phi;
  UNPROTECT(1);
  return result;
  
}