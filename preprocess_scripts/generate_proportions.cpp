// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//marginal <- function(dat, combos, ncol, colnames) {
//  d <- vector("list", length(dat))
//  for( i in 1:length(dat) ) {
//    cat("Iteration", i, "of", length(dat), ".\n")
//    x <- dat[[i]]
//    output <- matrix(FALSE, nrow=nrow(x), ncol=ncol)
//    for( j in 1:ncol(output) ) {
//      output[, j] <- as.logical(apply( x[, combos[[j]], drop=FALSE], 1, prod))
//    }
//    colnames(output) <- colnames
//    d[[i]] <- output
//  }
//  return(d)
//}

// [[Rcpp::export]]
List marginalRcpp(List dat, List combos, int ncol, CharacterVector colnames) {
  
  int dat_len = dat.size();
  List d(dat_len);
  for (int i=0; i < dat_len; ++i) {
    Rcout << "Iteration " << i+1 << " of " << dat_len << "." << std::endl;
    mat x = as<mat>(dat[i]);
    int nrow = x.n_rows;
    LogicalMatrix output(nrow, ncol);
    output.attr("dimnames") = List::create( R_NilValue, colnames );
    for (int j=0; j < ncol; ++j) {
      LogicalVector to_insert(nrow);
      uvec vecs_to_take = as<uvec>( combos[j] );
      vecs_to_take = vecs_to_take - 1; // 0-based indexing
      mat tmp = x.cols(vecs_to_take);
      for (int k=0; k < nrow; ++k) {
        to_insert[k] = prod( tmp.row(k) ) > 0;
      }
      output(_, j) = to_insert;
    }
    d[i] = output;
  }
  return d;
}

//joint <- function(dat, combos, ncol, colnames) {
//  
//  ## convert dat to binary
//  dat[] <- lapply(dat, function(x) {
//    colApply(x, as.logical, drop=FALSE)
//  })
//  
//  d <- vector("list", length(dat))
//  for( i in 1:length(dat) ) {
//    cat("Iteration", i, "of", length(dat), ".\n")
//    x <- dat[[i]]
//    output <- matrix(FALSE, nrow=nrow(x), ncol=ncol)
//    for( j in 1:ncol(output) ) {
//      check <- rep(FALSE, ncol(x))
//      check[ combos[[j]] ] <- TRUE
//      output[, j] <- apply(x, 1, function(x) { all(x == check) })
//    }
//    colnames(output) <- colnames
//    d[[i]] <- output
//  }
//  return(d)
//}

// [[Rcpp::export]]
List jointRcpp(List dat, List combos, int ncol, CharacterVector colnames) {
  
  int dat_len = dat.size();
  List d(dat_len);
  for (int i=0; i < dat_len; ++i) {
    Rcout << "Iteration " << i+1 << " of " << dat_len << "." << std::endl;
    LogicalMatrix x = as<LogicalMatrix>( dat[i] );
    int nrow = x.nrow();
    int ncol_x = x.ncol();
    LogicalMatrix output(nrow, ncol);
    for (int j=0; j < ncol; ++j) {
      LogicalVector check(ncol_x);
      IntegerVector tmp = combos[j];
      for (int k=0; k < tmp.size(); ++k) {
        check[ tmp[k]-1 ] = TRUE;
      }
      
      // debugging
//      Rcout << "(";
//      for (int b=0; b < check.size()-1; ++b) {
//        Rcout << check[b] << ", ";
//      }
//      Rcout << check[check.size()-1] << ")" << std::endl;
      
      for (int k=0; k < nrow; ++k) {
        for (int m=0; m < ncol_x; ++m) {
          if (x(k, m) != check[m]) {
            output(k, j) = FALSE;
            break;
          }
          output(k, j) = TRUE;
        }
      }
      
      // dimnames
      output.attr("dimnames") = List::create( R_NilValue, colnames );
    }
    d[i] = output;
  }
  return d;
}
