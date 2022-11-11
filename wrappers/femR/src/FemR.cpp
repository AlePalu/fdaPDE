#include<RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include<fdaPDE/Core.h>

// [[Rcpp::export]]
Rcpp::List rcpp_hello_world() {

    Rcpp::CharacterVector x = Rcpp::CharacterVector::create( "foo", "bar" )  ;
    Rcpp::NumericVector y   = Rcpp::NumericVector::create( 0.0, 1.0 ) ;
    Rcpp::List z            = Rcpp::List::create( x, y ) ;

    return z ;
}
