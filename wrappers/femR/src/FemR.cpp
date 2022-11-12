#include<RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include<fdaPDE/Core.h>
using fdaPDE::core::ScalarField;

#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;
using fdaPDE::core::MESH::neighboring_structure;

#include <fdaPDE/core/FEM/PDE.h>
using fdaPDE::core::FEM::PDE;

#include <fdaPDE/core/FEM/Evaluator.h>
using fdaPDE::core::FEM::Evaluator;
using fdaPDE::core::FEM::Raster;

#include <fdaPDE/core/utils/CompileTime.h>
#include <fdaPDE/core/FEM/basis/BasisTable.h>
using fdaPDE::core::FEM::BASIS_TABLE;

/*
#include <fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h>
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingAdvection;
using fdaPDE::core::FEM::SpaceVaryingReaction;
*/

template<unsigned int M, unsigned int N, unsigned int R>
Mesh<M,N,R> create_mesh(const Rcpp::NumericMatrix& points, const Rcpp::IntegerMatrix& edges,
            const Rcpp::IntegerMatrix& elements, const Rcpp::IntegerMatrix& neigh, const Rcpp::IntegerMatrix boundary){
            
    Mesh<M,N,R> mesh_(Rcpp::as<DMatrix<double>>(points), 
                      Rcpp::as<DMatrix<int>>(edges), 
                      Rcpp::as<DMatrix<int>>(elements), 
                      Rcpp::as<DMatrix<int>>(neigh), 
                      Rcpp::as<DMatrix<int>>(boundary));            
    return mesh_;
}

// [[Rcpp::export]]
Rcpp::List rcpp_hello_world(const Rcpp::NumericMatrix& points, const Rcpp::IntegerMatrix& edges,
            const Rcpp::IntegerMatrix& elements, const Rcpp::IntegerMatrix& neigh, const Rcpp::IntegerMatrix boundary, 
            const double& diffusion, const Rcpp::NumericMatrix& transport, const double& reaction) {
    
    auto mesh_ = create_mesh<2,2,1>(points, edges, elements, neigh, boundary);
    auto L = diffusion*Laplacian() + Gradient(Rcpp::as<DVector<double>>(transport)) + reaction*Identity();
    DMatrix<double> u; // forcing term u
    u.resize(mesh_.elements()*3, 1);
    u.fill(0);
    
    // define differential problem
    PDE problem(mesh_, L, u);
    problem.init();
    problem.solve();
    
    Rcpp::List z = Rcpp::List::create( Rcpp::Named("solution") = Rcpp::NumericMatrix(Rcpp::wrap((*problem.solution())))) ;

    return z ;
}
