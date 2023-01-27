#include<RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include<utility>

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

#include <fdaPDE/core/FEM/operators/BilinearFormExpressions.h>
using fdaPDE::core::FEM::DefaultOperator;

template<unsigned int M, unsigned int N, unsigned int R, typename F>         
class R_PDE {

    private:
        using BilinearForm = typename std::remove_reference<F>::type;
        
        Mesh<M,N,R> mesh_;
        PDE<M,N,R, BilinearForm, DMatrix<double>> pde_;
    
        DMatrix<double> solution_;
    
    public:    
        R_PDE(const Rcpp::List& R_Mesh):mesh_(Rcpp::as<DMatrix<double>>(R_Mesh["nodes"]), Rcpp::as<DMatrix<int>>( R_Mesh["edges"]), Rcpp::as<DMatrix<int>>( R_Mesh["elements"]), 
                                              Rcpp::as<DMatrix<int>>(R_Mesh["neigh"]), Rcpp::as<DMatrix<int>>(R_Mesh["boundary"])), pde_(mesh_){}
        
        void set_dirichletBC(const Rcpp::NumericVector& R_BC){ 
            pde_.setDirichletBC(Rcpp::as<DMatrix<double>>(R_BC));
        };
        
        void set_forcingTerm(const Rcpp::NumericVector& R_forcingTerm){
            pde_.setForcing( Rcpp::as<DMatrix<double>>(R_forcingTerm));
        };
        
        void set_PDEparameters(const Rcpp::List& R_PDEparameters){
            
            double mu_ = Rcpp::as<double>(R_PDEparameters["diffusion"]);
            SVector<M> beta_ = Rcpp::as<DMatrix<double>>(R_PDEparameters["transport"]);
            double c_ = Rcpp::as<double>(R_PDEparameters["reaction"]);
            BilinearForm bilinearForm = mu_*Laplacian() + Gradient(beta_) + c_*Identity();
            pde_.setBilinearForm(bilinearForm);  
        };
        
        DMatrix<double> get_quadrature_nodes() const {
            return pde_.integrator().quadratureNodes(mesh_);
        };
        
        DMatrix<double> get_dofs_coordinates() const {
            return mesh_.dofCoords();
        };
        
        DMatrix<double> force() const{
            return pde_.force();
        };
        
        Rcpp::List solve(){
            
            pde_.init();
            pde_.solve();
            
            DMatrix<double> solution   = pde_.solution();
            SpMatrix<double> Stiffness = pde_.R1();
            SpMatrix<double> Mass = pde_.R0();
         
            return Rcpp::List::create(Rcpp::Named("solution")     = solution,
                                      Rcpp::Named("Stiffness")    = Stiffness,
                                      Rcpp::Named("Mass")         = Mass);
        }; 
};

Rcpp::NumericMatrix interval2network(const Rcpp::NumericMatrix& interval){
    Rcpp::NumericMatrix output(interval.size(), 2);
    output(Rcpp::_,0) = interval;
    return output;
};

SpMatrix<int> fromTriplets(const Rcpp::IntegerMatrix& R_neigh, const Rcpp::NumericMatrix& interval){
    int nnodes = interval.size();
    std::vector<Eigen::Triplet<int>> triplets;
    triplets.reserve(R_neigh.rows());

    for(int i=0; i<R_neigh.rows(); ++i){
        // realign indexes (we assume index coming from mesh generator to be greater or equal to 1, C++ starts count from 0)
        triplets.push_back( Eigen::Triplet<int>( R_neigh(i,0)-1, R_neigh(i,1)-1, R_neigh(i,2)) );
    }

    SpMatrix<int> neigh(nnodes, nnodes);
    neigh.setFromTriplets(triplets.begin(), triplets.end());
    return(neigh);
};

template<unsigned int R, typename F>         
class R_PDE<1,2,R,F> {

    private:
        using BilinearForm = typename std::remove_reference<F>::type;
        
        Mesh<1,2,R> mesh_;
        PDE<1,2,R, BilinearForm, DMatrix<double>> pde_;
    
        DMatrix<double> solution_;
    
    public:    
        R_PDE(const Rcpp::List& R_Mesh):mesh_(Rcpp::as<DMatrix<double>>(interval2network(R_Mesh["nodes"])), 
                                              Rcpp::as<DMatrix<int>>( R_Mesh["edges"]), 
                                              Rcpp::as<DMatrix<int>>( R_Mesh["elements"]), 
                                              fromTriplets(R_Mesh["neigh"], R_Mesh["nodes"]), 
                                              Rcpp::as<DMatrix<int>>(R_Mesh["boundary"])), pde_(mesh_){}
        
        void set_dirichletBC(const Rcpp::NumericVector& R_BC){ 
            pde_.setDirichletBC(Rcpp::as<DMatrix<double>>(R_BC));
        };
        
        void set_forcingTerm(const Rcpp::NumericVector& R_forcingTerm){
            pde_.setForcing( Rcpp::as<DMatrix<double>>(R_forcingTerm));
        };
        
        void set_PDEparameters(const Rcpp::List& R_PDEparameters){
            
            double mu_ = Rcpp::as<double>(R_PDEparameters["diffusion"]);
            
            double c_ = Rcpp::as<double>(R_PDEparameters["reaction"]);
            BilinearForm bilinearForm = mu_*Laplacian() + c_*Identity(); 
            pde_.setBilinearForm(bilinearForm);  
        };
        
        DMatrix<double> get_quadrature_nodes() const {
            return pde_.integrator().quadratureNodes(mesh_);
        };
        
        DMatrix<double> get_dofs_coordinates() const {
            return mesh_.dofCoords();
        };
        
        DMatrix<double> force() const{
            return pde_.force();
        };
        
        Rcpp::List solve(){
            
            pde_.init();
            pde_.solve();
            
            DMatrix<double> solution   = pde_.solution();
            SpMatrix<double> Stiffness = pde_.R1();
            SpMatrix<double> Mass = pde_.R0();
            
            return Rcpp::List::create(Rcpp::Named("solution")     = solution,
                                      Rcpp::Named("Stiffness")    = Stiffness,
                                      Rcpp::Named("Mass")         = Mass);
        }; 
};

template<unsigned int M> using PDE_isotropic = decltype( std::declval<double>()*std::declval<Laplacian<DefaultOperator>>() + 
                                                         std::declval<Gradient<SVector<M>>>() + 
                                                         std::declval<double>()*std::declval<Identity<DefaultOperator>>() ); 
template<unsigned int M, unsigned int N, unsigned int R> using R_PDE_isotropic = R_PDE<M,N,R, PDE_isotropic<M>>;

using PDE_Diff_Reac_isotropic = decltype( std::declval<double>()*std::declval<Laplacian<DefaultOperator>>() + 
                                std::declval<double>()*std::declval<Identity<DefaultOperator>>() ); 

using R_PDE_isotropic_1D_ORDER_1 = R_PDE_isotropic<1,1,1>;
using R_PDE_isotropic_1D_ORDER_2 = R_PDE_isotropic<1,1,2>;

using R_PDE_isotropic_1_5D_ORDER_1 = R_PDE<1, 2, 1, PDE_Diff_Reac_isotropic>;
using R_PDE_isotropic_1_5D_ORDER_2 = R_PDE<1, 2, 2, PDE_Diff_Reac_isotropic>;

using R_PDE_isotropic_2D_ORDER_1 = R_PDE_isotropic<2,2,1>;
using R_PDE_isotropic_2D_ORDER_2 = R_PDE_isotropic<2,2,2>;

RCPP_MODULE(PDE_1D_isotropic_ORDER_1) {
  Rcpp::class_<R_PDE_isotropic_1D_ORDER_1>("PDE_1D_isotropic_ORDER_1")
    .constructor<Rcpp::List>()
    .method("get_quadrature_nodes", &R_PDE_isotropic_1D_ORDER_1::get_quadrature_nodes)
    .method("set_dirichletBC",      &R_PDE_isotropic_1D_ORDER_1::set_dirichletBC)
    .method("set_forcingTerm",      &R_PDE_isotropic_1D_ORDER_1::set_forcingTerm)
    .method("set_PDEparameters",    &R_PDE_isotropic_1D_ORDER_1::set_PDEparameters)
    .method("get_dofs_coordinates", &R_PDE_isotropic_1D_ORDER_1::get_dofs_coordinates)
    .method("force",                &R_PDE_isotropic_1D_ORDER_1::force)
    .method("solve",                &R_PDE_isotropic_1D_ORDER_1::solve)
    ;
}

RCPP_MODULE(PDE_1D_isotropic_ORDER_2) {
  Rcpp::class_<R_PDE_isotropic_1D_ORDER_2>("PDE_1D_isotropic_ORDER_2")
    .constructor<Rcpp::List>()
    .method("get_quadrature_nodes", &R_PDE_isotropic_1D_ORDER_2::get_quadrature_nodes)
    .method("set_dirichletBC",      &R_PDE_isotropic_1D_ORDER_2::set_dirichletBC)
    .method("set_forcingTerm",      &R_PDE_isotropic_1D_ORDER_2::set_forcingTerm)
    .method("set_PDEparameters",    &R_PDE_isotropic_1D_ORDER_2::set_PDEparameters)
    .method("get_dofs_coordinates", &R_PDE_isotropic_1D_ORDER_2::get_dofs_coordinates)
    .method("force",                &R_PDE_isotropic_1D_ORDER_2::force)
    .method("solve",                &R_PDE_isotropic_1D_ORDER_2::solve)
    ;
}

RCPP_MODULE(PDE_1_5D_isotropic_ORDER_1) {
  Rcpp::class_<R_PDE_isotropic_1_5D_ORDER_1>("PDE_1_5D_isotropic_ORDER_1")
    .constructor<Rcpp::List>()
    .method("get_quadrature_nodes", &R_PDE_isotropic_1_5D_ORDER_1::get_quadrature_nodes)
    .method("set_dirichletBC",      &R_PDE_isotropic_1_5D_ORDER_1::set_dirichletBC)
    .method("set_forcingTerm",      &R_PDE_isotropic_1_5D_ORDER_1::set_forcingTerm)
    .method("set_PDEparameters",    &R_PDE_isotropic_1_5D_ORDER_1::set_PDEparameters)
    .method("get_dofs_coordinates", &R_PDE_isotropic_1_5D_ORDER_1::get_dofs_coordinates)
    .method("force",                &R_PDE_isotropic_1_5D_ORDER_1::force)
    .method("solve",                &R_PDE_isotropic_1_5D_ORDER_1::solve)
    ;
}

RCPP_MODULE(PDE_1_5D_isotropic_ORDER_2) {
  Rcpp::class_<R_PDE_isotropic_1_5D_ORDER_2>("PDE_1_5D_isotropic_ORDER_2")
    .constructor<Rcpp::List>()
    .method("get_quadrature_nodes", &R_PDE_isotropic_1_5D_ORDER_2::get_quadrature_nodes)
    .method("set_dirichletBC",      &R_PDE_isotropic_1_5D_ORDER_2::set_dirichletBC)
    .method("set_forcingTerm",      &R_PDE_isotropic_1_5D_ORDER_2::set_forcingTerm)
    .method("set_PDEparameters",    &R_PDE_isotropic_1_5D_ORDER_2::set_PDEparameters)
    .method("get_dofs_coordinates", &R_PDE_isotropic_1_5D_ORDER_2::get_dofs_coordinates)
    .method("force",                &R_PDE_isotropic_1_5D_ORDER_2::force)
    .method("solve",                &R_PDE_isotropic_1_5D_ORDER_2::solve)
    ;
}

RCPP_MODULE(PDE_2D_isotropic_ORDER_1) {
  Rcpp::class_<R_PDE_isotropic_2D_ORDER_1>("PDE_2D_isotropic_ORDER_1")
    .constructor<Rcpp::List>()
    .method("get_quadrature_nodes", &R_PDE_isotropic_2D_ORDER_1::get_quadrature_nodes)
    .method("set_dirichletBC",      &R_PDE_isotropic_2D_ORDER_1::set_dirichletBC)
    .method("set_forcingTerm",      &R_PDE_isotropic_2D_ORDER_1::set_forcingTerm)
    .method("set_PDEparameters",    &R_PDE_isotropic_2D_ORDER_1::set_PDEparameters)
    .method("get_dofs_coordinates", &R_PDE_isotropic_2D_ORDER_1::get_dofs_coordinates)
    .method("force",                &R_PDE_isotropic_2D_ORDER_1::force)
    .method("solve",                &R_PDE_isotropic_2D_ORDER_1::solve)
    ;
}

RCPP_MODULE(PDE_2D_isotropic_ORDER_2) {
  Rcpp::class_<R_PDE_isotropic_2D_ORDER_2>("PDE_2D_isotropic_ORDER_2")
    .constructor<Rcpp::List>()
    .method("get_quadrature_nodes", &R_PDE_isotropic_2D_ORDER_2::get_quadrature_nodes)
    .method("set_dirichletBC",      &R_PDE_isotropic_2D_ORDER_2::set_dirichletBC)
    .method("set_forcingTerm",      &R_PDE_isotropic_2D_ORDER_2::set_forcingTerm)
    .method("set_PDEparameters",    &R_PDE_isotropic_2D_ORDER_2::set_PDEparameters)
    .method("get_dofs_coordinates", &R_PDE_isotropic_2D_ORDER_2::get_dofs_coordinates)
    .method("force",                &R_PDE_isotropic_2D_ORDER_2::force)
    .method("solve",                &R_PDE_isotropic_2D_ORDER_2::solve)
    ;
}



