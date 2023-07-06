// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <fdaPDE/core/utils/Symbols.h>
#include <fdaPDE/core/FEM/PDE.h>
using fdaPDE::core::FEM::PDEInterface;
using fdaPDE::core::FEM::PDE;
#include <fdaPDE/core/FEM/operators/Laplacian.h>
using fdaPDE::core::FEM::Laplacian;
#include <fdaPDE/core/MESH/Mesh.h>
using fdaPDE::core::MESH::Mesh;

// this file contains the R wrapper for the SRPDE model

// better to use PDE_Type instead of int, but for now... ok. Need to find better solution
enum PDE_Type { LAPLACIAN, DIFFUSION_TRANSPORT_REACTION };

// M : local dimension, N : embedding dimension, R : FEM order
template <unsigned int M, unsigned int N, unsigned int R>
class R_PDE {
private:
  // internal data
  Mesh<M,N,R> domain_;
  PDEInterface<DMatrix<double>>* pde_;
public:
  // constructor
  R_PDE(const Rcpp::List& R_Mesh,
	int pde_type,
	const Rcpp::Nullable<Rcpp::List>& pde_parameters) :
    // initialize domain
    domain_(Rcpp::as<DMatrix<double>>(R_Mesh["nodes"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["edges"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["elements"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["neigh"]),
	    Rcpp::as<DMatrix<int>>   (R_Mesh["boundary"])){

    // define at run-time the pde object based on requested type
    switch(pde_type) {
    case 1:
      { // L = Laplacian
	auto L = Laplacian();
	pde_ = new PDE<M,N,R,decltype(L),DMatrix<double>>(domain_, L);
      }
      break;
    case 2:
      { // L = div(K*grad(f)) + b*grad(f) + c*f
	if(pde_parameters.isNotNull()) {
	  // extract parameters from pack
	  Rcpp::List pde_parameters_(pde_parameters);
	  SMatrix<M> K = Rcpp::as<DMatrix<double>>(pde_parameters_["diffusion"]);
	  SVector<M> b = Rcpp::as<DMatrix<double>>(pde_parameters_["transport"]);
	  double c = Rcpp::as<double>(pde_parameters_["reaction"]);
	  // define bilinear form
	  auto L = Laplacian(K) + Gradient(b) + c*Identity();
	  pde_ = new PDE<M,N,R,decltype(L),DMatrix<double>>(domain_, L);
	} else {
	  throw std::runtime_error("pde parameters not supplied.");
	}
      }
      break;
    }
    
  };
  
  // setters
  void set_dirichlet_bc(const DMatrix<double>& data) { pde_->setDirichletBC(data); }
  void set_force(const DMatrix<double>& data) { pde_->setForcing(data); }
  // getters
  DMatrix<double> get_quadrature_nodes() const { return pde_->quadratureNodes(); };
  DMatrix<double> get_dofs_coordinates() const { return domain_.dofCoords(); };

  SpMatrix<double> R0() const { return pde_->R0(); }
  SpMatrix<double> R1() const { return pde_->R1(); }
  DMatrix<double> solution() const { return pde_->solution(); }
  
  void init() { return pde_->init(); }
  void solve() { return pde_->solve(); }
  
  // destructor
  ~R_PDE() { delete pde_; }
};


// Rcpp module

// should remove FEM order from mesh... this would avoid to duplicate rcpp modules
// for fem order (I don't beleve it is possible to avoid it as template parameter...)
typedef R_PDE<2,2,1> PDE_2D_ORDER_1;
RCPP_MODULE(PDE_2D_ORDER_1) {
  Rcpp::class_<PDE_2D_ORDER_1>("PDE_2D_ORDER_1")
    .constructor<Rcpp::List, int, Rcpp::Nullable<Rcpp::List>>()
    // getters
    .method("get_quadrature_nodes", &PDE_2D_ORDER_1::get_quadrature_nodes)
    .method("get_dofs_coordinates", &PDE_2D_ORDER_1::get_dofs_coordinates)
    .method("get_mass",             &PDE_2D_ORDER_1::R0)
    .method("get_stiff",            &PDE_2D_ORDER_1::R1)
    .method("solution",             &PDE_2D_ORDER_1::solution)
    // setters
    .method("set_dirichlet_bc",     &PDE_2D_ORDER_1::set_dirichlet_bc)
    .method("set_force",            &PDE_2D_ORDER_1::set_force)
    // init and solve
    .method("init",                 &PDE_2D_ORDER_1::init)
    .method("solve",                &PDE_2D_ORDER_1::solve);
}

typedef R_PDE<2,2,2> PDE_2D_ORDER_2;
RCPP_MODULE(PDE_2D_ORDER_2) {
  Rcpp::class_<PDE_2D_ORDER_2>("PDE_2D_ORDER_2")
    .constructor<Rcpp::List, int, Rcpp::Nullable<Rcpp::List>>()
    // getters
    .method("get_quadrature_nodes", &PDE_2D_ORDER_2::get_quadrature_nodes)
    .method("get_dofs_coordinates", &PDE_2D_ORDER_2::get_dofs_coordinates)
    .method("get_mass",             &PDE_2D_ORDER_2::R0)
    .method("get_stiff",            &PDE_2D_ORDER_2::R1)
    .method("solution",             &PDE_2D_ORDER_2::solution)
    // setters
    .method("set_dirichlet_bc",     &PDE_2D_ORDER_2::set_dirichlet_bc)
    .method("set_force",            &PDE_2D_ORDER_2::set_force)
    // init and solve
    .method("init",                 &PDE_2D_ORDER_2::init)
    .method("solve",                &PDE_2D_ORDER_2::solve);
}
