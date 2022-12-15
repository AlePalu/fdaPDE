#ifndef __STRPDE_H__
#define __STRPDE_H__

#include <memory>
#include <type_traits>
// CORE imports
#include "../../core/utils/Symbols.h"
#include "../../core/FEM/PDE.h"
using fdaPDE::core::FEM::PDEBase;
#include "../space_time/SplineBasis.h"
#include "../../core/NLA/SparseBlockMatrix.h"
using fdaPDE::core::NLA::SparseBlockMatrix;
#include "../../core/NLA/SMW.h"
using fdaPDE::core::NLA::SMW;
#include "../../core/NLA/KroneckerProduct.h"
using fdaPDE::core::NLA::SparseKroneckerProduct;
using fdaPDE::core::NLA::Kronecker;
// calibration module imports
#include "../../calibration/iGCV.h"
using fdaPDE::calibration::iGCV;
// regression module imports
#include "RegressionBase.h"
using fdaPDE::models::RegressionBase;

namespace fdaPDE{
namespace models{

  // base class for STRPDE model
  template <typename PDE, typename TimeRegularization, Sampling SamplingDesign> class STRPDE;

  template <typename PDE, Sampling SamplingDesign>
  class STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign>
    : public RegressionBase<STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign>>/*, public iGCV*/ {
    // compile time checks
    static_assert(std::is_base_of<PDEBase, PDE>::value);
  private:
    typedef SpaceTimeSeparableTag TimeRegularization;
    typedef RegressionBase<STRPDE<PDE, TimeRegularization, SamplingDesign>> Base;
    
    // system matrix of non-parametric problem (2N x 2N matrix)
    //     | -\Psi^T*D*W*\Psi  \lambda*R1^T |
    // A = |                                |
    //     |    \lambda*R1     \lambda*R0   |
    SpMatrix<double> A_{};
    // right hand side of problem's linear system (1 x 2N vector)
    //     | -\Psi^T*D*Q*y |
    // b = |               |, Q = W(I-H), H = X*(X^T*W*X)^{-1}*X^T*W
    //     |   \lambda*u   |
    DVector<double> b_{};

    // problem solution
    DMatrix<double> f_{};    // estimate of the spatial field (1 x N vector)
    DMatrix<double> g_{};    // PDE misfit
    DMatrix<double> beta_{}; // estimate of the coefficient vector (1 x q vector)
  public:
    // import commonly defined symbols from base
    IMPORT_REGRESSION_SYMBOLS;
    using Base::lambdaS; // smoothing parameter in space
    using Base::lambdaT; // smoothing parameter in time
    using Base::Pt;      // time penalization matrix: [Pt_]_{ij} = \int_{[0,T]} (\phi_i)_tt*(\phi_j)_tt
    using Base::Rt;      // time mass matrix: [Rt_]_{ij} = \int_{[0,T]} \phi_i*\phi_j
    
    // constructor
    STRPDE() = default;
    template <typename... SamplingData>
    STRPDE(const PDE& pde, const DMatrix<double>& time, const SamplingData&... s)
      : RegressionBase<STRPDE<PDE, SpaceTimeSeparableTag, SamplingDesign>>(pde, time, s...) {};
    
    // iStatModel interface implementation
    virtual void solve(); // finds a solution to the smoothing problem

    // iRegressionModel interface implementation
    virtual DMatrix<double> fitted();
    virtual double predict(const DVector<double>& covs, const std::size_t loc) const;
    // getters to problem solution
    virtual const DMatrix<double>& f() const { return f_; };
    virtual const DMatrix<double>& g() const { return g_; };
    virtual const DMatrix<double>& beta() const { return beta_; };
    
    // iGCV interface implementation
    // virtual const DMatrix<double>& T(); // T = \Psi^T*Q*\Psi + \lambda*(R1^T*R0^{-1}*R1)
    // virtual const DMatrix<double>& Q(); // Q = W(I - H) = W - W*X*(X^T*W*X)^{-1}X^T*W
    // returns the euclidian norm of y - \hat y
    // virtual double norm(const DMatrix<double>& obs, const DMatrix<double>& fitted) const;
    
    virtual ~STRPDE() = default;
  };

  // compile time informations related to the model
  template <typename PDE_, Sampling SamplingDesign>
  struct model_traits<STRPDE<PDE_, SpaceTimeSeparableTag, SamplingDesign>> {
    typedef PDE_ PDE;
    typedef SpaceTimeSeparableTag RegularizationType;
    typedef SplineBasis<3> TimeBasis; // use cubic B-splines
    static constexpr Sampling sampling = SamplingDesign;
  };

#include "STRPDE.tpp"
}}


#endif // __STRPDE_H__
