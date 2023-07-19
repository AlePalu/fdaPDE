// solution in case of fixed \lambda
template <typename PDE, typename RegularizationType,
          typename SamplingDesign, typename lambda_selection_strategy>
void FPCA_CS<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(fixed_lambda)
{
    // Data
    DMatrix<double> X = data().template get<double>(OBSERVATIONS_BLK);

    // Penalty matrix
    // std::cout << "Penalty matrix" << std::endl;
    if (mass_lumping_)
    {
        unsigned int K = n_basis();
        SpMatrix<double> invR0_;
        invR0_.resize(K, K);
        invR0_.reserve(K);
        for (std::size_t i = 0; i < K; ++i)
        {
            invR0_.insert(i, i) = 1 / R0().col(i).sum();
        }
        invR0_.makeCompressed();
        P_ = R1() * invR0_ * R1();
    }
    else
    {
        fdaPDE::SparseLU<SpMatrix<double>> invR0_;
        invR0_.compute(R0());
        P_ = R1() * invR0_.solve(R1());
    }

    // Resolution
    // std::cout << "RSVD" << std::endl;
    RSVD rsvd(X, lambda()[0], n_pc_, Psi(not_nan()), P_);

    // for (std::size_t i = 0; i < n_pc_; i++)
    // {
    //     // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda

    //     rsvd.solve();

    //     loadings_.col(i) = rsvd.loadings().transpose();
    //     double f_n_norm = std::sqrt(loadings_.col(i).dot(R0() * loadings_.col(i)));
    //     loadings_.col(i) = loadings_.col(i) / f_n_norm;
    //     scores_.col(i) = rsvd.scores() * f_n_norm;

    //     // subtract computed PC from data
    //     X -= scores_.col(i) * loadings_.col(i).transpose();
    //     rsvd.set_data(X);
    // }

    rsvd.solve();

    // std::cout << "Loadings" << std::endl;
    loadings_ = Psi(not_nan()) * (rsvd.loadings()).transpose();
    // std::cout << "Scores" << std::endl;
    scores_ = rsvd.scores();

    // std::cout << "Normalization" << std::endl;
    for (std::size_t i = 0; i < n_pc_; i++)
    {
        double f_n_norm;
        if constexpr (is_space_only<decltype(*this)>::value)
            f_n_norm = std::sqrt(loadings_.col(i).dot(R0() * loadings_.col(i)));
        else
            f_n_norm = loadings_.col(i).norm();
        // std::cout << f_n_norm << std::endl;
        loadings_.col(i) /= f_n_norm;
        scores_.col(i) *= f_n_norm;
    }

    return;
}

// best \lambda for PC choosen according to GCV index
template <typename PDE, typename RegularizationType,
          typename SamplingDesign, typename lambda_selection_strategy>
void FPCA_CS<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(gcv_lambda_selection)
{
    /*
    ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
    BlockFrame<double, int> data_ = data();
    // wrap GCV into a ScalarField accepted by OPT module
    const std::size_t n_lambda = n_smoothing_parameters<RegularizationType>::value;
    ScalarField<n_lambda> f;
    f = [&pe, &data_](const SVector<n_lambda> &p) -> double
    {
        // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f) fixed \lambda = p
        pe.compute(data_, p);
        return pe.gcv(); // return GCV at convergence
    };
    GridOptimizer<n_lambda> opt; // optimization algorithm
    // Principal Components computation
    for (std::size_t i = 0; i < n_pc_; i++)
    {
        opt.optimize(f, lambdas()); // select optimal \lambda for i-th PC
        // compute and store results given estimated optimal \lambda
        pe.compute(data_, opt.optimum());
        loadings_.col(i) = pe.f_n() / pe.f_n_norm();
        scores_.col(i) = pe.s() * pe.f_n_norm();
        // subtract computed PC from data
        data_.get<double>(OBSERVATIONS_BLK) -= scores_.col(i) * loadings_.col(i).transpose();
    }
    return;
    */
}

// best \lambda for PC choosen according to K-fold CV strategy, uses the reconstruction error on test set as CV score
template <typename PDE, typename RegularizationType,
          typename SamplingDesign, typename lambda_selection_strategy>
void FPCA_CS<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve_(kcv_lambda_selection)
{
    /*
    ProfilingEstimation<decltype(*this)> pe(*this, tol_, max_iter_);
    // number of smoothing parameters
    const std::size_t n_lambda = n_smoothing_parameters<RegularizationType>::value;

    // routine executed by the CV-engine to produce the model score
    std::function<double(DVector<double>, BlockFrame<double, int>, BlockFrame<double, int>)> cv_score =
        [&pe, this](const DVector<double> &lambda,
                    const BlockFrame<double, int> &train_df,
                    const BlockFrame<double, int> &test_df) -> double
    {
        // get references to train and test sets
        const DMatrix<double> &X_test = test_df.get<double>(OBSERVATIONS_BLK);

        SVector<n_lambda> p(lambda.data());
        // find vectors s,f minimizing \norm_F{Y - s^T*f}^2 + (s^T*s)*P(f)
        pe.compute(train_df, p);
        // compute reconstruction error and scores (X_test * f_n)/(\norm{f_n} + \lambda*P(f)) on test set
        if (this->has_nan())
        {
            auto nan_pattern = X_test.array().isNaN(); // missingness pattern
            std::size_t n = nan_pattern.count();       // number of not-NaN points in test set

            // scores vector
            DVector<double> s = (nan_pattern.select(0, X_test)) * pe.f_n();
            double pen = p[0] * (pe.g().dot(R0() * pe.g())); // \lambda*P(f)
            for (std::size_t i = 0; i < s.rows(); ++i)
            {
                // compute \norm{f_n} for i-th subject
                double sum_f = (X_test.row(i)).array().isNaN().select(0, pe.f_n().transpose()).squaredNorm();
                s[i] /= (sum_f + pen); // normalize i-th score
            }

            // evaluate reconstruction error on test set
            double err = nan_pattern.select(0, X_test - s * pe.f_n().transpose()).squaredNorm();
            return err / n;
        }
        else
        {
            // scores vector
            DVector<double> s = (X_test * pe.f_n()) / (pe.f_n().squaredNorm() + p[0] * (pe.g().dot(R0() * pe.g())));
            // evaluate reconstruction error on test set
            return (X_test - s * pe.f_n().transpose()).squaredNorm() / X_test.size();
        }
    };

    // define K-fold algorithm
    KFoldCV cv(10); // allow user-defined number of folds!
    std::vector<DVector<double>> lambdas_;
    lambdas_.reserve(lambdas().size());
    for (const auto &l : lambdas())
        lambdas_.emplace_back(Eigen::Map<const DVector<double>>(l.data(), n_lambda, 1));
    // Principal Components computation
    BlockFrame<double, int> data_ = data();
    for (std::size_t i = 0; i < n_pc_; i++)
    {
        cv.compute(lambdas_, data_, cv_score, false); // select optimal smoothing level
        // compute and store results given estimated optimal \lambda
        pe.compute(data_, cv.optimum());
        loadings_.col(i) = pe.f_n() / pe.f_n_norm();
        scores_.col(i) = pe.s() * pe.f_n_norm();
        // subtract computed PC from data
        data_.get<double>(OBSERVATIONS_BLK) -= pe.s() * pe.f_n().transpose();
    }
    */
    return;
}

// finds solution to FPCA_CS problem, dispatch to solver depending on \lambda selection criterion
template <typename PDE, typename RegularizationType,
          typename SamplingDesign, typename lambda_selection_strategy>
void FPCA_CS<PDE, RegularizationType, SamplingDesign, lambda_selection_strategy>::solve()
{
    // pre-allocate space
    loadings_.resize(X().cols(), n_pc_);
    scores_.resize(X().rows(), n_pc_);

    // dispatch to desired solution strategy
    solve_(lambda_selection_strategy());
    return;
}
