#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h"
using fdaPDE::core::FEM::SpaceVaryingAdvection;
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
#include "core/MESH/Mesh.h"
#include "../fdaPDE/models/functional/fSRPDE.h"
using fdaPDE::models::FSRPDE;
#include "../fdaPDE/models/SamplingDesign.h"

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

#include <string>
#include <fstream>
#include <filesystem>

namespace Test_fSRPDE
{
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
}

/* test 1
   domain:       unit square [0,1] x [0,1]
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
 */
TEST(FSRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes)
{
    // define domain and regularizing PDE
    MeshLoader<Mesh2D<>> domain("unit_square");
    auto L = Laplacian();
    DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements() * 3, 1);
    PDE problem(domain.mesh, L, u); // definition of regularizing PDE

    // define statistical model
    FSRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem);

    // set lambda
    double lambda = std::pow(0.1, 4);
    model.setLambdaS(lambda);

    // load data from .csv files
    std::string test_directory = "data/models/FSRPDE/2D_test/";
    CSVReader<double> reader{};
    CSVFile<double> XFile; // observation file
    XFile = reader.parseFile(test_directory + "X.csv");
    DMatrix<double> X = XFile.toEigen();
    CSVFile<double> bFile; // observation file
    bFile = reader.parseFile(test_directory + "b.csv");
    DMatrix<double> b = bFile.toEigen();

    // set model data
    BlockFrame<double, int> df;
    df.insert(OBSERVATIONS_BLK, X);
    df.insert(DESIGN_MATRIX_BLK, b);
    model.setData(df);
    // alternative
    // model.setData(X);

    // solve smoothing problem
    model.init();
    model.solve();

    //  **  export results  ** //

    // std::string results_directory = test_directory + "results/";

    // if (!std::filesystem::exists(results_directory))
    //     std::filesystem::create_directory(results_directory);

    // std::ofstream nodesfile;
    // nodesfile.open(results_directory + "f.csv");
    // DMatrix<double> computedF = model.f();
    // nodesfile << computedF.format(Test_fSRPDE::CSVFormat);
    // nodesfile.close();
}