#include <cstddef>
#include <gtest/gtest.h> // testing framework
#include <unsupported/Eigen/SparseExtra>

#include "../fdaPDE/core/utils/Symbols.h"
#include "../fdaPDE/core/utils/IO/CSVReader.h"
#include "../fdaPDE/core/FEM/PDE.h"
using fdaPDE::core::FEM::PDE;
#include "../fdaPDE/core/FEM/operators/SpaceVaryingFunctors.h"
using fdaPDE::core::FEM::SpaceVaryingDiffusion;
using fdaPDE::core::FEM::SpaceVaryingAdvection;
#include "core/MESH/Mesh.h"
#include "../fdaPDE/models/regression/SRPDE.h"
using fdaPDE::models::SRPDE;
#include "../fdaPDE/models/regression/SQRPDE.h"
using fdaPDE::models::SQRPDE;

#include "../fdaPDE/models/regression/GSRPDE.h"
using fdaPDE::models::GSRPDE;

#include "../fdaPDE/models/SamplingDesign.h"
#include "../fdaPDE/calibration/GCV.h"
using fdaPDE::calibration::GCV;
using fdaPDE::calibration::ExactGCV;
using fdaPDE::calibration::ExactEDF;
using fdaPDE::calibration::StochasticEDF;
#include "../fdaPDE/core/OPT/optimizers/GridOptimizer.h"
using fdaPDE::core::OPT::GridOptimizer;
#include "../fdaPDE/core/OPT/optimizers/Newton.h"
using fdaPDE::core::OPT::NewtonOptimizer;

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

// for time and memory performances
#include <chrono>
#include <iomanip>
using namespace std::chrono;
#include <unistd.h>
#include <fstream>


// /* test 1
//    domain:       unit square [0,1] x [0,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//  */
// TEST(GCV_SRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes_GridExact) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("unit_square_coarse");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   // define statistical model
//   SRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes>  model(problem);
  
//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile("data/models/SRPDE/2D_test6/y.csv");
//   DMatrix<double> y = yFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -6.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV function and optimize
//   GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();
  
//   // expected values for q + Tr[S]
//   std::vector<double> expected_edfs = {
//     432.0526271699909557, 425.5738291798880368, 414.9490914902157783, 398.3650445980014752, 374.2000509470916541,
//     341.8926575588438936, 302.6569434589166576, 259.4124363611769581, 215.8693404067796564, 175.3273544830321384,
//     139.8641263839342344, 110.2252857831315538, 86.2049347912456341
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );
  
//   // expected value of gcv(\lambda)
//   std::vector<double> expected_gcvs = {
//     0.0400234253534814, 0.0397604316008810, 0.0394198736081135, 0.0391060717299170, 0.0391008105529502,
//     0.0398771699638851, 0.0419425162273298, 0.0456080838211829, 0.0509765872233825, 0.0581903455597975,
//     0.0675560664911411, 0.0793326714839587, 0.0934793416959190
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

//   // check optimal lambda
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[4][0]) );
// }

// /* test 2
//    domain:       unit square [0,1] x [0,1] (coarse)
//    sampling:     locations = nodes
//    penalization: simple laplacian
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
//  */
// TEST(GCV_SRPDE, Test2_Laplacian_NonParametric_GeostatisticalAtNodes_GridStochastic) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("unit_square_coarse");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   // define statistical model
//   SRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes>  model(problem);
  
//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile("data/models/SRPDE/2D_test6/y.csv");
//   DMatrix<double> y = yFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -6.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV calibrator
//   std::size_t seed = 476813;
//   GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 1000, seed);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();

//   // expected values for q + Tr[S]
//   std::vector<double> expected_edfs = {
//     432.0405848764650045, 425.5529143800780503, 414.9133694755933561, 398.3056506911505608, 374.1052725793284139,
//     341.7500896866504831, 302.4588492694787192, 259.1628602033034667, 215.5870958975179690, 175.0398481074188624,
//     139.5961352594555080, 109.9912391645740399, 86.0073046035169710 
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)
//   std::vector<double> expected_gcvs = {
//     0.0399159071906597, 0.0396528360942592, 0.0393119874506983, 0.0389973432053000, 0.0389900907451345,
//     0.0397626888527567, 0.0418226582836935, 0.0454829731994149, 0.0508490093083462, 0.0580646045230030,
//     0.0674359858974278, 0.0792205235204336, 0.0933752877387227
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

//   // check optimal lambda
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[4][0]) );
// }

// /* test 3
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//  */
// TEST(GCV_SRPDE, Test3_Laplacian_SemiParametric_GeostatisticalAtLocations_GridExact) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("c_shaped");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   // define statistical model
//   CSVReader<double> reader{};
//   // load locations where data are sampled
//   CSVFile<double> locFile;
//   locFile = reader.parseFile("data/models/SRPDE/2D_test2/locs.csv");
//   DMatrix<double> loc = locFile.toEigen();

//   // define statistical model
//   SRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem);
//   model.set_spatial_locations(loc);
  
//   // load data from .csv files
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile  ("data/models/SRPDE/2D_test2/z.csv");
//   DMatrix<double> y = yFile.toEigen();
//   CSVFile<double> XFile; // design matrix
//   XFile = reader.parseFile  ("data/models/SRPDE/2D_test2/X.csv");
//   DMatrix<double> X = XFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK,  y);
//   df.insert(DESIGN_MATRIX_BLK, X);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -3.0; x <= 3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV calibrator
//   GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//   GridOptimizer<1> opt;
  
//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();

//   // expected values for q + Tr[S]
//   std::vector<double> expected_edfs = {
//     109.0637499040006588, 93.2229476278395026, 78.0169165468049499, 64.2577092364702054, 52.3413662612859625,
//      42.3404295978456844, 34.1086127469056137, 27.3954638966633190, 21.9573643522161781, 17.6119343172742404,
//      14.2191625526008725, 11.6385106137848773,  9.7178138350592107,  8.3098853216415982,  7.2842171925382466,
//       6.5281267228208399,  5.9498450377698173,  5.4843127845976953,  5.0921215492618579,  4.7522541472780739,
//       4.4551261509828262,  4.1946327397326151,  3.9615593551368464,  3.7465168923635441,  3.5489031263616644
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)
//   std::vector<double> expected_gcvs = {
//     2.1771484290002601, 2.0318749168846373, 1.9229411153254741, 1.8376099983427667, 1.7683487774209516,
//     1.7103300463087854, 1.6607101851375488, 1.6184324970362900, 1.5835055354976433, 1.5560348751711195,
//     1.5355021206088371, 1.5209371474466595, 1.5116339095670082, 1.5071872569869649, 1.5067140204400502,
//     1.5086445845177732, 1.5116269898476760, 1.5153441560837564, 1.5210973148544640, 1.5339264639588630,
//     1.5682200504833339, 1.6582436157224418, 1.8695409528491944, 2.2884361440774796, 2.9586790292370440
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

//   // check optimal lambda
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[14][0]) );
// }

// /* test 4
//    domain:       c-shaped
//    sampling:     locations != nodes
//    penalization: simple laplacian
//    covariates:   yes
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
//  */
// TEST(GCV_SRPDE, Test4_Laplacian_SemiParametric_GeostatisticalAtLocations_GridStochastic) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("c_shaped");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   // define statistical model
//   CSVReader<double> reader{};
//   // load locations where data are sampled
//   CSVFile<double> locFile;
//   locFile = reader.parseFile("data/models/SRPDE/2D_test2/locs.csv");
//   DMatrix<double> loc = locFile.toEigen();

//   // define statistical model
//   SRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem);
//   model.set_spatial_locations(loc);
  
//   // load data from .csv files
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile  ("data/models/SRPDE/2D_test2/z.csv");
//   DMatrix<double> y = yFile.toEigen();
//   CSVFile<double> XFile; // design matrix
//   XFile = reader.parseFile  ("data/models/SRPDE/2D_test2/X.csv");
//   DMatrix<double> X = XFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK,  y);
//   df.insert(DESIGN_MATRIX_BLK, X);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -3.0; x <= 3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
  
//   // define GCV calibrator
//   std::size_t seed = 66546513;
//   GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();

//   // expected values for q + Tr[S] (approximated)
//   std::vector<double> expected_edfs = {
//     106.2709775815650062, 91.2781878154271595, 76.8216314734006431, 63.6368232668402314, 52.0874830058305562,
//     42.2619783483174274,  34.0612617815139203, 27.2926608538180062, 21.7650509877889817, 17.3309390327612363,
//     13.8629820650085431,  11.2197970258596094,  9.2470064712368139,  7.7975062657152465,  6.7408531721908309,
//     5.9632650796519151,    5.3718179354249518,  4.9012488716219842,  4.5124258626382030,  4.1837705547537087,
//     3.9038225985944632,    3.6643245426561473,  3.4547903696710578,  3.2650834389510122,  3.0930648940657273
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)
//   std::vector<double> expected_gcvs = {
//     1.9933325455157656, 1.9339516163142583, 1.8744400183178187, 1.8164174958658492, 1.7608058393404040,
//     1.7082461763764853, 1.6595618959678919, 1.6161177888901777, 1.5794269592455561, 1.5503494275186216,
//     1.5285490034180831, 1.5129761785234856, 1.5028470225931159, 1.4977383948589298, 1.4967621840119401,
//     1.4983351571813046, 1.5010945226728243, 1.5047240675196598, 1.5105230912061263, 1.5234895155928254,
//     1.5578890056928336, 1.6477498970873763, 1.8582485673557088, 2.2753184385488714, 2.9426362338294938
//   };  
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

//   // check optimal lambda
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[14][0]) );
// }



// /* test 5
//    domain:       unit square [0,1] x [0,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid exact
//  */
// TEST(GCV_SRPDE, Test5_CostantCoefficientsPDE_NonParametric_GeostatisticalAtNodes_GridExact) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("unit_square_coarse");

//   // non unitary diffusion tensor
//   SMatrix<2> K;
//   K << 1,0,0,4;
//   auto L = Laplacian(K); // anisotropic diffusion
  
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   // define statistical model
//   SRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem);
  
//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile("data/models/SRPDE/2D_test5/y.csv");
//   DMatrix<double> y = yFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -6.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV calibrator
//   GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();

//   // expected values for q + Tr[S]
//   std::vector<double> expected_edfs = {
//     392.6223660935990551, 366.9539394059336246, 334.2897566350799252, 296.6918970050114126, 257.1195078756908288,
//     218.2528178022757857, 181.8976825480112609, 149.0838748306272237, 120.3659604356742108,  95.9678116102015224, 
//      75.7957057139429935,  59.4825236929897869,  46.5013015724528103
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)  
//   std::vector<double> expected_gcvs = {
//     0.4401244042041300, 0.4129157240364682, 0.3816562739955814, 0.3494124806447785, 0.3189922357544511,
//     0.2922854246010580, 0.2703666651130021, 0.2539506267229450, 0.2437506624632450, 0.2405042252882551,
//     0.2449930879669074, 0.2586192765896194, 0.2850382613491984
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );  

//   // check optimal lambda
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[9][0]) );  
// }

// /* test 6
//    domain:       unit square [0,1] x [0,1]
//    sampling:     locations = nodes
//    penalization: costant coefficients PDE
//    covariates:   no
//    BC:           no
//    order FE:     1
//    GCV optimization: grid stochastic
//  */
// TEST(GCV_SRPDE, Test6_CostantCoefficientsPDE_NonParametric_GeostatisticalAtNodes_GridStochastic) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("unit_square_coarse");

//   // non unitary diffusion tensor
//   SMatrix<2> K;
//   K << 1,0,0,4;
//   auto L = Laplacian(K); // anisotropic diffusion
  
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   // define statistical model
//   SRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem);
  
//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile("data/models/SRPDE/2D_test5/y.csv"); // load file for unit square coarse!
//   DMatrix<double> y = yFile.toEigen();

//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -6.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV calibrator
//   std::size_t seed = 4564168;
//   GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();
  
//   // expected values for q + Tr[S]
//   std::vector<double> expected_edfs = {
//     392.6398736612755442, 367.0053593585058707, 334.3960328473687582, 296.8543902422412657, 257.2935036770556962,
//     218.3463104599331643, 181.8107241424921483, 148.7611248424874759, 119.8187842113097616,  95.2531545500184080,
//      74.9790711890389048,  58.6159568240670694,  45.6216944905244262
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)  
//   std::vector<double> expected_gcvs = {
//     0.4404431338340472, 0.4134898057270393, 0.3824176190452916, 0.3502006998245064, 0.3195967826322993, 
//     0.2925309384160517, 0.2701852786717965, 0.2533900079360062, 0.2429208445928721, 0.2395110106902127, 
//     0.2439010923516004, 0.2574484277620766, 0.2837714099478065
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );  

//   // check optimal lambda
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[9][0]) );  
// }

// /* test 7
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid exact
//  */
// TEST(GCV_SRPDE, Test7_NonCostantCoefficientsPDE_NonParametric_Areal_GridExact) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("quasi_circle");

//   // load PDE coefficients data
//   CSVReader<double> reader{};
//   CSVFile<double> diffFile; // diffusion tensor
//   diffFile = reader.parseFile("data/models/SRPDE/2D_test4/K.csv");
//   DMatrix<double> diffData = diffFile.toEigen();
//   CSVFile<double> adveFile; // transport vector
//   adveFile = reader.parseFile("data/models/SRPDE/2D_test4/b.csv");
//   DMatrix<double> adveData = adveFile.toEigen();
  
//   // define non-constant coefficients
//   SpaceVaryingDiffusion<2> diffCoeff;
//   diffCoeff.setData(diffData);
//   SpaceVaryingAdvection<2> adveCoeff;
//   adveCoeff.setData(adveData);

//   auto L = Laplacian(diffCoeff.asParameter()) + Gradient(adveCoeff.asParameter());
  
//   // load non-zero forcing term
//   CSVFile<double> forceFile; // transport vector
//   forceFile = reader.parseFile("data/models/SRPDE/2D_test4/force.csv");
//   DMatrix<double> u = forceFile.toEigen();
  
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  
//   // define statistical model
//   CSVReader<int> int_reader{};
//   CSVFile<int> arealFile; // incidence matrix for specification of subdomains
//   arealFile = int_reader.parseFile("data/models/SRPDE/2D_test4/incidence_matrix.csv");
//   DMatrix<int> areal = arealFile.toEigen();

//   double lambda = std::pow(0.1, 3);
//   SRPDE<decltype(problem), fdaPDE::models::Areal> model(problem);
//   model.set_spatial_locations(areal);
  
//   // load data from .csv files
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile("data/models/SRPDE/2D_test4/z.csv");
//   DMatrix<double> y = yFile.toEigen();
  
//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -6.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV calibrator
//   GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();
  
//   // expected values for q + Tr[S]
//   std::vector<double> expected_edfs = {
//     6.8256470417592974, 6.7027363460970051, 6.5065235219756685, 6.2115932502007984, 5.8013394965992671,
//     5.2785614057558456, 4.6679680479060641, 4.0086692884299344, 3.3453578380021134, 2.7250862599643653,
//     2.1926047404607063, 1.7772541551903058, 1.4821993609506263
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)  
//   std::vector<double> expected_gcvs = {
//     24.1288371239571298, 23.3331732915000316, 22.1169607470995366, 20.4274688418074035, 18.3836687947684716,
//     16.3298382592033455, 14.7125701729663341, 13.8561091376951602, 13.8455543587040868, 14.5546522388024950,
//     15.7134846701210957, 16.9881160092813985, 18.1037530114315466
//   };
//   for(std::size_t i = 0; i < expected_gcvs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );  
  
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[8][0]) );
// }

// /* test 8
//    domain:       quasicircular domain
//    sampling:     areal
//    penalization: non-costant coefficients PDE
//    covariates:   no
//    BC:           yes
//    order FE:     1
//    GCV optimization: grid stochastic
//  */
// TEST(GCV_SRPDE, Test8_NonCostantCoefficientsPDE_NonParametric_Areal_GridStochastic) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("quasi_circle");

//   // load PDE coefficients data
//   CSVReader<double> reader{};
//   CSVFile<double> diffFile; // diffusion tensor
//   diffFile = reader.parseFile("data/models/SRPDE/2D_test4/K.csv");
//   DMatrix<double> diffData = diffFile.toEigen();
//   CSVFile<double> adveFile; // transport vector
//   adveFile = reader.parseFile("data/models/SRPDE/2D_test4/b.csv");
//   DMatrix<double> adveData = adveFile.toEigen();
  
//   // define non-constant coefficients
//   SpaceVaryingDiffusion<2> diffCoeff;
//   diffCoeff.setData(diffData);
//   SpaceVaryingAdvection<2> adveCoeff;
//   adveCoeff.setData(adveData);

//   auto L = Laplacian(diffCoeff.asParameter()) + Gradient(adveCoeff.asParameter());
  
//   // load non-zero forcing term
//   CSVFile<double> forceFile; // transport vector
//   forceFile = reader.parseFile("data/models/SRPDE/2D_test4/force.csv");
//   DMatrix<double> u = forceFile.toEigen();
  
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  
//   // define statistical model
//   CSVReader<int> int_reader{};
//   CSVFile<int> arealFile; // incidence matrix for specification of subdomains
//   arealFile = int_reader.parseFile("data/models/SRPDE/2D_test4/incidence_matrix.csv");
//   DMatrix<int> areal = arealFile.toEigen();

//   double lambda = std::pow(0.1, 3);
//   SRPDE<decltype(problem), fdaPDE::models::Areal> model(problem);
//   model.set_spatial_locations(areal);
  
//   // load data from .csv files
//   CSVFile<double> yFile; // observation file
//   yFile = reader.parseFile("data/models/SRPDE/2D_test4/z.csv");
//   DMatrix<double> y = yFile.toEigen();
  
//   // set model data
//   BlockFrame<double, int> df;
//   df.insert(OBSERVATIONS_BLK, y);
//   model.setData(df);
//   model.init(); // init model

//   // define grid of lambda values
//   std::vector<SVector<1>> lambdas;
//   for(double x = -6.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
//   // define GCV calibrator
//   std::size_t seed = 438172;
//   GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
//   GridOptimizer<1> opt;

//   ScalarField<1, decltype(GCV)> obj(GCV);
//   opt.optimize(obj, lambdas); // optimize gcv field
//   SVector<1> best_lambda = opt.optimum();

//   std::vector<double> expected_edfs = {
//     6.831047684097598, 6.712363557176861, 6.523569437275610, 6.241165650575311, 5.850412256000765,
//     5.354051019972657, 4.772559035126907, 4.136996797446094, 3.484589099615463, 2.860439790874401,
//     2.313681605547809, 1.880681758930778, 1.569948378164815
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

//   // expected value of gcv(\lambda)  
//   std::vector<double> expected_gcvs = {
//     25.6960716321667775, 24.9212379291962769, 23.7278907064885978, 22.0506431657932183, 19.9866647378768612,
//     17.8620876450079429, 16.1266862014873773, 15.1260832893288590, 14.9640061955584951, 15.5220163816125858,
//     16.5359284471161168, 17.6814853738732261, 18.6935898436080734
//   };
//   for(std::size_t i = 0; i < expected_edfs.size(); ++i)
//     EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );  

//   // check optimal lambda is correct
//   EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[8][0]) );  
// }


// GCV test for SQRPDE
// -------------------

/* test 9
   domain:       unit square [0,1] x [0,1] (coarse)
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: grid exact
 */
TEST(GCV_SQRPDE, Test9_Laplacian_NonParametric_GeostatisticalAtNodes_GridExact) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  double alpha = 0.1; 
  SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SQRPDE/2D_test1/z.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -8.0; x <= -5.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV function and optimize
  GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
  GridOptimizer<1> opt;

  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();
  
  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    440.1795193868704, 440.4548360048558, 437.5031093759197, 429.1941462350546, 
    413.3066572959513, 386.2208699550068, 348.1542714973974, 300.2020494486881, 
    240.5712416385806, 189.5029669915334, 147.4448979828643, 108.8464304505313, 81.71407293458549
  };
  for(std::size_t i = 0; i < expected_edfs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );
  
  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.1843258585917735, 0.08717619210568393, 0.03371374125087186, 0.01427976715611628, 
    0.006597947357249371, 0.003597353032388866, 0.002330792403062457, 0.001694801815432476, 
    0.001592253180401067, 0.001861788915609857, 0.002347903955410432, 0.003298553296342958, 0.004354027656093313
  };
  for(std::size_t i = 0; i < expected_edfs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[8][0]) );
     
}


/* test 10
   domain:       unit square [0,1] x [0,1] (coarse)
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: grid stochastic
 */
TEST(GCV_SQRPDE, Test10_Laplacian_NonParametric_GeostatisticalAtNodes_GridStochastic) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  double alpha = 0.1; 
  SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SQRPDE/2D_test1/z.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -8.0; x <= -5.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV calibrator
  std::size_t seed = 438172;
  GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 1000, seed);
  GridOptimizer<1> opt;

  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();
  
  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    440.1833466448017, 440.4602065737022, 437.5351076472637, 429.2612919993422, 
    413.4118759555859, 386.3430407536638, 348.2603122745472, 300.2762033244442, 
    240.7352965858459, 189.7061087038693, 147.75784423317, 109.3592217525395, 82.24828527401776
  };
  for(std::size_t i = 0; i < expected_edfs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );
   
  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.1860575986995367, 0.08891950641690292, 0.03433930914854689, 0.01444359577098376, 
    0.006648371303191344, 0.003613452818438541, 0.002336125617948957, 0.001696588424925723, 
    0.001594862966352748, 0.001864800209314017, 0.002352917967170364, 0.003308761794494414, 0.004367004358329372
  };
  for(std::size_t i = 0; i < expected_edfs.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[8][0]) );
     
}


/* test 11
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   GCV optimization: grid exact
 */
TEST(GCV_SQRPDE, Test11_Laplacian_SemiParametric_GeostatisticalAtLocations_GridExact) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/SQRPDE/2D_test2/locs.csv");
  DMatrix<double> loc = locFile.toEigen();

  // define statistical model
  double alpha = 0.9;
  SQRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem, alpha);
  model.set_spatial_locations(loc);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/SQRPDE/2D_test2/z.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/SQRPDE/2D_test2/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK,  y);
  df.insert(DESIGN_MATRIX_BLK, X);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -5.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV calibrator
  GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
  GridOptimizer<1> opt;
  
  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();

  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    132.9609167899539, 115.6188521061802, 96.98364501070189, 80.85716309709031, 
    67.0846022368506, 55.63929214088965, 45.06689308847258, 36.37777081975081, 30.23149565775936
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  
 
  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.01098722004315636, 0.009488475187585613, 0.008989810196085437, 0.008962959206273264, 
    0.008730799662312284, 0.008491268873056851, 0.008233546988020136, 0.008276226659556974, 0.00844007718180541
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[6][0]) );
 
}

/* test 12
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   GCV optimization: grid stochastic
 */
TEST(GCV_SQRPDE, Test12_Laplacian_SemiParametric_GeostatisticalAtLocations_GridStochastic) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  CSVReader<double> reader{};
  // load locations where data are sampled
  CSVFile<double> locFile;
  locFile = reader.parseFile("data/models/SQRPDE/2D_test2/locs.csv");
  DMatrix<double> loc = locFile.toEigen();

  // define statistical model
  double alpha = 0.9;
  SQRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem, alpha);
  model.set_spatial_locations(loc);
  
  // load data from .csv files
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/SQRPDE/2D_test2/z.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/SQRPDE/2D_test2/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK,  y);
  df.insert(DESIGN_MATRIX_BLK, X);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -5.0; x <= -3.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10, x)));
  
  // define GCV calibrator
  std::size_t seed = 66546513;
  GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 1000, seed);
  GridOptimizer<1> opt;

  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();

  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    132.8436017261902, 115.4993621499878, 96.88454935643702, 80.68461366031806, 
    66.80508876832693, 55.33985331850447, 44.76154075634351, 36.23318157143425, 30.27186405859233
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  


  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.01094769035605072, 0.009461009789507948, 0.008972198355574178, 0.008936612854569186, 
    0.008693636994372356, 0.008455660510399613, 0.008200766341706319, 0.008261438477268663, 0.008444140345302014
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[6][0]) );

}


/* test 13
   domain:       unit square [0,1] x [0,1]
   sampling:     locations = nodes
   penalization: costant coefficients PDE
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: grid exact
 */
TEST(GCV_SQRPDE, Test13_ConstantCoefficientsPDE_NonParametric_GeostatisticalAtNodes_GridExact) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");

  // non unitary diffusion tensor
  SMatrix<2> K;
  K << 1,0,0,4;
  auto L = Laplacian(K); // anisotropic diffusion
  
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
double alpha = 0.1; 
SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SQRPDE/2D_test3/z.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -7.0; x <= -5.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV calibrator
  GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
  GridOptimizer<1> opt;

  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();

  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    276.0886362314304, 234.8633711874632, 198.1724603886276, 164.8759545543741, 
    135.94089902281, 106.7993570612985, 79.51151604434595, 49.91561596383931, 27.4680328065335
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.003040397358024496, 0.002580598024481436, 0.002330573603960109, 0.002194742107880386, 
    0.002119607936620425, 0.00215815254420581, 0.00312027986655731, 0.008227789543693067, 0.03548656981038163
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[4][0]) );


}


/* test 14
   domain:       unit square [0,1] x [0,1]
   sampling:     locations = nodes
   penalization: costant coefficients PDE
   covariates:   no
   BC:           no
   order FE:     1
   GCV optimization: grid stochastic
 */
TEST(GCV_SQRPDE, Test14_ConstantCoefficientsPDE_NonParametric_GeostatisticalAtNodes_GridStochastic) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square_coarse");

  // non unitary diffusion tensor
  SMatrix<2> K;
  K << 1,0,0,4;
  auto L = Laplacian(K); // anisotropic diffusion
  
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  double alpha = 0.1; 
  SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile("data/models/SQRPDE/2D_test3/z.csv");
  DMatrix<double> y = yFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK, y);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -7.0; x <= -5.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV calibrator
  std::size_t seed = 438172;
  GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 1000, seed);
  GridOptimizer<1> opt;
  
  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();

  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    276.6602875384942, 235.4466919705868, 198.6662820414302, 165.3284639426754, 
    136.2693101328658, 107.0030730475357, 79.44094556923778, 50.04651002203742, 27.10440207647036
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.003061586028406282, 0.002595265288894978, 0.002340081652151136, 0.002201953273269074, 
    0.002124179040987324, 0.002160786007074692, 0.003119061929181641, 0.008233299912979723, 0.03542424327730619
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[4][0]) );

}


/* test 15
   domain:       c-shaped
   sampling:     areal
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   GCV optimization: grid exact
 */
TEST(GCV_SQRPDE, Test15_Laplacian_SemiParametric_Areal_GridExact) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped_areal");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  double alpha = 0.5;

  CSVReader<int> int_reader{};
  CSVFile<int> arealFile; // incidence matrix for specification of subdomains
  arealFile = int_reader.parseFile("data/models/SQRPDE/2D_test4/incidence_matrix.csv");
  DMatrix<int> areal = arealFile.toEigen();

  SQRPDE<decltype(problem), fdaPDE::models::Areal> model(problem, alpha);
  model.set_spatial_locations(areal);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/SQRPDE/2D_test4/z.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/SQRPDE/2D_test4/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK,  y);
  df.insert(DESIGN_MATRIX_BLK, X);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -4.0; x <= -1.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV calibrator
  GCV<decltype(model), ExactEDF<decltype(model)>> GCV(model);
  GridOptimizer<1> opt;
  
  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();

  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    15.5665805247258, 14.57206645250153, 12.07731267732268, 10.81888046946576, 
    9.132700712396884, 8.905677185957018, 7.543159447321961, 6.728320474488777, 
    6.102418202565106, 4.577356496044276, 4.056257902018482, 3.456352438552989, 2.524958204290085
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.07292942774413834, 0.06259585948366753, 0.04253643967049681, 0.03840749538881675, 
    0.03146488693048492, 0.03196533654328652, 0.02750310175047853, 0.02637187613955527, 
    0.03089986126232853, 0.04728237874711664, 0.09558382319047626, 0.2586941727343947, 0.5440417366341342
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[7][0]) );

}

/* test 16
   domain:       c-shaped
   sampling:     areal
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
   GCV optimization: grid stochastic
 */
TEST(GCV_SQRPDE, Test16_Laplacian_SemiParametric_Areal_GridStochastic) {
  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("c_shaped_areal");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  PDE problem(domain.mesh, L, u); // definition of regularizing PDE

  // define statistical model
  double alpha = 0.5;

  CSVReader<int> int_reader{};
  CSVFile<int> arealFile; // incidence matrix for specification of subdomains
  arealFile = int_reader.parseFile("data/models/SQRPDE/2D_test4/incidence_matrix.csv");
  DMatrix<int> areal = arealFile.toEigen();

  SQRPDE<decltype(problem), fdaPDE::models::Areal> model(problem, alpha);
  model.set_spatial_locations(areal);
  
  // load data from .csv files
  CSVReader<double> reader{};
  CSVFile<double> yFile; // observation file
  yFile = reader.parseFile  ("data/models/SQRPDE/2D_test4/z.csv");
  DMatrix<double> y = yFile.toEigen();
  CSVFile<double> XFile; // design matrix
  XFile = reader.parseFile  ("data/models/SQRPDE/2D_test4/X.csv");
  DMatrix<double> X = XFile.toEigen();

  // set model data
  BlockFrame<double, int> df;
  df.insert(OBSERVATIONS_BLK,  y);
  df.insert(DESIGN_MATRIX_BLK, X);
  model.setData(df);
  model.init(); // init model

  // define grid of lambda values
  std::vector<SVector<1>> lambdas;
  for(double x = -4.0; x <= -1.0; x +=0.25) lambdas.push_back(SVector<1>(std::pow(10,x)));
  
  // define GCV calibrator
  std::size_t seed = 438172;
  GCV<decltype(model), StochasticEDF<decltype(model)>> GCV(model, 100, seed);
  GridOptimizer<1> opt;

  ScalarField<1, decltype(GCV)> obj(GCV);
  opt.optimize(obj, lambdas); // optimize gcv field
  SVector<1> best_lambda = opt.optimum();

  // expected values for q + Tr[S]
  std::vector<double> expected_edfs = {
    15.67597974463673, 14.72866022659128, 12.33514848654242, 11.15275933010177, 
    9.444495374763706, 9.149456717766643, 7.443647450043284, 7.096969631945466, 
    6.638288840960764, 5.167892242081829, 4.0219954504947, 4.698181394281461, 3.348548069419884
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_edfs[i], GCV.edfs()[i]) );  

  // expected value of gcv(\lambda)
  std::vector<double> expected_gcvs = {
    0.07666639032156336, 0.06637012422164439, 0.0454463152958035, 0.04136105304760033, 
    0.03335119734921556, 0.03341780425789594, 0.02706889231931407, 0.02790033132991809, 
    0.03342803502267503, 0.05112239755457207, 0.09517433215647723, 0.3023869790965932, 0.5991898507520871
  };
  for(std::size_t i = 0; i < lambdas.size(); ++i)
    EXPECT_TRUE( almost_equal(expected_gcvs[i], GCV.values()[i]) );

  // check optimal lambda
  EXPECT_TRUE( almost_equal(best_lambda[0], lambdas[6][0]) );

}


