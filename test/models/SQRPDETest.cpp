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
#include "../fdaPDE/models/regression/SQRPDE.h"
using fdaPDE::models::SQRPDE;
#include "../fdaPDE/models/SamplingDesign.h"
#include "../../fdaPDE/models/regression/Distributions.h"

#include "../utils/MeshLoader.h"
using fdaPDE::testing::MeshLoader;
#include "../utils/Constants.h"
using fdaPDE::testing::DOUBLE_TOLERANCE;
#include "../utils/Utils.h"
using fdaPDE::testing::almost_equal;

#include<fstream>
#include<iostream>

#include <chrono>
using namespace std::chrono;


/* test 1
   domain:       unit square
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1
   Data generation: as described in data/models/SRPDE/test1/README.md 

 */

TEST(SQRPDE, Test1_Laplacian_NonParametric_GeostatisticalAtNodes) {

  
  // Parameters 
  const std::string TestNumber = "1"; 
  std::string data_macro_strategy_type = "matern_data"; 
  std::string data_strategy_type = "F"; 

  
  // Marco
  std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  // Ilenia 
  // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 
  
  double tol_weights = 0.000001; 
  double tol_FPIRLS = 0.000001;


  // define domain and regularizing PDE
  MeshLoader<Mesh2D<>> domain("unit_square");
  auto L = Laplacian();
  DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);

  
  bool massLumping_system = false;
  bool massLumping_GCV = false; 

  PDE problem(domain.mesh, L, u, massLumping_system); // definition of regularizing PDE

  unsigned int M = 10; 
  std::vector<double> alphas = {0.05, 0.25, 0.75, 0.95}; 

  for(double alpha : alphas){

    std::cout << "-----------------------------------------------------------" << std::endl; 
    std::cout << "alpha =  " << alpha << std::endl;

    unsigned int alpha_int = alpha*100; 
    const std::string alpha_string = std::to_string(alpha_int); 
    
    // define the statistical model 
    SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);

    for(unsigned int m = 1; m <= M; ++m){

      std::string path_solutions = R_path + "/R/Our/data/Test_" + 
                      TestNumber + "/alpha_" + alpha_string + "/" + data_macro_strategy_type + "/strategy_"  + 
                      data_strategy_type + "/our/sim_" + std::to_string(m);   
                  
      // load data from .csv files
      CSVReader<double> reader{};
      CSVFile<double> yFile; // observation file
      yFile = reader.parseFile(path_solutions + "/z.csv");             
      DMatrix<double> y = yFile.toEigen();


      // set model data
      BlockFrame<double, int> df;
      df.insert(OBSERVATIONS_BLK,  y);
      model.setData(df);

      CSVFile<double> lambdaCSV; 

      // Use optimal lambda to avoid possible numerical issues
      double lambda;               // read from C++
      // read from C++
      std::ifstream fileLambda(path_solutions + "/LambdaCpp.csv");
      if (fileLambda.is_open()){
        fileLambda >> lambda; 
        fileLambda.close();
      }
      model.setLambdaS(lambda);    //read from C++

      // Set
      model.setMassLumpingGCV(massLumping_GCV); 

      // solve smoothing problem
      auto t0 = high_resolution_clock::now();
      model.init();     
      model.solve();
      auto t1 = high_resolution_clock::now();
      std::chrono::duration<double> delta_time = t1 - t0;


      // Save C++ solution 
      DMatrix<double> computedF = model.f();
      const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
      std::ofstream filef(path_solutions + "/fCpp.csv");
      if (filef.is_open()){
        filef << computedF.format(CSVFormatf);
        filef.close();
      }


      DMatrix<double> computedFn = model.Psi()*model.f();
      const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
      std::ofstream filefn(path_solutions + "/fnCpp.csv");
      if (filefn.is_open()){
        filefn << computedFn.format(CSVFormatfn);
        filefn.close();
      }

      // // Duration 
      // std::ofstream myfileTime(path_solutions + "/Time_Cpp.csv");
      // if (myfileTime.is_open()){
      //   myfileTime << std::setprecision(16) << delta_time.count() << "\n";
      //   myfileTime.close();
      // }


    }
  }

} 


// TEST(SQRPDE, Test_massLumping_Laplacian_NonParametric_GeostatisticalAtLocations) {

  
//   // Parameters 
//   const std::string TestNumber = "lumping"; 
//   double alpha = 0.5; 
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int); 

//   // Marco
//   std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 
  
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("unit_square_32");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);

//   bool massLumping_system = false;
//   bool massLumping_GCV = false; 
//   std::string mass_type; 
//   if(!massLumping_system & !massLumping_GCV)
//     mass_type = "FALSE";
//   if(massLumping_system & massLumping_GCV)
//     mass_type = "TRUE";

//   PDE problem(domain.mesh, L, u, massLumping_system); // definition of regularizing PDE
//   SQRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem, alpha);

//   unsigned int M = 10;
//   std::vector<double> lambdas;
//   for(double x = -7.0; x <= -6.72; x +=0.01) lambdas.push_back(std::pow(10,x)); 

//   DMatrix<double> RMSE_mat; 
//   RMSE_mat.conservativeResize(M, lambdas.size());

//   unsigned int count_lambda = 0; 
  
//   for(auto lambda : lambdas){

//     count_lambda ++; 
//     model.setLambdaS(lambda); 
//     DVector<double> RMSE_vector(M); 

//     for(unsigned int m = 1; m <= M; ++m){

//       std::cout << "------------------------------------------- " << std::endl; 
//       std::cout << "Simulation " << m << std::endl;     

//       std::string path_solutions = R_path + "/R/Our/data/Test_" + 
//                                   TestNumber + "/alpha_" + alpha_string + "/lump" + 
//                                   mass_type + "/sim_" + std::to_string(m);  
                              
//       // load data from .csv files
//       CSVReader<double> reader{};
//       CSVFile<double> yFile; // observation file
//       yFile = reader.parseFile(path_solutions + "/z.csv");             
//       DMatrix<double> y = yFile.toEigen();
//       CSVFile<double> fTrueFile; // true solution file
//       fTrueFile = reader.parseFile(path_solutions + "/f_true.csv");             
//       DMatrix<double> fTrue = fTrueFile.toEigen();

//       // load locations where data are sampled
//       CSVFile<double> locFile;
//       locFile = reader.parseFile(R_path + "/R/Our/data/Test_" + TestNumber + "/alpha_" + alpha_string + 
//                                 "/locs.csv");
//       DMatrix<double> loc = locFile.toEigen();
//       model.set_spatial_locations(loc);

//       // set model data
//       BlockFrame<double, int> df;
//       df.insert(OBSERVATIONS_BLK,  y);
//       model.setData(df);

//       // Setter
//       model.setMassLumpingGCV(massLumping_GCV); 

//       // solve smoothing problem
//       auto t0 = high_resolution_clock::now();
//       model.init();     
//       model.solve();
//       auto t1 = high_resolution_clock::now();
//       std::chrono::duration<double> delta_time = t1 - t0;
      
//       DMatrix<double> computedF = model.f();
//       DMatrix<double> computedFn = model.Psi()*model.f();

//       // RMSE
//       double RMSE = 0.; 
//       for(auto i = 0; i < computedF.size(); ++i){
//         RMSE += std::pow(computedF(i,0) - fTrue(i,0), 2);   
//       }
//       RMSE = sqrt(RMSE / computedF.size()); 
//       RMSE_vector[m-1] = RMSE; 

//       // // Duration 
//       // std::ofstream myfileTime(path_solutions + "/Time_Cpp.csv");
//       // if (myfileTime.is_open()){
//       //   myfileTime << std::setprecision(16) << delta_time.count() << "\n";
//       //   myfileTime.close();
//       // }

//       // f 
//       const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filef(path_solutions + "/fCpp.csv");
//       if (filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//       }

//       // // fn 
//       // const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       // std::ofstream filefn(path_solutions + "/fnCpp.csv");
//       // if (filefn.is_open()){
//       //   filefn << computedFn.format(CSVFormatfn);
//       //   filefn.close();
//       // }

//     }  

//     RMSE_mat.col(count_lambda-1) = RMSE_vector; 

//   }

//   // Save C++ RMSE 
//   const static Eigen::IOFormat CSVFormatRMSE(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//   std::ofstream fileRMSE(R_path + "/R/Our/data/Test_" + 
//                     TestNumber + "/alpha_" + alpha_string + "/lump" + mass_type + "/RMSE_Cpp.csv");
//   if (fileRMSE.is_open()){
//     fileRMSE << RMSE_mat.format(CSVFormatRMSE);
//     fileRMSE.close();
//   }


// } 



/* test 2
   domain:       c-shaped
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
 */


// TEST(SQRPDE, Test2_Laplacian_SemiParametric_GeostatisticalAtLocations) {

//   const std::string TestNumber = "2";
//   double alpha = 0.5;  
//   std::string alpha_string = "50";

//   double tol_weights = 1e-6;
//   std::string tol_weights_string = "1e-6";
//   double tol_FPIRLS = 1e-6;
   
//   // Marco
//   std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 

//   // Ilenia 
//   // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 

//   // std::vector<int> seq_n = {4729, 5502};  // {250, 497, 997, 2003, 3983} ; 
//   // std::vector<std::string> seq_n_string = {"4729", "5502"};  // {"250", "497", "997", "2003", "3983"};

//   std::vector<int> seq_N = {597, 1020, 1520, 2003, 2563, 3040, 3503, 4012};  
//   std::vector<std::string> seq_N_string = {"597", "1020", "1520", "2003", "2563", "3040", "3503", "4012"};
  
//   unsigned int M = 1;

//   CSVFile<double> yFile; 
//   CSVFile<double> XFile;
//   DMatrix<double> y;  
//   DMatrix<double> X;
//   CSVFile<double> lambdaCSV; 


//   for(int m = 1; m <= M ; ++m){

//     std::cout << "Simulation m: " << m << std::endl; 
//     CSVReader<double> reader{};

//     //for(int n = 0; n < seq_n.size(); ++n){
//       for(int k = 0; k < seq_N_string.size(); ++k){

//         auto t0 = high_resolution_clock::now();

//         // define domain and regularizing PDE
//         MeshLoader<Mesh2D<>> domain("multiple_c_shaped/mesh_" +  seq_N_string[k]);  
//         auto L = Laplacian();
//         DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//         PDE problem(domain.mesh, L, u); // definition of regularizing PDE
//         // define statistical model
//         SQRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem, alpha);

//         auto t1 = high_resolution_clock::now();
//         std::chrono::duration<double> delta_model = t1 - t0;

//         // load locations where data are sampled
//         CSVFile<double> locFile;
//         locFile = reader.parseFile(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/locs.csv");  // "/sim_M/n_" + seq_n_string[n] + "/sim_" + std::to_string(m) + "/locs.csv");
//         DMatrix<double> loc = locFile.toEigen();



//         DMatrix<double> lambda;  // from R
//         //double lambda;   // from Cpp

//         // Read lambda from R
//         lambdaCSV = reader.parseFile(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/N_" + seq_N_string[k] + "/LambdaR" + ".csv");      
//         lambda = lambdaCSV.toEigen();
        
//         // // Read from Cpp
//         // std::ifstream fileLambda(R_path + "/R/Our/data/Test_" 
//         //                   + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//         //                   "/compare_nodes" + "/N_" + seq_N_string[k] + "/LambdaCpp_" + alpha_string + ".csv");  
//         // if(fileLambda.is_open()){
//         //   fileLambda >> lambda; 
//         //   fileLambda.close();
//         // }
         

//         // load data from .csv files
//         yFile = reader.parseFile(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/z.csv");   // "/sim_M/n_" + seq_n_string[n] + "/sim_" + std::to_string(m) + "/z.csv"); 
//         y = yFile.toEigen(); 
//         XFile = reader.parseFile(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/X.csv"); 
//         X = XFile.toEigen();


//         // Set model data

//         auto t2 = high_resolution_clock::now();

//         model.setLambdaS(lambda(0,0));  // read from R 
//         //model.setLambdaS(lambda);    // read from C++

//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK, y);
//         df.insert(DESIGN_MATRIX_BLK, X);
        
//         model.set_spatial_locations(loc);
//         model.setData(df);
//         model.setTolerances(tol_weights, tol_FPIRLS); 

//         // Solve smoothing problem
//         model.init();
//         model.solve();

//         auto t3 = high_resolution_clock::now();
//         std::chrono::duration<double> delta_fit = t3 - t2;

//         std::cout << "Duration: "  << delta_model.count() + delta_fit.count() << "seconds" << std::endl;

//         // Save time 
//         std::ofstream myfile(R_path + "/R/Our/data/Test_" 
//                               + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                               "/compare_nodes" + "/N_" + seq_N_string[k] +
//                               "/Time_Cpp_optR.csv");
//         myfile << std::setprecision(16) << delta_model.count() + delta_fit.count() << "\n" ;


//         // Save C++ solution 
//         DMatrix<double> computedFn = model.Psi()*model.f();
//         const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filefn(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/N_" + seq_N_string[k] + "/fnCpp_optR_" + alpha_string + ".csv");

//         if (filefn.is_open()){
//         filefn << computedFn.format(CSVFormatfn);
//         filefn.close();
//         }

//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream filef(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/N_" + seq_N_string[k] +  "/fCpp_optR_" + alpha_string + ".csv");

//         if (filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//         }

//         DVector<double> computedBeta = model.beta();
//         const static Eigen::IOFormat CSVFormat_beta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::ofstream file_beta(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/N_" + seq_N_string[k] + "/betaCpp_optR_" + alpha_string + ".csv");
//         if (file_beta.is_open()){
//         file_beta << computedBeta.format(CSVFormat_beta);
//         file_beta.close();
//         }


//         double J = model.J_final_sqrpde();
//         std::ofstream fileJ(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/N_" + seq_N_string[k] + "/JCpp_optR.csv");
//         if (fileJ.is_open()){
//           fileJ << J;
//           fileJ.close();
//         }

//         std::size_t niter = model.niter_sqrpde();
//         std::ofstream filen(R_path + "/R/Our/data/Test_" 
//                           + TestNumber + "/alpha_" + alpha_string  + "/tol_weights_" + tol_weights_string + 
//                           "/compare_nodes" + "/N_" + seq_N_string[k] + "/niterCpp_optR.csv");
//         if (filen.is_open()){
//           filen << niter;
//           filen.close();
//         }
//       }


//     //}


//   }

// }




/* test 3
   domain:       unit square
   sampling:     locations = nodes
   penalization: costant coefficients PDE
   covariates:   no
   BC:           no
   order FE:     1

 */

// TEST(SQRPDE, Test3_Laplacian_NonParametric_GeostatisticalAtNodes) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("unit_square");
//   const std::string TestNumber = "3"; 

//   // non unitary diffusion tensor
//   SMatrix<2> K;
//   K << 1, 0., 0., 4;
//   auto L = Laplacian(K); // anisotropic diffusion

//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   double alpha = 0.9; 
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int); 
  
  
//   SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);

//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   std::string data_macro_strategy_type = "skewed_data"; 
//   std::string data_strategy_type = "E"; 

//   // Marco
//   std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  
//   // Ilenia 
//   // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 
  

//   std::vector<double> seq_tol_weights = {0.000001}; 
//   std::vector<std::string> seq_tol_weights_string = {"1e-06"}; 
//   std::vector<double> seq_tol_FPIRLS = {0.000001};  
//   std::vector<std::string> seq_tol_FPIRLS_string = {"1e-06"}; 

//   std::string lin_sys_solver = "Chol";    // depends on the "symmetry" option in R 
//   std::string stopping_type = "our";

//   bool readFromR = false; 
//   unsigned int M = 10;      // number of simulations
//   std::string GCV_type = "Exact"; 

// for(int m = 1; m <= M; ++m ){
//   for(int i = 0; i < seq_tol_weights.size(); ++i ){

//     std::string tol_weights_string = seq_tol_weights_string[i];
//     double tol_weights = seq_tol_weights[i]; 

//       for(int j = 0; j < seq_tol_FPIRLS.size(); ++j){

//         std::string tol_FPIRLS_string = seq_tol_FPIRLS_string[j]; 
//         double tol_FPIRLS = seq_tol_FPIRLS[j];  

//         std::string path_solutions = R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/" + data_macro_strategy_type + "/strategy_"  + data_strategy_type + 
//                 "/" + stopping_type + "/tol_weights_" + tol_weights_string + "/tol_FPIRLS_" + tol_FPIRLS_string + 
//                 "/" + lin_sys_solver + "/sim_" + std::to_string(m); 

//         std::string path_GCV = path_solutions + "/GCV/" + GCV_type;  

//         yFile = reader.parseFile(path_solutions  + "/z.csv");             
//         DMatrix<double> y = yFile.toEigen();


//         // set model data
//         BlockFrame<double, int> df;
//         df.insert(OBSERVATIONS_BLK,  y);  
//         model.setData(df);


//         CSVFile<double> lambdaCSV; 

//         if(readFromR){
//           // Read from R
//           DMatrix<double> lambda;
//           lambdaCSV = reader.parseFile(path_solutions + 
//                           "/LambdaR.csv");    
          
//           lambda = lambdaCSV.toEigen();
//           model.setLambdaS(lambda(0,0));
//         } else{
//           double lambda;  
//           // Read from Cpp
//           std::ifstream fileLambda(path_solutions +
//                           "/LambdaCpp.csv");  
//           if (fileLambda.is_open()){
//           fileLambda >> lambda; 
//           fileLambda.close();
//           }   
//           model.setLambdaS(lambda);                    
//         }
//         // solve smoothing problem
//         model.init();     
//         model.solve();

//         // Save C++ solution 
//         DMatrix<double> computedF = model.f();
//         const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//         std::string opt_string; 
//         if(readFromR){
//           opt_string = "_optR_"; 
//         } else{
//           opt_string = "_"; 
//         }
//         std::ofstream filef(path_solutions + 
//                         "/fnCpp.csv");

//         if (filef.is_open()){
//           filef << computedF.format(CSVFormatf);
//           filef.close();
//         }

//       }
//     }
// }


// }


/* test 4
   domain:       c-shaped
   sampling:     areal
   penalization: laplacian
   covariates:   yes
   BC:           no
   order FE:     1
 */
// TEST(SQRPDE, Test4_NonParametric_Areal) {
//  // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  
  
//   // Ilenia 
//   std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared";

//   double alpha = 0.1; 
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int);
//   const std::string TestNumber = "4"; 

//   std::string path_solutions = R_path + "/R/Our/data/Test_" + 
//           TestNumber + "/alpha_" + alpha_string + "/Test_GSRPDE_Palu/dati_sqrpde" ;
//   // define domain and regularizing PDE
//   MeshLoader<Mesh2D<>> domain("c_shaped_areal");
//   CSVReader<double> reader{};
  
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
  
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE
  
//   // define statistical model
//   CSVReader<int> int_reader{};
//   CSVFile<int> arealFile; // incidence matrix for specification of subdomains
//   arealFile = int_reader.parseFile(path_solutions + "/incidence_matrix.csv");
//   DMatrix<int> areal = arealFile.toEigen();

  
//   SQRPDE<decltype(problem), fdaPDE::models::Areal> model(problem, alpha);

//   unsigned int M = 1;

//   for(int m = 1; m <= M; ++m){

//     // Read from Cpp
//     double lambda;   
//     std::ifstream fileLambda(path_solutions + "/sim_" + std::to_string(m) + "/GCV/Exact/LambdaCpp.csv");
//     if (fileLambda.is_open()){
//       fileLambda >> lambda; 
//       fileLambda.close();
//     }


//     model.setLambdaS(lambda);
//     model.set_spatial_locations(areal);
    
//     // load data from .csv files
//     CSVFile<double> yFile; // observation file
//     yFile = reader.parseFile(path_solutions + "/sim_" + std::to_string(m) + "/z.csv");
//     DMatrix<double> y = yFile.toEigen();

//     CSVFile<double> XFile; // design matrix
//     XFile = reader.parseFile  (path_solutions + "/sim_" + std::to_string(m) + "/X.csv");
//     DMatrix<double> X = XFile.toEigen();

//     // set model data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     df.insert(DESIGN_MATRIX_BLK, X);
//     model.setData(df);
    
//     // solve smoothing problem
//     model.init();
//     model.solve();

//     // Save C++ solution 
//     // DMatrix<double> computedF = model.f();
//     // const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     // std::ofstream filef(path_solutions + "/sim_" + std::to_string(m) + "/fCpp.csv");
//     // if (filef.is_open()){
//     //   filef << computedF.format(CSVFormatf);
//     //   filef.close();
//     // }

//     // DMatrix<double> computedFn = model.Psi()*model.f();
//     // const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     // std::ofstream filefn(path_solutions + "/sim_" + std::to_string(m) + "/fnCpp.csv");
//     // if (filefn.is_open()){
//     //   filefn << computedFn.format(CSVFormatfn);
//     //   filefn.close();
//     // }

//     // DVector<double> computedBeta = model.beta();
//     // const static Eigen::IOFormat CSVFormat_beta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     // std::ofstream file_beta(path_solutions + "/sim_" + std::to_string(m) + "/betaCpp.csv");
//     // if (file_beta.is_open()){
//     //   file_beta << computedBeta.format(CSVFormat_beta);
//     //   file_beta.close();
//     // }

//   }


// }
//   // }


/* test 5
   domain:       unit square
   sampling:     locations  != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1

 */

// TEST(SQRPDE, Test5_Laplacian_SemiParametric_GeostatisticalAtLocations) {

//   double alpha = 0.5; 
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int);
//   const std::string TestNumber = "5"; 

//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   CSVFile<double> XFile; // covariates file
//   CSVFile<double> locFile; // locations file

//   std::string data_macro_strategy_type = "skewed_data"; 
//   std::string data_strategy_type = "E"; 

//   // Marco
//   std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  
//   // Ilenia 
//   // std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 
  
//   double tol_weights = 0.000001; 
//   double tol_FPIRLS = 0.000001;      
//   std::string GCV_type = "Exact";  

//   std::string lin_sys_solver = "Cholesky";   // Woodbury Cholesky
//   std::string lin_sys_solver_abbrv = "Chol";  // Chol Wood

//   std::string GCV_lin_sys_solver = "Woodbury";     // to be chosen only when GCV_type = "Stochastic"

//   CSVFile<double> lambdaCSV; 
//   double lambda;                        // to read from Cpp

//   DMatrix<double> loc; 
//   DMatrix<double> X; 
//   DMatrix<double> y; 

//   bool massLumping_system = true;
//   bool massLumping_GCV = true; 
//   std::string mass_type; 
//   if(!massLumping_system & !massLumping_GCV)
//     mass_type = "FF";
//   if(!massLumping_system & massLumping_GCV)
//     mass_type = "FT"; 
//   if(massLumping_system & massLumping_GCV)
//     mass_type = "TT";

//   unsigned int launch_sim = 1; 

//   for(int nsim = launch_sim; nsim <= launch_sim; ++nsim){

//       auto t0 = high_resolution_clock::now();
      
//       // define domain and regularizing PDE
//       MeshLoader<Mesh2D<>> domain("unit_square_71");  
//       auto L = Laplacian();
//       DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);
//       PDE problem(domain.mesh, L, u, massLumping_system); // definition of regularizing PDE
//       SQRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem, alpha);

//       auto t1 = high_resolution_clock::now();
//       std::chrono::duration<double> delta_model = t1 - t0;


//       // Read data

//       std::string path_solutions = R_path + "/R/Our/data/Test_" + 
//           TestNumber + "/alpha_" + alpha_string + "/" + data_macro_strategy_type + 
//           "/for_slides_lumping&solvers" + "/system_solver_" + lin_sys_solver_abbrv + 
//           "/lump" + mass_type  + "/sim_" + std::to_string(nsim); 

//       std::string path_GCV = path_solutions + "/GCV/" + GCV_type; 

//       yFile = reader.parseFile(path_solutions + "/z.csv");             
//       y = yFile.toEigen();

//       XFile = reader.parseFile(path_solutions + "/X.csv");             
//       X = XFile.toEigen();

//       locFile = reader.parseFile(path_solutions + "/locs.csv");
//       loc = locFile.toEigen();
//       model.set_spatial_locations(loc);


//       //Use optimal lambda to avoid possible numerical issues
//       //Read from Cpp
      
//       if(GCV_type == "Stochastic"){
//         std::ifstream fileLambda(path_solutions + "/LambdaCpp_" + GCV_lin_sys_solver + ".csv");
//         if(fileLambda.is_open()){
//           fileLambda >> lambda; 
//           fileLambda.close();
//         }
//       }       
//       if(GCV_type == "Exact"){
//         std::ifstream fileLambda(path_solutions + "/LambdaCpp.csv");
//         if(fileLambda.is_open()){
//           fileLambda >> lambda; 
//           fileLambda.close();
//         }
//       }

//       auto t2 = high_resolution_clock::now();

//       // set model data, locations, lambdas and tolerances  
//       model.set_spatial_locations(loc);
//       BlockFrame<double, int> df;
//       df.insert(OBSERVATIONS_BLK, y);
//       df.insert(DESIGN_MATRIX_BLK, X);
//       model.setData(df);

//       model.setLambdaS(lambda);         // to read from Cpp
//       model.setTolerances(tol_weights, tol_FPIRLS); 
//       model.setMassLumpingGCV(massLumping_GCV); 
//       model.setLinearSystemType(lin_sys_solver); 


//       // Solve smoothing problem
//       model.init();     
//       model.solve();
//       auto t3 = high_resolution_clock::now();
//       std::chrono::duration<double> delta_fit = t3 - t2;

//       std::cout << "Duration: "  << delta_model.count() + delta_fit.count() << "seconds" << std::endl;

//       std::ofstream myfile(path_solutions +  "/Time_Cpp.csv");
//       myfile << std::setprecision(16) << delta_model.count() + delta_fit.count() << "\n";

//       // Save C++ solution 
//       DMatrix<double> computedF = model.f();
//       const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filef(path_solutions + "/fCpp.csv");
//       if (filef.is_open()){
//         filef << computedF.format(CSVFormatf);
//         filef.close();
//       }

//       DMatrix<double> computedFn = model.Psi()*model.f();
//       const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream filefn(path_solutions + "/fnCpp.csv");
//       if (filefn.is_open()){
//         filefn << computedFn.format(CSVFormatfn);
//         filefn.close();
//       }

//       DVector<double> computedBeta = model.beta();
//       const static Eigen::IOFormat CSVFormat_beta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//       std::ofstream file_beta(path_solutions + "/betaCpp.csv");
//       if (file_beta.is_open()){
//         file_beta << computedBeta.format(CSVFormat_beta);
//         file_beta.close();
//       }

//       double J = model.J_final_sqrpde();
//       std::ofstream fileJ(path_solutions + "/JCpp.csv");
//       if (fileJ.is_open()){
//         fileJ << J;
//         fileJ.close();
//       }

//       std::size_t niter = model.niter_sqrpde();
//       std::ofstream filen(path_solutions + "/niterCpp.csv");
//       if (filen.is_open()){
//         filen << niter;
//         filen.close();
//       }

//   }
  
// }





/* test 7
   domain:       unit sphere
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
 */

// TEST(SQRPDE, Test7_Laplacian_NonParametric_GeostatisticalAtLocations) {
//   // define domain and regularizing PDE
//   MeshLoader<Mesh3D<>> domain("unit_sphere");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*6, 1);   // *6 ---> VERIFICARE 
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   double alpha = 0.1; 
//   const std::string alpha_string = "10"; 
//   const std::string TestNumber = "7"; 

//   SQRPDE<decltype(problem), fdaPDE::models::GeoStatLocations> model(problem, alpha);

//   // Marco
//   std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  
//   // Ilenia 
//   //std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 
  

//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; // observation file
//   CSVFile<double> XFile; 
//   CSVFile<double> locFile;
//   CSVFile<double> lambdaCSV;
//   DMatrix<double> y;  
//   DMatrix<double> X;
//   DMatrix<double> loc;
//   double lambda; 

//   unsigned int M = 30;

//   // load locations where data are sampled
//   locFile = reader.parseFile(R_path + "/R/Our/data/Test_" 
//                       + TestNumber + "/alpha_" + alpha_string + "/locs.csv");

//   loc = locFile.toEigen();
//   model.set_spatial_locations(loc); 

//   double tol_weights = 0.000001; 
//   std::string tol_weights_string = "1e-06"; 

//   double tol_FPIRLS = 0.000001;
//   std::string tol_FPIRLS_string = "1e-06"; 
//   model.setTolerances(tol_weights, tol_FPIRLS); 

//   for(std::size_t m=1; m<=M; m++) {

//     yFile = reader.parseFile(R_path + "/R/Our/data/Test_" + 
//                       TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/z.csv");             
//     y = yFile.toEigen();

//     XFile = reader.parseFile(R_path + "/R/Our/data/Test_" + 
//                     TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/X.csv");             
//     X = XFile.toEigen();

//     // set model data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     df.insert(DESIGN_MATRIX_BLK, X);

//     model.setData(df);

//     // Read from Cpp
//     std::ifstream fileLambda(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/LambdaCpp.csv");
//     if (fileLambda.is_open()){
//       fileLambda >> lambda; 
//       fileLambda.close();
//     }

//     model.setLambdaS(lambda);

//     // solve smoothing problem
//     model.init();     
//     model.solve();


//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filef(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/fCpp.csv");

//     if (filef.is_open()){
//       filef << computedF.format(CSVFormatf);
//       filef.close();
//     }

//     DMatrix<double> computedFn = model.Psi()*model.f();
//     const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filefn(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/fnCpp.csv");

//     if (filefn.is_open()){
//       filefn << computedFn.format(CSVFormatfn);
//       filefn.close();
//     }

//     DVector<double> computedBeta = model.beta();
//     const static Eigen::IOFormat CSVFormat_beta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream file_beta(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/betaCpp.csv");
//     if (file_beta.is_open()){
//       file_beta << computedBeta.format(CSVFormat_beta);
//       file_beta.close();
//     }


//     double J = model.J_final_sqrpde();
//     std::ofstream fileJ(R_path + "/R/Our/data/Test_" 
//                 + TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/JCpp.csv");
//     if (fileJ.is_open()){
//       fileJ << J;
//       fileJ.close();
//     }

//     std::size_t niter = model.niter_sqrpde();
//     std::ofstream filen(R_path + "/R/Our/data/Test_" 
//                 + TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/niterCpp.csv");
//     if (filen.is_open()){
//       filen << niter;
//       filen.close();
//     }

//   }


// }


/* test 8
   domain:       Hub
   sampling:     locations = nodes
   penalization: simple laplacian
   covariates:   yes
   BC:           no
   order FE:     1
 */

// TEST(SQRPDE, Test8_Laplacian_SemiParametric_GeostatisticalAtNodes) {
//   MeshLoader<SurfaceMesh<>> domain("surface_fine");

//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*6, 1);
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE


//   double alpha = 0.1; 
//   const std::string alpha_string = "10"; 
//   const std::string TestNumber = "8"; 

//   SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);

//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; 
//   CSVFile<double> XFile;
//   DMatrix<double> y;  
//   DMatrix<double> X;
//   double lambda; 

//   // Marco
//   // std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
  
//   // Ilenia 
//   std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared"; 
  
//   double tol_weights = 0.000001; 
//   std::string tol_weights_string = "1e-06";

//   double tol_FPIRLS = 0.000001; 
//   std::string tol_FPIRLS_string = "1e-06"; 


//   unsigned int M = 10;

//   for(std::size_t m=1; m<=M; m++) {

//     yFile = reader.parseFile(R_path + "/R/Our/data/Test_" + 
//                     TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/z.csv");             
//     y = yFile.toEigen();

//     XFile = reader.parseFile(R_path + "/R/Our/data/Test_" + 
//                     TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/X.csv");             
//     X = XFile.toEigen();

//     // set model data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);
//     df.insert(DESIGN_MATRIX_BLK, X);

//     model.setData(df);

//     // Read from Cpp
//     std::ifstream fileLambda(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/LambdaCpp.csv");
//     if (fileLambda.is_open()){
//       fileLambda >> lambda; 
//       fileLambda.close();
//     }

//     model.setLambdaS(lambda);

//     // solve smoothing problem
//     model.init();     
//     model.solve();


//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filef(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/fCpp.csv");

//     if (filef.is_open()){
//       filef << computedF.format(CSVFormatf);
//       filef.close();
//     }


//     DVector<double> computedBeta = model.beta();
//     const static Eigen::IOFormat CSVFormat_beta(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream file_beta(R_path + "/R/Our/data/Test_" + 
//                 TestNumber + "/alpha_" + alpha_string + "/sim_" + std::to_string(m) + "/betaCpp.csv");
//     if (file_beta.is_open()){
//       file_beta << computedBeta.format(CSVFormat_beta);
//       file_beta.close();
//     }


//   }


// }


/* test 9
   domain:       linear network
   sampling:     locations != nodes
   penalization: simple laplacian
   covariates:   no
   BC:           no
   order FE:     1

 */

// TEST(SQRPDE, Test9_Laplacian_NonParametric_GeostatisticalAtNodes) {

//   const std::string TestNumber = "9"; 

//   // define domain and regularizing PDE
//   MeshLoader<NetworkMesh<>> domain("network");
//   auto L = Laplacian();
//   DMatrix<double> u = DMatrix<double>::Zero(domain.mesh.elements()*3, 1);  
//   PDE problem(domain.mesh, L, u); // definition of regularizing PDE

//   double alpha = 0.9; 
//   unsigned int alpha_int = alpha*100; 
//   const std::string alpha_string = std::to_string(alpha_int);
//   SQRPDE<decltype(problem), fdaPDE::models::GeoStatMeshNodes> model(problem, alpha);

//   // Marco
//   std::string R_path = "/mnt/c/Users/marco/OneDrive - Politecnico di Milano/Corsi/Magistrale/Anno_II_Semestre_II/PACS_project_shared"; 
//   // Ilenia 
//   //std::string R_path = "/mnt/c/Users/ileni/OneDrive - Politecnico di Milano/PACS_project_shared";   

//   // load data from .csv files
//   CSVReader<double> reader{};
//   CSVFile<double> yFile; 
//   CSVFile<double> lambdaCSV;
//   DMatrix<double> y;  
//   double lambda; 

//   std::string GCV_type = "Exact"; 

//   unsigned int M = 10;
//   for(std::size_t m = 1; m <= M; m++) {

//     std::string path_solutions = R_path + "/R/Our/data/Test_" + 
//                     TestNumber + "/alpha_" + alpha_string + "/c_network/sim_" + std::to_string(m); 
//     std::string path_GCV =  path_solutions + "/GCV/" + GCV_type;   

//     yFile = reader.parseFile(path_solutions + "/z.csv");             
//     y = yFile.toEigen();

//     // set model data
//     BlockFrame<double, int> df;
//     df.insert(OBSERVATIONS_BLK, y);

//     model.setData(df);

//     // Read from Cpp
//     std::ifstream fileLambda(path_solutions + "/LambdaCpp.csv");
//     if (fileLambda.is_open()){
//       fileLambda >> lambda; 
//       fileLambda.close();
//     }

//     model.setLambdaS(lambda);

//     // solve smoothing problem
//     model.init();     
//     model.solve();


//     // Save C++ solution 
//     DMatrix<double> computedF = model.f();
//     const static Eigen::IOFormat CSVFormatf(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filef(path_solutions + "/fCpp.csv");
//     if (filef.is_open()){
//       filef << computedF.format(CSVFormatf);
//       filef.close();
//     }

//     DMatrix<double> computedFn = model.Psi()*model.f();
//     const static Eigen::IOFormat CSVFormatfn(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
//     std::ofstream filefn(path_solutions + "/fnCpp.csv");
//     if (filefn.is_open()){
//       filefn << computedFn.format(CSVFormatfn);
//       filefn.close();
//     }


//   }


// }

// MANCA HUB

