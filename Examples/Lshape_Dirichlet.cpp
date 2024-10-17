#include <Eigen/Sparse>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "BoundaryMesh.hpp"
#include "L2ProjectionToP1.hpp"
#include "L2WithP1Basis.hpp"
#include "LshapeMesh.hpp"
#include "LshapeMeshAdaptive.hpp"
#include "NeumannTraceInterpolation.hpp"
#include "buildK.hpp"
#include "buildM.hpp"
#include "buildV.hpp"
#include "buildW.hpp"
#include "evaluateK.hpp"
#include "evaluateV.hpp"
#include "gaussQuadrature.h"
#include "geometry.hpp"
#include "matplotlibcpp.h"
#include "potentials.hpp"

namespace plt = matplotlibcpp;

int main(int, char **) {
  // Exponents for adaptive meshes
  int NExps = 5;
  // alpha = 1 corresponds to a uniform mesh refinement
  double alphamin = 1;
  double alphamax = 2;

  Eigen::VectorXd Exponents =
      Eigen::VectorXd::LinSpaced(NExps, alphamin, alphamax);

  std::cout << "Exponents = \n" << Exponents.transpose() << std::endl;

  // Number of refinement levels
  int N = 5;
  Eigen::MatrixXd hdata = Eigen::MatrixXd::Zero(N, NExps);
  Eigen::MatrixXd err_data = Eigen::MatrixXd::Zero(N, NExps);

  // Looping over refinement levels
  for (int k = 0; k < N; ++k) {
    // Looping over mesh exponent levels
    for (int i = 0; i < NExps; ++i) {
      printf("Refinement level: %d, exponent number: %d \r", k, i);
      fflush(stdout);
      BoundaryMesh mesh = createLshapeMeshAdaptive(k + 1, Exponents(i));

      // Assembling the Galerkin matrices
      // Single layer Galerkin matrix
      Eigen::MatrixXd V;
      computeV(V, mesh, 0);

      // Double layer Galerkin matrix
      Eigen::MatrixXd K;
      computeK(K, mesh, 0);

      // Mass matrices
      // Test space (rows) = P1, Trial space (columns) = P1
      Eigen::SparseMatrix<double> M11(mesh.numVertices(), mesh.numVertices());
      computeM11(M11, mesh);
      // Test space (rows) = P0, Trial space (columns) = P1
      Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
      computeM01(M01, mesh);

      // L2 projection of Dirichlet data
      Eigen::VectorXd g_projected = L2ProjectionToP1(u, mesh);

      // RHS
      Eigen::VectorXd rhs_projected = (0.5 * M01 + K) * g_projected;
      Eigen::VectorXd sol_projected = V.lu().solve(rhs_projected);
      Eigen::VectorXd sol_ex = NeumannTraceInterpolation(gradu, mesh);

      err_data(k, i) = std::sqrt((sol_ex - sol_projected).transpose() * V *
                                 (sol_ex - sol_projected));

      // get the mean meshwidth
      double htot = 0;
      for (int j = 0; j < mesh.numElements(); ++j) {
        // identify element's vertices
        int aidx = mesh.getElementVertex(j, 0);
        int bidx = mesh.getElementVertex(j, 1);
        const Eigen::Vector2d &a = mesh.getVertex(aidx);
        const Eigen::Vector2d &b = mesh.getVertex(bidx);

        htot += (b - a).squaredNorm();
      }
      hdata(k, i) = sqrt(htot / mesh.numElements());
    }
  }
  printf("Computed errors at all refinement and exponent levels \n");
  std::cout << "hvals = \n" << hdata << std::endl;
  std::cout << "errs = \n" << err_data << std::endl;

  // Least squares fit to the computed points
  // Least squares fit parameters
  Eigen::MatrixXd lsqParams(2, NExps);

  // Plotting the different errors
  for (int i = 0; i < NExps; ++i) {
    // System matrix for least squares
    Eigen::MatrixXd A(N, 2);
    // log of hdata
    Eigen::VectorXd logh =
        hdata.col(i).unaryExpr([](double x) { return std::log(x); });

    // Populating the system matrix
    A << logh, Eigen::VectorXd::Ones(N);

    // Initializing the SVD object
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Log of err_data
    Eigen::VectorXd logerr =
        err_data.col(i).unaryExpr([](double x) { return std::log(x); });

    // Using SVD to compute the least squares solution
    lsqParams.col(i) = svd.solve(logerr);

    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << "Exp: " << Exponents(i) << " rate: " << lsqParams(0, i);

    // loglog plot of the error
    plt::loglog(hdata.col(i), err_data.col(i), {{"label", ss.str()}});

    // Plotting the least squares fit
    // log of fitted errors
    Eigen::VectorXd logerrfit = A * lsqParams.col(i);
    // fitted errors
    Eigen::VectorXd errfit =
        logerrfit.unaryExpr([](double x) { return std::exp(x); });

    // Log log plot of the fitted error
    plt::loglog(hdata.col(i), errfit, {{"linestyle", "--"}, {"color", "gray"}});
  }

  plt::xlabel("meshwidth h");
  plt::ylabel("Error");
  plt::legend();

  plt::show();

  return 0;
}
