#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "BoundaryMesh.hpp"
#include "L2ProjectionToP1.hpp"
#include "L2WithP1Basis.hpp"
#include "NeumannTraceInterpolation.hpp"
#include "buildK.hpp"
#include "buildM.hpp"
#include "buildV.hpp"
#include "buildW.hpp"
#include "capacitorMesh.hpp"
#include "evaluateK.hpp"
#include "evaluateV.hpp"
#include "gaussQuadrature.h"
#include "geometry.hpp"
#include "potentials.hpp"

int main() {
  // Number of refinement stages
  int N = 10;

  // Vectors for holding meshwidths and errors
  Eigen::VectorXd hdata(N), err_interp_data(N), err_proj_data(N);

  std::cout << std::setw(10) << "Meshwidth" << std::setw(20)
            << "Err. interpolation" << std::setw(20) << "Err. projected"
            << std::setw(20) << "Err. ratio interp." << std::setw(20)
            << "Err. ratio proj." << std::endl;

  // Going over the refinement stages
  for (int k = 0; k < N; ++k) {
    // Capacitor mesh
    BoundaryMesh mesh = createCapacitorMesh(k + 1);

    // Assembling Discrete BIOs
    // Single layer
    Eigen::MatrixXd V;
    computeV(V, mesh, 0);

    // Double layer
    Eigen::MatrixXd K;
    computeK(K, mesh, 0);

    // Mass matrix
    // Trial space (rows) = P0, Test space (columns) = P1
    Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
    computeM01(M01, mesh);

    // L2 projection of Dirichlet data
    Eigen::VectorXd g_projected = L2ProjectionToP1(u2, mesh);

    // Nodal interpolation for Dirichlet data, nodal values as coefficient
    Eigen::VectorXd g_interpolated = Eigen::VectorXd::Zero(mesh.numVertices());
    for (int i = 0; i < mesh.numVertices(); ++i) {
      g_interpolated(i) = u2(mesh.getVertex(i));
    }

    // Constructing the rhs
    Eigen::VectorXd rhs_interpolated = (0.5 * M01 + K) * g_interpolated;
    Eigen::VectorXd rhs_projected = (0.5 * M01 + K) * g_projected;

    // Solving the linear system
    Eigen::VectorXd sol_interpolated = V.lu().solve(rhs_interpolated);
    Eigen::VectorXd sol_projected = V.lu().solve(rhs_projected);

    // Neumann trace via interpolation
    Eigen::VectorXd sol_ex = NeumannTraceInterpolation(gradu2, mesh);

    // Evaluation point for the potential
    Eigen::Vector2d eval_pt(0.2, 0.13);

    // Point evaluation of the single layer and double layer potentials
    Eigen::VectorXd sl_proj, sl_interp, dl_proj, dl_interp;
    evaluateV(sl_proj, mesh, sol_projected, eval_pt.transpose(), 0);
    evaluateV(sl_interp, mesh, sol_interpolated, eval_pt.transpose(), 0);
    evaluateK(dl_proj, mesh, g_projected, eval_pt.transpose(), 0);
    evaluateK(dl_interp, mesh, g_interpolated, eval_pt.transpose(), 0);

    // Representation formula
    double eval_proj = sl_proj(0) - dl_proj(0);
    double eval_interp = sl_interp(0) - dl_interp(0);

    double exact_potential = u2(eval_pt);

    // Computing error
    err_interp_data(k) = std::fabs(exact_potential - eval_interp);
    err_proj_data(k) = std::fabs(exact_potential - eval_proj);

    // Meshwidth
    hdata(k) = (mesh.getVertex(1) - mesh.getVertex(0)).norm();

    // Printing a table of errors
    if (k == 0) {
      std::cout << std::setw(10) << hdata(k) << std::setw(20)
                << err_interp_data(k) << std::setw(20) << err_proj_data(k)
                << std::endl;
    } else {
      std::cout << std::setw(10) << hdata(k) << std::setw(20)
                << err_interp_data(k) << std::setw(20) << err_proj_data(k)
                << std::setw(20) << err_interp_data(k - 1) / err_interp_data(k)
                << std::setw(20) << err_proj_data(k - 1) / err_proj_data(k)
                << std::endl;
    }
  }

  return 0;
}
