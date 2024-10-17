#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "BoundaryMesh.hpp"
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

int main(int, char **) {
  // Number of refinement stages
  int N = 10;

  // Vectors for holding meshwidths and errors
  Eigen::VectorXd hdata(N), err_direct_data(N), err_indirect_data(N);

  std::cout << std::setw(10) << "meshwidth" << std::setw(20) << "Err. indirect"
            << std::setw(20) << "Err. direct" << std::setw(20)
            << "Err. ratio indirect" << std::endl;

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
    // Test space (rows) = P0, Trial space (columns) = P1
    Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
    computeM01(M01, mesh);

    // Nodal interpolation for Dirichlet data
    Eigen::VectorXd g_interpolated(mesh.numVertices());
    for (int i = 0; i < mesh.numVertices(); ++i) {
      g_interpolated(i) = u1(mesh.getVertex(i));
    }

    // RHS of linear system
    Eigen::VectorXd rhs_indirect = M01 * g_interpolated;
    Eigen::VectorXd rhs_direct = (0.5 * M01 + K) * g_interpolated;

    // Solving the linear system
    Eigen::VectorXd sol_indirect = V.lu().solve(rhs_indirect);
    Eigen::VectorXd sol_direct = V.lu().solve(rhs_direct);

    // Exact neumann trace using interpolation
    Eigen::VectorXd sol_ex = NeumannTraceInterpolation(gradu1, mesh);

    // Evaluation point for the potential
    Eigen::Vector2d eval_pt(0.2, 0.13);

    // Point evaluation of the single layer and double layer potentials
    Eigen::VectorXd sl_direct, sl_indirect, dl_direct;
    evaluateV(sl_direct, mesh, sol_direct, eval_pt.transpose(), 0);
    evaluateV(sl_indirect, mesh, sol_indirect, eval_pt.transpose(), 0);
    evaluateK(dl_direct, mesh, g_interpolated, eval_pt.transpose(), 0);

    // Representation formula
    double eval_direct = sl_direct(0) - dl_direct(0);
    double eval_indirect = sl_indirect(0);

    double exact_potential = u1(eval_pt);

    // Computing error
    err_indirect_data(k) = std::fabs(exact_potential - eval_indirect);
    err_direct_data(k) = std::fabs(exact_potential - eval_direct);

    hdata(k) = (mesh.getVertex(1) - mesh.getVertex(0)).norm();
    if (k == 0) {
      std::cout << std::setw(10) << hdata(k) << std::setw(20)
                << err_indirect_data(k) << std::setw(20) << err_direct_data(k)
                << std::endl;
    } else {
      std::cout << std::setw(10) << hdata(k) << std::setw(20)
                << err_indirect_data(k) << std::setw(20) << err_direct_data(k)
                << std::setw(20)
                << err_indirect_data(k - 1) / err_indirect_data(k) << std::endl;
    }
  }

  return 0;
}
