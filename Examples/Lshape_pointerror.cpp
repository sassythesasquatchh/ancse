#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "BoundaryMesh.hpp"
#include "L2ProjectionToP1.hpp"
#include "L2WithP1Basis.hpp"
#include "LshapeMesh.hpp"
#include "NeumannTraceInterpolation.hpp"
#include "buildK.hpp"
#include "buildM.hpp"
#include "buildV.hpp"
#include "buildW.hpp"
#include "evaluateK.hpp"
#include "evaluateV.hpp"
#include "gaussQuadrature.h"
#include "geometry.hpp"
#include "potentials.hpp"

int main(int, char **) {
  int N = 10;
  Eigen::VectorXd hdata(N), err_proj(N), err_interp(N);
  std::cout << std::setw(15) << "h" << std::setw(20) << "Error Proj."
            << std::setw(20) << "Error Interp." << std::setw(20)
            << "Err. ratio Proj." << std::setw(20) << "Err. ratio Interp."
            << std::endl;
  for (int k = 0; k < N; ++k) {
    // Creating the LShapeMesh
    BoundaryMesh mesh = createLshapeMesh(k + 1);

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

    // Nodal interpolation for Dirichlet data, nodal values as coefficient
    Eigen::VectorXd g_interpolated = Eigen::VectorXd::Zero(mesh.numVertices());

    for (int i = 0; i < mesh.numVertices(); ++i) {
      g_interpolated(i) = u(mesh.getVertex(i));
    }

    // RHS
    Eigen::VectorXd rhs_interpolated = (0.5 * M01 + K) * g_interpolated;
    Eigen::VectorXd rhs_projected = (0.5 * M01 + K) * g_projected;

    Eigen::VectorXd sol_interpolated = V.lu().solve(rhs_interpolated);
    Eigen::VectorXd sol_projected = V.lu().solve(rhs_projected);

    Eigen::VectorXd sol_ex = NeumannTraceInterpolation(gradu, mesh);

    Eigen::Vector2d eval_pt(0.2, 0.13);

    // Point evaluation of the potential
    Eigen::VectorXd sl_proj, sl_interp, dl_proj, dl_interp;

    // Computing single and double layers potentials
    evaluateV(sl_proj, mesh, sol_projected, eval_pt.transpose(), 0);
    evaluateV(sl_interp, mesh, sol_interpolated, eval_pt.transpose(), 0);
    evaluateK(dl_proj, mesh, g_projected, eval_pt.transpose(), 0);
    evaluateK(dl_interp, mesh, g_interpolated, eval_pt.transpose(), 0);

    double exact_potential = u(eval_pt);

    // Computing error with the representation formula
    err_proj(k) = std::fabs(exact_potential - (sl_proj(0) - dl_proj(0)));
    err_interp(k) = std::fabs(exact_potential - (sl_interp(0) - dl_interp(0)));

    hdata(k) = (mesh.getVertex(1) - mesh.getVertex(0)).norm();
    if (k == 0) {
      std::cout << std::setw(15) << hdata(k) << std::setw(20) << err_proj(k)
                << std::setw(20) << err_interp(k) << std::endl;
    } else {
      std::cout << std::setw(15) << hdata(k) << std::setw(20) << err_proj(k)
                << std::setw(20) << err_interp(k) << std::setw(20)
                << err_proj(k - 1) / err_proj(k) << std::setw(20)
                << err_interp(k - 1) / err_interp(k) << std::endl;
    }
  }

  return 0;
}
