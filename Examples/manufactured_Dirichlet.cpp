#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "BoundaryMesh.hpp"
#include "buildK.hpp"
#include "buildM.hpp"
#include "buildV.hpp"
#include "buildW.hpp"
#include "evaluateV.hpp"
#include "gaussQuadrature.h"
#include "geometry.hpp"

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Vector type
typedef Eigen::VectorXd Vector;

int main(int, char **) {

  BoundaryMesh mesh("capacitor");

  Eigen::MatrixXd V;
  computeV(V, mesh, 0);

  Eigen::MatrixXd K;
  computeK(K, mesh, 0);

  Eigen::SparseMatrix<double> M11(mesh.numVertices(), mesh.numVertices());
  computeM11(M11, mesh);

  Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
  computeM01(M01, mesh);

  auto u2 = [](Eigen::Vector2d x) { return std::cosh(x(0)) * std::cos(x(1)); };

  // L2 projection of Dirichlet data
  // Compute the integral of g * beta for all basis functions beta in the space
  // P1 Can be done by simply using Gauss quadrature
  Eigen::VectorXd g_beta = Eigen::VectorXd::Zero(mesh.numVertices());

  unsigned order = 32;
  // Get the Gauss quadrature points and weights for the interval [-1,1]
  const double *X = getGaussPoints(order);
  const double *W = getGaussWeights(order);

  // Local to global assembly. Start traversing the panels
  for (int i = 0; i < mesh.numElements(); i++) {
    // identify element's vertices
    int aidx = mesh.getElementVertex(i, 0);
    int bidx = mesh.getElementVertex(i, 1);
    const Eigen::Vector2d &a = mesh.getVertex(aidx);
    const Eigen::Vector2d &b = mesh.getVertex(bidx);

    double h = (b - a).norm();

    // Parametrization of the panel
    auto gamma = [&](double t) { return 0.5 * (b + a) + 0.5 * t * (b - a); };

    double integral1 = 0, integral2 = 0;
    // Quadrature loop
    for (int pt = 0; pt < order; ++pt) {
      double t = X[pt];
      integral1 += 0.5 * h * W[pt] * u2(gamma(t)) * (1 - t) * 0.5;
      integral2 += 0.5 * h * W[pt] * u2(gamma(t)) * (1 + t) * 0.5;
    }

    g_beta(aidx) += integral1;
    g_beta(bidx) += integral2;
  }

  // Projecting the dirichlet data
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(M11); // Analyze the sparsity pattern of A
  solver.factorize(M11);      // Factorize A
  Eigen::VectorXd g_projected = solver.solve(g_beta); // Solve the system Ax = b

  // Eigen::VectorXd g_projected = M11.lu().solve(g_beta);

  // Nodal interpolation for Dirichlet data, nodal values as coefficient
  Eigen::VectorXd g_interpolated = Eigen::VectorXd::Zero(mesh.numVertices());

  for (int i = 0; i < mesh.numVertices(); ++i) {
    std::cout << "vtx = " << mesh.getVertex(i).transpose() << std::endl;
    g_interpolated(i) = u2(mesh.getVertex(i));
  }

  Eigen::VectorXd rhs_interpolated = (0.5 * M01 + K) * g_interpolated;
  Eigen::VectorXd rhs_projected = (0.5 * M01 + K) * g_projected;

  Eigen::VectorXd sol_interpolated = V.lu().solve(rhs_interpolated);
  Eigen::VectorXd sol_projected = V.lu().solve(rhs_projected);

  Eigen::VectorXd sol_ex = 0 * sol_interpolated;

  // Exact values for the Neumann Data
  auto gradu2 = [](Eigen::Vector2d x) {
    Eigen::Vector2d out;
    out << std::sinh(x(0)) * std::cos(x(1)), -std::cosh(x(0)) * std::sin(x(1));
    return out;
  };

  for (int i = 0; i < mesh.numElements(); ++i) {
    Eigen::Vector2d a, b;
    std::tie(a, b) = mesh.getElementVertices(i);
    Eigen::Vector2d x = 0.5 * (a + b);
    Eigen::Vector2d tgt = b - a;
    Eigen::Vector2d nrm(tgt(1), -tgt(0));
    nrm /= nrm.norm();
    sol_ex(i) = gradu2(x).dot(nrm);
  }

  Eigen::Vector2d eval_pt(0.2, 0.13);

  std::cout << "exact Neumann data = \n" << sol_ex.transpose() << std::endl;
  std::cout << "Solved Neumann data (proj) = \n"
            << sol_projected.transpose() << std::endl;
  std::cout << "Solved Neumann data (interp) = \n"
            << sol_interpolated.transpose() << std::endl;

  // Point evaluation of the potential
  Eigen::VectorXd eval_proj, eval_interp;

  evaluateV(eval_proj, mesh, sol_projected, eval_pt.transpose(), 0);
  evaluateV(eval_interp, mesh, sol_interpolated, eval_pt.transpose(), 0);

  double exact_potential = u2(eval_pt);

  std::cout << "proj pot = " << eval_proj << std::endl;
  std::cout << "interp pot = " << eval_interp << std::endl;
  std::cout << "exact pot = " << exact_potential << std::endl;

  return 0;
}
