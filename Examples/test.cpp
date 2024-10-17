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
#include "geometry.hpp"

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Vector type
typedef Eigen::VectorXd Vector;

int main(int, char **) {

  BoundaryMesh mesh("Lshape");

  Eigen::MatrixXd W;
  computeW(W, mesh, 0.001);

  Eigen::MatrixXd V;
  computeV(V, mesh, 0.001);

  Eigen::MatrixXd K;
  computeK(K, mesh, 0.001);

  Eigen::SparseMatrix<double> M(mesh.numVertices(), mesh.numVertices());
  computeM11(M, mesh);

  Eigen::SparseMatrix<double> M2(mesh.numElements(), mesh.numVertices());
  computeM01(M2, mesh);

  auto u1 = [](Eigen::Vector2d x) { return x(0); };

  auto u2 = [](Eigen::Vector2d x) { return std::cosh(x(0)) * std::cos(x(1)); };

  // Nodal interpolation for Dirichlet data, nodal values as coefficient
  Eigen::VectorXd dirichlet_coeff = Eigen::VectorXd::Zero(mesh.numVertices());

  for (int i = 0; i < mesh.numVertices(); ++i) {
    std::cout << "vtx = " << mesh.getVertex(i).transpose() << std::endl;
    dirichlet_coeff(i) = u2(mesh.getVertex(i));
  }

  std::cout << "Dirichlet coeff = \n" << dirichlet_coeff << std::endl;

  Eigen::VectorXd rhs = (0.5 * M2 + K) * dirichlet_coeff;

  Eigen::VectorXd sol = V.lu().solve(rhs);
  std::cout << "solved Neumann data = \n" << sol.transpose() << std::endl;

  Eigen::VectorXd sol_ex = 0 * sol;

  // Exact values for the Neumann Data
  auto gradu1 = [](Eigen::Vector2d x) {
    Eigen::Vector2d out(1, 0);
    return out;
  };

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

  return 0;
}
