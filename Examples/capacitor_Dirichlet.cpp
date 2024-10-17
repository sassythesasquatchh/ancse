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
#include "capacitorMesh.hpp"
#include "geometry.hpp"
#include "potentials.hpp"

int main(int, char **) {
  // Create the capacitor mesh
  BoundaryMesh mesh = createCapacitorMesh(1);

  // Assembling the Galerkin matrices
  // Single layer Galerkin matrix
  Eigen::MatrixXd V;
  computeV(V, mesh, 0);
  // Double layer Galerkin matrix
  Eigen::MatrixXd K;
  computeK(K, mesh, 0);

  // Assembling mass matrix with test space P0 and trial space P1
  // Rows <-> test space, columns <-> trial space
  Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
  computeM01(M01, mesh);

  // Nodal interpolation for Dirichlet data, nodal values as coefficient
  Eigen::VectorXd dirichlet_coeff = Eigen::VectorXd::Zero(mesh.numVertices());

  for (int i = 0; i < mesh.numVertices(); ++i) {
    // Computing the nodal values
    dirichlet_coeff(i) = u1(mesh.getVertex(i));
  }

  // Assembling RHS and solving the system
  Eigen::VectorXd rhs = (0.5 * M01 + K) * dirichlet_coeff;
  // Solving the linear system
  Eigen::VectorXd sol = V.lu().solve(rhs);
  std::cout << "solved Neumann data = \n" << sol.transpose() << std::endl;

  // Calculating Neumann trace from the known gradu1
  Eigen::VectorXd sol_ex = NeumannTraceInterpolation(gradu1, mesh);
  std::cout << "exact Neumann data = \n" << sol_ex.transpose() << std::endl;

  return 0;
}
