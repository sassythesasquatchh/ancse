#ifndef L2PROJECTIONTOP1
#define L2PROJECTIONTOP1

#include "BoundaryMesh.hpp"
#include "L2WithP1Basis.hpp"
#include "buildM.hpp"
#include "gaussQuadrature.h"

template <typename FUNC>
Eigen::VectorXd L2ProjectionToP1(FUNC g, BoundaryMesh& mesh) {
  // L2 projection of Dirichlet data g
  // Idea: projection g_h can be written in terms of P1 basis functions and satisfies
  // the orthogonality relation: (g - g_h, b_i) = 0 for all i.
  // Collecting coefficients of g_h^j in terms of the basis b_j gives the linear system
  // \sum_{j=1..N} g_h^j (b_j,b_i) = (g, b_i)
  // (b_j, b_i) is the matrix M11

  // Computing (g, b_i) for all basis functions b_i in the space P1
  // Can be done by simply using Gauss quadrature
  Eigen::VectorXd g_b = L2WithP1Basis(g, mesh);

  // Assembling mass matrix
  // Test space (rows) = P1, Trial space (columns) = P1
  Eigen::SparseMatrix<double> M11(mesh.numVertices(), mesh.numVertices());
  computeM11(M11, mesh);

  // Projecting the dirichlet data by solving the linear system
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(M11);
  solver.factorize(M11);
  // Computing the coefficients g_h^j by solving the linear system
  Eigen::VectorXd g_projected = solver.solve(g_b);

  return g_projected;
}

#endif