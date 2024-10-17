#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "BoundaryMesh.hpp"
#include "L2ProjectionToP1.hpp"
#include "L2WithP1Basis.hpp"
#include "NeumannTraceInterpolation.hpp"
#include "buildHypsingStabilization.hpp"
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
  // Capacitor mesh
  BoundaryMesh mesh = createCapacitorMesh(1);

  // Assembling Galerkin matrices
  // Hypersingular operator (unstabilized, singular)
  Eigen::MatrixXd Wmat;
  computeW(Wmat, mesh, 0);

  // Stabilization matrix for Hypersingular operator
  Eigen::MatrixXd WmatStab;
  buildHypsingStabilization(WmatStab, mesh);
  // Adding the stabilization to the Hypersingular operator
  WmatStab += Wmat;

  // Double layer operator
  Eigen::MatrixXd K;
  computeK(K, mesh, 0);

  // Assembling mass matrix
  // Test space (rows) = P0, Trial space (columns) = P1
  Eigen::SparseMatrix<double> M01(mesh.numElements(), mesh.numVertices());
  computeM01(M01, mesh);

  // Projecting the dirichlet data
  Eigen::VectorXd g_projected = L2ProjectionToP1(u1, mesh);

  // Getting the Neumann Trace in the discrete space via interpolation
  Eigen::VectorXd neumann_trace = NeumannTraceInterpolation(gradu1, mesh);

  // RHS of the linear system
  Eigen::VectorXd rhs =
      0.5 * M01.transpose() * neumann_trace - K.transpose() * neumann_trace;

  // Solving the linear system
  Eigen::VectorXd dirichlet_sol = WmatStab.lu().solve(rhs);

  // Solving with conjugate gradient
  Eigen::ConjugateGradient<Eigen::MatrixXd> cg;
  cg.compute(Wmat);
  Eigen::VectorXd dirichlet_sol_cg = cg.solve(rhs);

  std::cout << "W g_projected = " << (Wmat * g_projected).transpose()
            << std::endl;
  std::cout << "W dirichlet_sol = " << (Wmat * dirichlet_sol).transpose()
            << std::endl;
  std::cout << "W dirichlet_sol_cg = " << (Wmat * dirichlet_sol_cg).transpose()
            << std::endl;

  return 0;
}
