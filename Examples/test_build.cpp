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

  Eigen::Vector2d a, b, c;
  a << 1.0, 4.5;
  b << 2., 2;
  c << 5.1, 3.;
  std::cout << distancePointToSegment(a, b, c) << std::endl;

  BoundaryMesh mesh("Lshape");

  Eigen::MatrixXd W;
  computeW(W, mesh, 0.001);

  const std::string fnameW = "W_Lshape_Hilbert.dat";
  std::ofstream outW(fnameW.c_str());
  outW << std::setprecision(18) << W;
  outW.close();

  Eigen::MatrixXd V;
  computeV(V, mesh, 0.001);

  const std::string fnameV = "V_Lshape_Hilbert.dat";
  std::ofstream outV(fnameV.c_str());
  outV << std::setprecision(18) << V;
  outV.close();

  Eigen::MatrixXd K;
  computeK(K, mesh, 0.001);

  const std::string fnameK = "K_Lshape_Hilbert.dat";
  std::ofstream outK(fnameK.c_str());
  outK << std::setprecision(18) << K;
  outK.close();

  Eigen::SparseMatrix<double> M(mesh.numVertices(), mesh.numVertices());
  computeM11(M, mesh);
  const std::string fnameM = "M_Lshape_Hilbert.dat";
  std::ofstream outM(fnameM.c_str());
  outM << std::setprecision(18) << M;
  outM.close();

  Eigen::SparseMatrix<double> M2(mesh.numElements(), mesh.numVertices());
  computeM01(M2, mesh);
  const std::string fnameM2 = "M2_Lshape_Hilbert.dat";
  std::ofstream outM2(fnameM2.c_str());
  outM2 << std::setprecision(18) << M2;
  outM2.close();

  std::cout << " Done " << std::endl;

  // loaded them in matlab and compared with Hilbert original implementation.
  // the norm of the difference is zero :)
}
