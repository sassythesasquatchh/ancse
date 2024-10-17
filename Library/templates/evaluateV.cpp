///////////////////////////////////////////////////////////////////////////////
/// \file evaluateV.hpp
/// \brief This file contains the function evaluateV that evaluates the single
///        layer potential operator V tilde on any number of evaluation points
///        within the domain \f$Omega\f$.
///
///  This file contains only the implementation. For extensive documentation
///  consult the corresponding header-file.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////

#include "evaluateV.hpp"

#include <iostream>

#include "constants.hpp"

namespace mysolution {
void evaluateV(Eigen::VectorXd& Vphi_x, const BoundaryMesh& mesh,
               const Eigen::VectorXd& phi, const Eigen::MatrixXd& x,
               double eta) {
  int nX = x.rows();
  int nE = mesh.numElements();
  // Initialize output vector
  Vphi_x.resize(nX);
  Vphi_x.setZero();

  // Get quadrature points and weights
  const double* qp = getGaussPoints(GAUSS_ORDER);
  const double* qw = getGaussWeights(GAUSS_ORDER);

  /* For each boundary element Ej = [a,b], we calculate (a+b)/2
   * and (b-a)/2, because we need these vectors multiple times below.
   */
  Eigen::VectorXd lengthE(nE);
  Eigen::MatrixXd mEl(nE, 2);
  Eigen::MatrixXd dEl(nE, 2);
  // traverse the elements
  for (int j = 0; j < nE; ++j) {
    // get vertices indices and coordinates for Ej=[a,b]
    Eigen::Vector2d a, b;
    std::tie(a, b) = mesh.getElementVertices(j);
    // fill the vectors
    lengthE(j) = (b - a).norm();
    mEl.row(j) = 0.5 * (a + b);
    dEl.row(j) = 0.5 * (a - b);
  }

  for (int i = 0; i < nX; ++i) {
    // save current point for readibility
    const Eigen::Vector2d& xi = x.row(i);

    for (int j = 0; j < nE; ++j) {
      // get vertices indices and coordinates for Ej=[a,b]
      int aidx = mesh.getElementVertex(j, 0);
      int bidx = mesh.getElementVertex(j, 1);
      const Eigen::Vector2d& a = mesh.getVertex(aidx);
      const Eigen::Vector2d& b = mesh.getVertex(bidx);

      // check admissibility
      double dist_x_Ej = distancePointToSegment(xi, a, b);
      if (lengthE(j) <= eta * dist_x_Ej) {
        double sum = 0.;
        for (int k = 0; k < GAUSS_ORDER; ++k) {
          // transform point
          const Eigen::Vector2d& s = mEl.row(j) + qp[k] * dEl.row(j);
          sum += phi(j) * qw[k] * lengthE(j) * log((xi - s).squaredNorm());
          std::cout << " wii " << i << " , " << j << " , " << k << std::endl;
        }
        Vphi_x(i) -= sum;
      } else {
        Vphi_x(i) -= phi(j) * lengthE(j) *
                     slp(0, dEl.row(j), xi - mEl.row(j).transpose());
      }
    }
    Vphi_x[i] /= 8. * M_PI;
  }
}
}  // namespace mysolution
