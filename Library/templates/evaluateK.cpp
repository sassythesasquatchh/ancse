///////////////////////////////////////////////////////////////////////////////
/// \file evaluateK.cpp
/// \brief This file contains the function evaluateK that evaluates the double
///        layer potential operator K tilde on any number of evaluation points
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

#include "evaluateK.hpp"

#include "constants.hpp"

namespace mysolution {
void evaluateK(Eigen::VectorXd& Kfx, const BoundaryMesh& mesh,
               const Eigen::VectorXd& fh, const Eigen::MatrixXd& x,
               double eta) {
  int nX = x.rows();
  int nE = mesh.numElements();
  // Initialize output vector
  Kfx.resize(nX);
  Kfx.setZero();

  // Get quadrature points and weights
  const double* qp = getGaussPoints(GAUSS_ORDER);
  const double* qw = getGaussWeights(GAUSS_ORDER);

  // Auxiliary vector containing elements info.
  Eigen::VectorXd lengthE(nE);
  Eigen::MatrixXd mEl(nE, 2);
  Eigen::MatrixXd dEl(nE, 2);
  Eigen::MatrixXd nEl(nE, 2);

  // Traverse the elements
  for (int j = 0; j < nE; ++j) {
    // Get vertices indices and coordinates for Ej=[a,b]
    Eigen::Vector2d a, b;
    std::tie(a, b) = mesh.getElementVertices(j);
    // Fill the vectors
    lengthE(j) = (b - a).norm();
    mEl.row(j) = 0.5 * (a + b);
    dEl.row(j) = 0.5 * (b - a);
    nEl.row(j) = unitNormal(a, b);
  }

  // Traverse the evaluation points
  for (int i = 0; i < nX; ++i) {
    // Save current point for readibility
    const Eigen::Vector2d& xi = x.row(i);

    // Traverse elements
    for (int j = 0; j < nE; ++j) {
      // Get vertices indices and coordinates for Ej=[a,b]
      int aidx = mesh.getElementVertex(j, 0);
      int bidx = mesh.getElementVertex(j, 1);
      const Eigen::Vector2d& a = mesh.getVertex(aidx);
      const Eigen::Vector2d& b = mesh.getVertex(bidx);

      // Check admissibility
      double dist_xi_Ej = distancePointToSegment(xi, a, b);
      if (lengthE(j) <= dist_xi_Ej * eta) {
        // Integrate
        double sum = 0.;
        for (int k = 0; k < GAUSS_ORDER; ++k) {
          // Transform point
          const Eigen::Vector2d& s = mEl.row(j) + qp[k] * dEl.row(j);
          // Add contribution
          sum += qw[k] * (s - xi).dot(nEl.row(j)) *
                 (fh[aidx] + ((1 + qp[k]) * (fh(bidx) - fh(aidx)) / 2)) /
                 (s - xi).squaredNorm();
        }

        Kfx(i) -= lengthE(j) * sum;
      } else {
        // Compute integral analitically
        const Eigen::Vector2d& aux = mEl.row(j) - xi.transpose();
        if (fabs(aux.dot(nEl.row(j))) > EPS) {
          double commonFactor = aux.dot(nEl.row(j)) * lengthE(j);
          double prod1 = commonFactor * dlp(0, dEl.row(j), aux) / 2.;
          double prod2 = commonFactor * dlp(1, dEl.row(j), aux) / 2.;
          Kfx(i) -=
              (((fh(aidx) * prod1 + fh(bidx) * prod1) + fh(bidx) * prod2) -
               fh(aidx) * prod2);
        }
      }

    }  // end elements' for loop

    Kfx(i) /= (4. * M_PI);
  }  // end evaluation points' for loop
}
}  // namespace mysolution
