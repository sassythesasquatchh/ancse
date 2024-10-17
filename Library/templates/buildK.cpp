///////////////////////////////////////////////////////////////////////////////
/// \file buildK.cpp
/// \brief This file provides functions to read the input-parameters, which are
///        relevant to compute the Galerkin-Matrix K corresponding to the
///        double-layer potential. The matrix is given by
///  \f[ K_{ij} = -\frac{1}{2 \pi} \int_{Ei}\int_{supp \phi_j} \frac{<y-x,n>}
///               {\vert y-x \vert^2} \phi_j(y) ds_y ds_x.  \f]
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

#include "buildK.hpp"

#include <cmath>

#include "constants.hpp"
#include "doubleLayerPotential.hpp"

namespace mysolution {
/* SAM_LISTING_BEGIN_1 */
void computeK(Eigen::MatrixXd &K, const BoundaryMesh &mesh, double eta) {
  int nE = mesh.numElements();
  int nC = mesh.numVertices();
  // Matrix returned through reference: resize and initialize matrix
  K.resize(nE, nC);
  K.setZero();
  double I0 = 0.0, I1 = 0.0;

  // Outer loop over the panels
  for (int j = 0; j < nE; ++j) {
    // TODO:
    // get vertices indices and coordinates for panel \cob{$\pan_j = [\Ba,\Bb]$}

    // Inner loop over the panels
    for (int i = 0; i < nE; ++i) {
      // TODO:
      // get vertices indices and coordinates for panel \cob{$\pan_i = [\Bc,\Bd]$}

      // TODO:
      // Handle parallel panels !

      // TODO:
      // compute entries of $1\times2$ local interaction matrix using the function
      // computeKij. Refer to doubleLayerPotential cpp/hpp or the lecture document
      // for more details

      //computeKij();

      // TODO:
      // Local to global mapping to populate the Galerkin matrix K

    }  // endfor
  }  // endfor
}
/* SAM_LISTING_END_1 */

void computeK00(Eigen::MatrixXd &K, const BoundaryMesh &mesh, double eta) {
  int nE = mesh.numElements();
  // Matrix returned through reference: resize and initialize matrix
  K.resize(nE, nE);
  K.setZero();
  double I0 = 0.0, I1 = 0.0;

  // \com{outer loop}: traverse the panels
  for (int j = 0; j < nE; ++j) {
    // get vertices indices and coordinates for panel \cob{$\pan_j = [\Ba,\Bb]$}
    Eigen::Vector2d a, b;
    std::tie(a, b) = mesh.getElementVertices(j);

    // \com{inner loop}: traverse the panels
    for (int i = 0; i < nE; ++i) {
      // get vertices indices and coordinates for panel \cob{$\pan_i = [\Bc,\Bd]$}
      Eigen::Vector2d c, d;
      std::tie(c, d) = mesh.getElementVertices(i);
      // Zero contribution for parallel panels !
      double lindep1 = fabs((a - c)[0] * (b - a)[1] - (a - c)[1] * (b - a)[0]);
      double lindep2 = fabs((a - d)[0] * (b - a)[1] - (a - d)[1] * (b - a)[0]);

      if (lindep1 > EPS * (a - c).norm() || lindep2 > EPS * (a - d).norm()) {
        // compute entries of $1\times2$ interaction matrix
        // double I0=0.0, I1=0.0;
        computeKij(&I0, &I1, eta, a, b, c, d);
        // distribute values to matrix entries
        K(j, i) += 2 * I0;
      }  // endif

    }  // endfor
  }  // endfor
}
}  // namespace mysolution
