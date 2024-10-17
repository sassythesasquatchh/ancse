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
    int j_first_idx = mesh.getElementVertex(j, 0);
    int j_second_idx = mesh.getElementVertex(j, 1);
    const Eigen::Vector2d &j_first_vertex = mesh.getVertex(j_first_idx);
    const Eigen::Vector2d &j_second_vertex = mesh.getVertex(j_second_idx);
    // Inner loop over the panels
    for (int i = 0; i < nE; ++i) {
      // TODO:
      // get vertices indices and coordinates for panel \cob{$\pan_i = [\Bc,\Bd]$}
      int i_first_idx = mesh.getElementVertex(i, 0);
      int i_second_idx = mesh.getElementVertex(i, 1);
      const Eigen::Vector2d &i_first_vertex = mesh.getVertex(i_first_idx);
      const Eigen::Vector2d &i_second_vertex = mesh.getVertex(i_second_idx);
      
      // TODO:
      // Handle parallel panels !

      // Note that the double layer boundary integral operator is zero over parallel panels.
      // The is because the vector from a point on one panel to a point on the other panel is
      // always perpendicular to the normal vector of either panel (importantly here, the second
      // panel). Therefore the dot product in the numerator is 0 for all pairs of points in the set
      // panel X panel'

      // If the panels are parallel, then 
      // (j_first_vertex - i_first_vertex) = c*(j_second_vertex - j_first_vertex)
      // and
      // (j_first_vertex - i_second_vertex) = c*(j_second_vertex - j_first_vertex)
      // where c is some constant. Forming equations in the x and y coordinates and substituting
      // leads to the following equations, that are zero iff the panels are parallel.
      double lindep1 = fabs((j_first_vertex - i_first_vertex)[0] * (j_second_vertex - j_first_vertex)[1] - (j_first_vertex - i_first_vertex)[1] * (j_second_vertex - j_first_vertex)[0]);
      double lindep2 = fabs((j_first_vertex - i_second_vertex)[0] * (j_second_vertex - j_first_vertex)[1] - (j_first_vertex - i_second_vertex)[1] * (j_second_vertex - j_first_vertex)[0]);

      // TODO:
      // compute entries of $1\times2$ local interaction matrix using the function
      // computeKij. Refer to doubleLayerPotential cpp/hpp or the lecture document
      // for more details
      if (lindep1 > EPS * (j_first_vertex - i_first_vertex).norm() || 
          lindep2 > EPS * (j_first_vertex - i_second_vertex).norm()) 
      {

      // Compute Kij evaluates the double layer operator over the given elements
      computeKij(&I0, &I1, eta, j_first_vertex, j_second_vertex, i_first_vertex, i_second_vertex);

      // TODO:
      // Local to global mapping to populate the Galerkin matrix K

      // The interaction matrix is 2x1, because the the Dirichlet data is approximated
      // by the S^0_1 space with a tent function basis (support over two elements) and
      // the Neumann data is approximated by the S^-1_0 space with a characteristic
      // function basis (support over one element).

      // A linear combination is used, because I0 and I1 are calculated with unusual 
      // reference shape functions. The reference shape function for I0 is a constant
      // function of 1/2 over the reference element. The reference shape function 
      // for I1 is a zero mean function with gradient 1/2 over the reference element.
      // The typical piecewise linear reference shape functions can be reconstructed 
      // from these unusual reference shape functions via the linear combinations below.

      // Entry corresponding to the basis function of the j-th element and the 
      // local shape funciton of the first vertex of the i-th element
      // NB. This is the local shape function that slopes down from 1 at the first vertex
      // of i to 0 at the second vertex of i.
      K(j, i_first_idx) += I0 - I1;

      // Entry corresponding to the basis function of the j-th element and the 
      // local shape function of the second vertex of the i-th element
      // NB. This is the local shape function that slopes up from 0 at the first vertex
      // of i to 1 at the second vertex of i.
      K(j, i_second_idx) += I0 + I1;
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
