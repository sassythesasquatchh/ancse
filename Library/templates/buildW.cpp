///////////////////////////////////////////////////////////////////////////////
/// \file buildW.cpp
/// \brief This file provides functions for the calculation of the hypersingular
///  integral operator matrix W that is defined by
///  \f[ W_{jk} = < W \phi_j, \phi_k >  \f]
///  with the hat functions \f$\phi_i\f$ corresponding to a piece-wise affine
///  function that is linear on any element of the boundary mesh and that is 1
///  on the i-th node z_i and 0 on any other node.
///
///  The matrix W is in fact calculated using the following relation between W
///  and the simple-layer potential V:
///  \f[ <Wu, v> = <Vu', v'>   \f]
///  In this formula, u' and v' denote the arc-length derivative of u and v,
///  respectively.
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
#include "buildW.hpp"

#include <cmath>

#include "constants.hpp"
#include "singleLayerPotential.hpp"

namespace mysolution {
//------------------------------------------------------------------------------
void computeW(Eigen::MatrixXd& W, const BoundaryMesh& mesh, double eta) {
  // resize and initialize matrix
  int nE = mesh.numElements();
  int nC = mesh.numVertices();
  W.resize(nC, nC);
  W.setZero();

  // traverse the elements
  for (int i = 0; i < nE; ++i) {
    // get vertices indices and coordinates for Ei=[a,b]
    int aidx = mesh.getElementVertex(i, 0);
    int bidx = mesh.getElementVertex(i, 1);
    const Eigen::Vector2d& a = mesh.getVertex(aidx);
    const Eigen::Vector2d& b = mesh.getVertex(bidx);
    // traverse the elements
    for (int j = 0; j < nE; ++j) {
      // get vertices indices and coordinates for Ej=[c,d]
      int cidx = mesh.getElementVertex(j, 0);
      int didx = mesh.getElementVertex(j, 1);
      const Eigen::Vector2d& c = mesh.getVertex(cidx);
      const Eigen::Vector2d& d = mesh.getVertex(didx);
      // compute elements' contribution
      double tmp = computeWij(a, b, c, d, eta);
      // distribute among the matrix entries
      W(aidx, cidx) += tmp;
      W(aidx, didx) -= tmp;
      W(bidx, cidx) -= tmp;
      W(bidx, didx) += tmp;
    }
  }
}
}  // namespace mysolution
