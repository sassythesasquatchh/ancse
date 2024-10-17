///////////////////////////////////////////////////////////////////////////////
/// \file buildV.cpp
/// \brief This file provides functions for the computation of the single-layer
///  potential matrix V defined by
///  \f[ V_{jk} = < V \chi_j , \chi_k >  \f]
///  with characteristic functions \f$\chi_j\f$ corresponding to a line segment
///  Tj.
///  The single-layer potential is defined by
///  \f[ V \chi_j = -\frac{1}{2\pi} \int_{Ej} \log{ \vert x-y \vert } ds_y. \f]
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
#include <cmath>
#include "buildV.hpp"
#include "constants.hpp"
#include "singleLayerPotential.hpp"


//------------------------------------------------------------------------------
void computeV(Eigen::MatrixXd& V, const BoundaryMesh& mesh, double eta)
{
  assert(eta >= 0);

  // resize matrix
  int nE = mesh.numElements();
  V.resize(nE,nE);
  // outer loop traversing all panels 
  for (int i=0; i<nE; ++i)
  {
    // get endpoint indices and coordinates for $i$-th panel
    int aidx = mesh.getElementVertex(i,0);
    int bidx = mesh.getElementVertex(i,1);
    const Eigen::Vector2d& a = mesh.getVertex(aidx);
    const Eigen::Vector2d& b = mesh.getVertex(bidx);

    // inner loop through all panels
    for (int j=i; j<nE; ++j)
    {
      // get vertices indices and coordinates for Ej=[c,d]
      int cidx = mesh.getElementVertex(j,0);
      int didx = mesh.getElementVertex(j,1);
      const Eigen::Vector2d& c = mesh.getVertex(cidx);
      const Eigen::Vector2d& d = mesh.getVertex(didx);

      // compute contribution of a pair of panels
      double tmp = computeVij(a, b, c, d, eta);
      // distribute it among the matrix entries
      V(i,j) = tmp; V(j,i) = V(i,j);
    }
  }
  
}

