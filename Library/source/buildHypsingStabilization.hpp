#include <Eigen/Dense>

#include "buildM.hpp"

#ifndef BUILDHYPSINGSTABILIZATIONHPP
#define BUILDHYPSINGSTABILIZATIONHPP

/**
 *  Stabilization matrix for the hypersingular IE using S1-elements.
 *
 *  The hypersingular integral equation reads \f$ W u = (1/2-K') \phi \f$
 *  where \f$\phi\f$ is the known Neumann data for the solution u of
 *  \f[ -\Delta u = f  \qquad in \Omega \\
 *      du/dn = \phi \qquad on \Gamma = \partial \Omega \f]
 *
 *  The corresponding stiffness matrix W of the hypersingular integral
 *  operator is usually built for S1 boundary elements. It therefore has a
 *  non-trivial kernel. This function computes the stabilization matrix
 *  \f[ S_{jk} = (\int_{\Gamma} \phi_j ds) (\int_{\Gamma} \phi_k ds) \f]
 *  which may be used to obtain a discrete system \f$ (W+S)u = RHS \f$,
 *  which now allows for a uniquely determined numerical solution of the
 *  hypersingular integral equation.
 *
 *  @param[out] S  (nC x nC) Stabilization matrix.
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 */
void buildHypsingStabilization(Eigen::MatrixXd &S, const BoundaryMesh &mesh) {
  int nC = mesh.numVertices();
  Eigen::SparseMatrix<double> M(nC, nC);
  computeM11(M, mesh);
  Eigen::VectorXd aux(nC);
  aux.setOnes();
  Eigen::VectorXd c = M * aux;

  S.resize(nC, nC);
  S = c * c.transpose();
}

#endif
