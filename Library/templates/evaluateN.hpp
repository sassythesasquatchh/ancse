///////////////////////////////////////////////////////////////////////////////
/// \file evaluateN.hpp
/// \brief This file contains the function evaluateN that evaluates the Newton
///        potential on any number of evaluation points in \f$R^2\f$.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _EVALUATEN_HPPmysolution
#define _EVALUATEN_HPPmysolution

#include "newtonPotential.hpp"

namespace mysolution {
/**
 *  evaluateN evaluates the Newton potential N on a number of evaluation points.
 *
 *  @param[out] s vector
 *  @param[in] vertices  (nV x 2) matrix containing the coordinates of the
 *                       vertices of the domain mesh.
 *  @param[in] triangles  (nT x 3) matrix containing the indices of the vertices
 *                        conforming the triangles of the domain mesh.
 *  @param[in] f  nT vector such that phi_i=phi(z_i).
 *  @param[in] x  (nX x 2) matrix that contains the evaluation points.
 */
void evaluateN(Eigen::VectorXd& s,
               const Eigen::Matrix<double, Eigen::Dynamic, 2>& vertices,
               const Eigen::Matrix<int, Eigen::Dynamic, 3>& triangles,
               const Eigen::VectorXd& f, const Eigen::MatrixXd& x);
}  // namespace mysolution

#endif
