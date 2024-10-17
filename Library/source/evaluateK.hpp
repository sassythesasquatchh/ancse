///////////////////////////////////////////////////////////////////////////////
/// \file evaluateK.hpp
/// \brief This file contains the function evaluateK that evaluates the double
///        layer potential operator K tilde on any number of evaluation points
///        within the domain \f$Omega\f$.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _EVALUATEK_HPP
#define _EVALUATEK_HPP

#include "geometry.hpp"
#include "doubleLayerPotential.hpp"
extern "C" {
#include "gaussQuadrature.h"
}
#include "BoundaryMesh.hpp"


/**
 *  evaluateK evaluates the double layer potential operator K tilde on a number
 *  of evaluation points.
 *
 *  @param[out] Kfx vector
 *  @param[in] mesh 2D BoundaryMesh (initialized with vertices and elements).
 *  @param[in] fh nC vector such that fh_i=f(z_i).
 *  @param[in] x  (nX x 2) matrix that contains the evaluation points.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void evaluateK(Eigen::VectorXd& Kfx, const BoundaryMesh& mesh, const Eigen::VectorXd& fh,
               const Eigen::MatrixXd& x, double eta);

#endif
