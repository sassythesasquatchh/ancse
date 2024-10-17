////////////////////////////////////////////////////////////////////////////////
/// \file evaluateN.cpp
/// \brief This file contains the function evaluateN that evaluates the Newton
///        potential on any number of evaluation points in \f$R^2\f$.
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

#include "evaluateN.hpp"

#include "constants.hpp"

namespace mysolution {
void evaluateN(Eigen::VectorXd& s,
               const Eigen::Matrix<double, Eigen::Dynamic, 2>& vertices,
               const Eigen::Matrix<int, Eigen::Dynamic, 3>& triangles,
               const Eigen::VectorXd& f, const Eigen::MatrixXd& x) {
  int nX = x.rows();
  // Initialize output vector
  s.resize(nX);
  s.setZero();

  Eigen::MatrixXd nodes(3, 2);
  // run over evaluation points
  for (int m = 0; m < nX; ++m) {
    // run over triangles
    for (int n = 0; n < triangles.rows(); ++n) {
      nodes << vertices.row(triangles(n, 0)), vertices.row(triangles(n, 1)),
          vertices.row(triangles(n, 2));

      s(m) += f[n] * newtonPotential(nodes, x.row(m));
    }

    s(m) /= (-2 * M_PI);
  }
}
}  // namespace mysolution
