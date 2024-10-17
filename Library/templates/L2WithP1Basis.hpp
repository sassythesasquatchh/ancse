#ifndef L2WITHP1BASISmysolution
#define L2WITHP1BASISmysolution

#include <Eigen/Dense>

#include "BoundaryMesh.hpp"
#include "gaussQuadrature.hpp"

namespace mysolution {
template <typename FUNC>
Eigen::VectorXd L2WithP1Basis(FUNC f, BoundaryMesh& mesh) {
  Eigen::VectorXd out = Eigen::VectorXd::Zero(mesh.numVertices());
  unsigned order = 32;
  // Get the Gauss quadrature points and weights for the interval [-1,1]
  const double* X = getGaussPoints(order);
  const double* W = getGaussWeights(order);

  // Local to global assembly. Start traversing the panels
  for (int i = 0; i < mesh.numElements(); i++) {
    // identify element's vertices
    int aidx = mesh.getElementVertex(i, 0);
    int bidx = mesh.getElementVertex(i, 1);
    const Eigen::Vector2d& a = mesh.getVertex(aidx);
    const Eigen::Vector2d& b = mesh.getVertex(bidx);

    double h = (b - a).norm();

    // Parametrization of the panel
    auto gamma = [&](double t) { return 0.5 * (b + a) + 0.5 * t * (b - a); };

    double integral1 = 0, integral2 = 0;
    // Quadrature loop
    for (int pt = 0; pt < order; ++pt) {
      double t = X[pt];
      integral1 += 0.5 * h * W[pt] * f(gamma(t)) * (1 - t) * 0.5;
      integral2 += 0.5 * h * W[pt] * f(gamma(t)) * (1 + t) * 0.5;
    }

    out(aidx) += integral1;
    out(bidx) += integral2;
  }

  return out;
}
}  // namespace mysolution

#endif