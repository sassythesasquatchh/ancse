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
    int first_vertex_idx = mesh.getElementVertex(i, 0);
    int second_vertex_idx = mesh.getElementVertex(i, 1);
    const Eigen::Vector2d& first_vertex = mesh.getVertex(first_vertex_idx);
    const Eigen::Vector2d& second_vertex = mesh.getVertex(second_vertex_idx);

    double h = (second_vertex - first_vertex).norm();

    // Parametrization of the panel
    auto gamma = [&](double t) { return first_vertex+0.5*(t+1)*(second_vertex-first_vertex); };

    double integral1 = 0, integral2 = 0;
    // Quadrature loop
    for (int pt = 0; pt < order; ++pt) {
      // t is the quadrature point on the reference element, W[pt] is the corresponding quadrature weight
      double t = X[pt];

      // The following equations arise from the RHS integral.
      // To get the integral in this form, one must bring the integral
      // to the reference element through a reparametrization. 
      // This transformation introduces the factor 0.5 * h. 
      // By subsituting x = gamma(t) (where gamma(t) is the parametrisation of the 
      // element), we get f(gamma(t)) for the RHS function and 0.5*(1-t), 0.5*(1+t)
      // for the first and second local shape functions, respectively, on the reference element.

      // Integral of the reference local shape function that slopes down from 
      // one at the first vertex to zero at the second vertex
      integral1 += 0.5 * h * W[pt] * f(gamma(t)) * (1 - t) * 0.5;

      // Integral of the reference local shape function that slopes up from
      // zero at the first vertex to one at the second vertex
      integral2 += 0.5 * h * W[pt] * f(gamma(t)) * (1 + t) * 0.5;
    }

    // Add the contribution of the local shape functions to their corresponding
    // global vertex.
    out(first_vertex_idx) += integral1;
    out(second_vertex_idx) += integral2;
  }

  return out;
}
}  // namespace mysolution

#endif