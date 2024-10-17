#ifndef NEUMANNTRACEINTERPOLATIONmysolution
#define NEUMANNTRACEINTERPOLATIONmysolution

#include <Eigen/Dense>

#include "BoundaryMesh.hpp"
#include "geometry.hpp"

namespace mysolution {
template <typename GRADU>
Eigen::VectorXd NeumannTraceInterpolation(GRADU gradu, BoundaryMesh &mesh) {
  Eigen::VectorXd neumann_trace(mesh.numElements());

  for (int i = 0; i < mesh.numElements(); ++i) {
    Eigen::Vector2d a, b;
    std::tie(a, b) = mesh.getElementVertices(i);
    Eigen::Vector2d x = 0.5 * (a + b);
    Eigen::Vector2d nrm = unitNormal(a, b);
    neumann_trace(i) = gradu(x).dot(nrm);
  }
  return neumann_trace;
}
}  // namespace mysolution

#endif