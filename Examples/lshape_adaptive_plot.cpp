#include "BoundaryMesh.hpp"
#include "LshapeMeshAdaptive.hpp"
#include "matplotlibcpp.h"
#include <iostream>

namespace plt = matplotlibcpp;

int main(int, char **) {
  // Create the capacitor mesh
  const BoundaryMesh mesh = createLshapeMeshAdaptive(3, 2);

  Eigen::MatrixXd coords = mesh.getMeshVertices();
  Eigen::MatrixXd coords_extended(coords.rows() + 1, coords.cols());
  coords_extended.block(0, 0, coords.rows(), coords.cols()) = coords;
  coords_extended.row(coords.rows()) = coords.row(0);

  plt::plot(coords_extended.col(0), coords_extended.col(1));
  plt::scatter(coords.col(0), coords.col(1), 12);
  plt::axis("equal");
  plt::show();

  return 0;
}
