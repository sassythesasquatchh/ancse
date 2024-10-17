#ifndef LSHAPE_MESHmysolution
#define LSHAPE_MESHmysolution

#include <Eigen/Dense>

#include "BoundaryMesh.hpp"

namespace mysolution {
BoundaryMesh createLshapeMesh(int k) {
  Eigen::MatrixXd rotmat(2, 2);
  rotmat << std::cos(M_PI / 4), -std::sin(M_PI / 4), std::sin(M_PI / 4),
      std::cos(M_PI / 4);
  int eltsPerBigSide = std::pow(2, k + 1);
  int eltsPerSmallSide = std::pow(2, k);

  int totalElts = 4 * eltsPerBigSide;
  Eigen::MatrixXd vertices = Eigen::MatrixXd::Zero(totalElts, 2);
  Eigen::MatrixXi elements = Eigen::MatrixXi::Zero(totalElts, 2);

  elements.col(0) = Eigen::VectorXi::LinSpaced(totalElts, 0, totalElts - 1);
  elements.col(1) = Eigen::VectorXi::LinSpaced(totalElts, 1, totalElts);
  elements(totalElts - 1, 1) = 0;

  Eigen::Vector2d SW(-0.25, -0.25);
  Eigen::Vector2d SE(0.25, -0.25);
  Eigen::Vector2d NE(0.25, 0.25);
  Eigen::Vector2d NW(-0.25, 0.25);
  Eigen::Vector2d T(0, 0.25);
  Eigen::Vector2d C(0, 0);
  Eigen::Vector2d L(-0.25, 0);

  for (int i = 0; i < eltsPerBigSide; ++i) {
    // Bottom edge
    vertices.row(i) = SW + ((double)i) / (double)eltsPerBigSide * (SE - SW);
    // Right edge
    vertices.row(eltsPerBigSide + i) =
        SE + ((double)i) / (double)eltsPerBigSide * (NE - SE);
    // Top edge
    if (i < eltsPerSmallSide) {
      vertices.row(2 * eltsPerBigSide + i) =
          NE + ((double)i) / (double)eltsPerSmallSide * (T - NE);
    } else {
      vertices.row(2 * eltsPerBigSide + i) =
          T +
          ((double)(i - eltsPerSmallSide)) / (double)eltsPerSmallSide * (C - T);
    }

    if (i < eltsPerSmallSide) {
      vertices.row(3 * eltsPerBigSide + i) =
          C + ((double)i) / (double)eltsPerSmallSide * (L - C);
    } else {
      vertices.row(3 * eltsPerBigSide + i) =
          L + ((double)(i - eltsPerSmallSide)) / (double)eltsPerSmallSide *
                  (SW - L);
    }
  }
  vertices = (rotmat * vertices.transpose()).transpose();
  BoundaryMesh out(vertices, elements);
  return out;
}
}  // namespace mysolution

#endif