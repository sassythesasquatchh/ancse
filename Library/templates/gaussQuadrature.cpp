///////////////////////////////////////////////////////////////////////////////
/// \file gaussQuadrature.cpp
/// \brief This file provides functions to get the gauss points and gauss
///  weights for a certain order.
///
///  This file contains only the implementation. For extensive documentation
///  consult the corresponding header-file.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
/// VERSION: 3.1
/// (C) 2009-2013 HILBERT-Team '09, '10, '12
/// support + bug report:  hilbert@asc.tuwien.ac.at
///////////////////////////////////////////////////////////////////////////////
#include "gaussQuadrature.hpp"

#include <assert.h>
#include <stdlib.h>

#include "gaussTable.hpp"

namespace mysolution {
const double* getGaussPointsT(int points, int coordinate) {
  assert((coordinate == 0 || coordinate == 1) && (points == 7));

  if (coordinate == 0) {
    switch (points) {
      case 7:
        return gaussT7_x;
      default:
        return NULL;
    }
  } else {
    switch (points) {
      case 7:
        return gaussT7_y;
      default:
        return NULL;
    }
  }
}

const double* getGaussWeightsT(int points) {
  assert(points == 7);

  switch (points) {
    case 7:
      return gaussT7_w;
    default:
      return NULL;
  }
}

const double* getGaussPoints(int order) {
  assert(order == 2 || order == 4 || order == 8 || order == 16 || order == 32);

  switch (order) {
    case 2:
      return gauss2_points;
      /*    case 3:
      return gauss3_points; */
    case 4:
      return gauss4_points;
      /*    case 5:
      return gauss5_points;
    case 6:
      return gauss6_points;
    case 7:
      return gauss7_points; */
    case 8:
      return gauss8_points;
      /*    case 9:
      return gauss9_points;
    case 10:
      return gauss10_points;
    case 11:
      return gauss11_points;
    case 12:
      return gauss12_points;
    case 13:
      return gauss13_points;
    case 14:
      return gauss14_points;
    case 15:
      return gauss15_points; */
    case 16:
      return gauss16_points;
      /*    case 17:
      return gauss17_points;
    case 18:
      return gauss18_points;
    case 19:
      return gauss19_points;
    case 20:
      return gauss20_points;
    case 21:
      return gauss21_points;
    case 22:
      return gauss22_points;
    case 23:
      return gauss23_points;
    case 24:
      return gauss24_points;
    case 25:
      return gauss25_points;
    case 26:
      return gauss26_points;
    case 27:
      return gauss27_points;
    case 28:
      return gauss28_points;
    case 29:
      return gauss29_points;
    case 30:
      return gauss30_points;
    case 31:
      return gauss31_points; */
    case 32:
      return gauss32_points;
    default:
      return NULL;
  }
}

const double* getGaussWeights(int order) {
  assert(order == 2 || order == 4 || order == 8 || order == 16 || order == 32);

  switch (order) {
    case 2:
      return gauss2_weights;
      /*    case 3:
      return gauss3_weights; */
    case 4:
      return gauss4_weights;
      /*    case 5:
      return gauss5_weights;
    case 6:
      return gauss6_weights;
    case 7:
      return gauss7_weights; */
    case 8:
      return gauss8_weights;
      /*    case 9:
      return gauss9_weights;
    case 10:
      return gauss10_weights;
    case 11:
      return gauss11_weights;
    case 12:
      return gauss12_weights;
    case 13:
      return gauss13_weights;
    case 14:
      return gauss14_weights;
    case 15:
      return gauss15_weights; */
    case 16:
      return gauss16_weights;
      /*    case 17:
      return gauss17_weights;
    case 18:
      return gauss18_weights;
    case 19:
      return gauss19_weights;
    case 20:
      return gauss20_weights;
    case 21:
      return gauss21_weights;
    case 22:
      return gauss22_weights;
    case 23:
      return gauss23_weights;
    case 24:
      return gauss24_weights;
    case 25:
      return gauss25_weights;
    case 26:
      return gauss26_weights;
    case 27:
      return gauss27_weights;
    case 28:
      return gauss28_weights;
    case 29:
      return gauss29_weights;
    case 30:
      return gauss30_weights;
    case 31:
      return gauss31_weights; */
    case 32:
      return gauss32_weights;
    default:
      return NULL;
  }
}
}  // namespace mysolution
