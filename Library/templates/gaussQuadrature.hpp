///////////////////////////////////////////////////////////////////////////////
/// \file gaussQuadrature.hpp
/// \brief This file provides functions to get the gauss points and gauss
///  weights for a certain order.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
/// VERSION: 3.1
/// (C) 2009-2013 HILBERT-Team '09, '10, '12
/// support + bug report:  hilbert@asc.tuwien.ac.at
///////////////////////////////////////////////////////////////////////////////
#ifndef GAUSS_QUADRATURE_H_GUARDmysolution
#define GAUSS_QUADRATURE_H_GUARDmysolution

namespace mysolution {
/**
 *  The function returns quadrature points for x-direction (coordinate=0) or
 *  y-direction (coordinate=1).
 *
 *  @param[in] points  Only points=7 is permitted at the moment
 *  @param[in] coordinate  An integer of the set {0,1}
 *  @return A constant array of doubles containing the gauss points to the
 *          given array size (=points).
 */
const double* getGaussPointsT(int points, int coordinate);

/**
 *  The function returns quadrature weights for gauss quadrature On the
 *  reference triangle conv{(0,0),(1,0),(0,1)} according to the size of the
 *  rule (=points).
 *  @param[in] points  Inly points=7 is permitted at the moment
 *  @return A constant array of doubles containing the gauss weights to the
 *          given array size (=points).
 */
const double* getGaussWeightsT(int points);

/**
 *  The function returns quadrature points on the interval [-1,1] for the
 *  Gaussian quadrature rule of a given order. It does not create a copy of the
 *  gauss points, so modifications of the returned array are forbidden.
 *  @param[in] order  An integer of the set {2,4,8,16,32}.
 *  @return A constant array of doubles containing the gauss points to the
 *          given order. The array contains exactly "order" elements.
 */
const double* getGaussPoints(int order);

/**
 *  The function returns quadrature weights for the Gaussian quadrature rule of
 *  given order. It does not create a copy of the gauss weights, so
 *  modifications of that array are forbidden.
 *  @param[in] order  An integer of the set {2,4,8,16,32}.
 *  @return A constant array of doubles containing the gauss weights to the
 *          given order. The array contains exactly "order" elements.
 */
const double* getGaussWeights(int order);
}  // namespace mysolution

#endif
