///////////////////////////////////////////////////////////////////////////////
/// \file geometry.hpp
/// \brief This file provides functions to calculate the distance between a
///  segment and a point and two segments.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _GEOMETRY_HPP_GUARD_
#define _GEOMETRY_HPP_GUARD_

#include <cmath>
#include <cstdio>
#include <cassert>

#include <Eigen/Dense>
#include "constants.hpp"


/**
 *  This function treats the segment as a parametrized vector where the
 *  parameter t varies from 0 to 1. It finds the value of t that minimizes the
 *  distance from the point to the line.
 *
 *  If t is between 0.0 and 1.0, then the closest point lies on the segment.
 *  Otherwise the closest point is one of the segments end points.
 *
 *  @param[in] a,b,p  2d vectors containing the coordinates of these points.
 *  @return  Distance between the point p and the segment [a,b].
 */
double distancePointToSegment(const Eigen::Vector2d& p,
			      const Eigen::Vector2d& a,
			      const Eigen::Vector2d& b);

/**
 *  This function returns 1, in case the angle ABX is in (0,pi),i.e. the
 *  segment AB can be transformed into the segment AX using a rotation in
 *  mathematical positive sense (counter-clock wise).
 *
 *  It returns -1 in the other case, ABX in (pi,2pi). In case, a, b, and x are
 *  collinear, ABX is either 0 or  pi.
 *
 *  The function returns 0, in case x is in [a,b].
 *  If it's not, either -1 is returned, in case that x is closer to a or 1 is
 *  returned if x is closer to b.
 *
 *  This function is called by distanceSegmentToSegment multiple times to check
 *  whether two segments intersect.
 *
 *  @param[in] a,b,x  2d vectors containing the coordinates of these points.
 *  @return  An integer being either -1, 0 or 1.
*/
int ccw(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
	const Eigen::Vector2d& x);

/**
 *  The distance between two segments is either zero if they intersect or the
 *  distance from one of the segments' end points to the other line. "distance"
 *  first checks if the segments intersect. In case they do, the function
 *  returns 0.
 *  Otherwise it calculates the distance between each end point and the
 *  respective other segment and returns the minimum.
 *
 *  @param[in] a,b,c,d  2d vectors containing the coordinates of these points.
 *  @return  Distance between the segments [a,b] and [c,d].
 */
double distanceSegmentToSegment(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
				const Eigen::Vector2d& c, const Eigen::Vector2d& d);

/**
 *  Computes the unit normal of [a,b]
 *
 *  @param[in] a,b 2d vectors containing the coordinates of these points.
 *  @return  Normal vector.
 */
Eigen::Vector2d unitNormal(const Eigen::Vector2d& a, const Eigen::Vector2d& b);


/**
 *  Computes the 2D cross product of the vectors a and b
 *
 *  @param[in] a,b 2d vectors.
 *  @return 2D cross product
 */
double CrossProd2d(const Eigen::Vector2d& a, const Eigen::Vector2d& b);

#endif
