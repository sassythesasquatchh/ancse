///////////////////////////////////////////////////////////////////////////////
/// \file geometry.cpp
/// \brief This file provides functions to calculate the distance between a
///  segment and a point and two segments.
///
///  This file contains only the implementation. For documentation see
///  geometry.hpp
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#include "geometry.hpp"

namespace mysolution {
//-----------------------------------------------------------------------------
double distancePointToSegment(const Eigen::Vector2d& p,
                              const Eigen::Vector2d& a,
                              const Eigen::Vector2d& b) {
  double t = 0.;
  Eigen::Vector2d dab = b - a; /* Distance between a and b.  */
  Eigen::Vector2d r;           /* Distance between [a,b] and p. */
  r.setZero();

  if (fabs(dab[0]) <= EPS && fabs(dab[1]) <= EPS) {
    /* Case that the segment is shorter than sqrt(2)*EPS.*/
    r = p - a;
  } else {
    t = (p - a).dot(dab) / (dab.squaredNorm());

    if (t <= EPS) {
      r = p - a;
    } else if ((1 - t) <= EPS) {
      r = p - b;
    } else {
      r = p - a - t * dab;
    }
  }

  return r.norm();
}

//-----------------------------------------------------------------------------
int ccw(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
        const Eigen::Vector2d& x) {
  assert(a[0] != b[0] || a[1] != b[1]);

  Eigen::Vector2d dab = b - a;
  Eigen::Vector2d dax = x - a;

  /* dab[1] / dab[0] < dax[1] / dax[0] <=> dab[0] * dax[1] > dax[0] * dab[1] */
  /* dab[0] / dab[1] > dax[0] / dax[1] <=> the angle ABX is in (0,pi). */
  if (dab[0] * dax[1] > dax[0] * dab[1]) {
    return +1;
  } else if (dab[0] * dax[1] < dax[0] * dab[1]) {
    /* Contrary, this is equivalent to the fact that ABX is in (pi,2pi). */
    return -1;
  } else { /* Otherwise, a, b, x are colinear. */
    if (dab[0] * dax[0] < 0 || dab[1] * dax[1] < 0) { /* x is closer to a */
      return -1;
    } else if (dab.squaredNorm() < dax.squaredNorm()) {
      return 1; /* x is closer to b */
    } else {    /* x is in [a,b]. */
      return 0;
    }
  }
}

//-----------------------------------------------------------------------------
double distanceSegmentToSegment(const Eigen::Vector2d& a,
                                const Eigen::Vector2d& b,
                                const Eigen::Vector2d& c,
                                const Eigen::Vector2d& d) {
  double test = 0.;
  double best = (a - c).norm();

  /* The only possibility that two segments intersect: */
  if (ccw(a, b, c) * ccw(a, b, d) <= 0 && ccw(c, d, a) * ccw(c, d, b) <= 0) {
    best = 0.;
  } else { /* If segements do not intersect: */
    /* Calculate the distance from a to segment [c,d]: */
    test = distancePointToSegment(a, c, d);
    if ((test - best) <= EPS) {
      best = test;
    }
    /* Calculate the distance from b to segment [c,d]: */
    test = distancePointToSegment(b, c, d);
    if ((test - best) <= EPS) {
      best = test;
    }
    /* Calculate the distance from c to segment [a,b]: */
    test = distancePointToSegment(c, a, b);
    if ((test - best) <= EPS) {
      best = test;
    }
    /* Calculate the distance from d to segment [a,b]: */
    test = distancePointToSegment(d, a, b);
    if ((test - best) <= EPS) {
      best = test;
    }
  }

  return best;
}

//-----------------------------------------------------------------------------
Eigen::Vector2d unitNormal(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
  Eigen::Vector2d n;
  n << (b[1] - a[1]), -(b[0] - a[0]);
  n /= (b - a).norm();
  return n;
}

//-----------------------------------------------------------------------------
double CrossProd2d(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
  return a[0] * b[1] - a[1] * b[0];
}
}  // namespace mysolution
