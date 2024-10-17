///////////////////////////////////////////////////////////////////////////////
/// \file doubleLayerPotential.cpp
/// \brief This file provides functions to compute the Galerkin Matrix K given
///  by
///  \f[ K_{ij} = -\frac{1}{2\pi}\int_{Ei} \int_{supp \phi_j}\frac{<y-x,n>}{
///                |y-x|^2} \phi_j(y)ds_y ds_x \f]
///
///  This file contains only the implementation. For detailed documentation see
///  doubleLayerPotential.hpp
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
////////////////////////////////////////////////////////////////////////////////
#include "doubleLayerPotential.hpp"

#include <cassert>
#include <cmath>
#include <iostream>

#include "constants.hpp"
#include "gaussQuadrature.hpp"

namespace mysolution {
//------------------------------------------------------------------------------
// See MAI08, Section 2
//------------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
double dlp(int k, const Eigen::Vector2d& p, const Eigen::Vector2d& q) {
  // The full recursion is not implemented
  assert(k <= 2 && (k >= 0));

  double a = p.squaredNorm();    // a = <p,p>
  double b = 2 * p.dot(q);       // b = 2 <p,q>
  double c = q.squaredNorm();    // c = <q,q>
  double D = 4 * a * c - b * b;  // Discriminant
  double root_D = 0.;
  double G0 = 0., G1 = 0.;

  assert(D >= -EPS * 4 * a * c);  // In exact arithmetic, D >= 0
  if (D > EPS * 4 * a * c) {
    root_D = sqrt(D);
  } else {
    D = 0.0;
  }

  // TODO: Implement the function to give the correct output

  return 0;
}
/* SAM_LISTING_END_1 */

//-----------------------------------------------------------------------------
void computeKij(double* I0, double* I1, double eta, const Eigen::Vector2d& a,
                const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                const Eigen::Vector2d& d) {
  int swap = 0;
  double hi = (b - a).squaredNorm();  // hi = norm(b-a) squared
  double hj = (d - c).squaredNorm();  // hj = norm(d-c) squared

  if ((hi - hj) / hi > EPS) {
    swap = 1; /* swap Ej and Ei */
    std::swap(hi, hj);
  }

  if (eta == 0) /* fully analytic computation */
  {
    if (swap == 0)
      computeKijAnalytic(I0, I1, a, b, c, d);
    else
      computeKijSwappedAnalytic(I0, I1, a, b, c, d);
  } else /* fully analytic or semianalytic */
  {
    if (distanceSegmentToSegment(a, b, c, d) * eta >= sqrt(hi)) {
      if (swap == 0)
        computeKijSemianalytic(I0, I1, a, b, c, d);
      else
        computeKijSwappedSemianalytic(I0, I1, a, b, c, d);
    } else {
      if (swap == 0)
        computeKijAnalytic(I0, I1, a, b, c, d);
      else
        computeKijSwappedAnalytic(I0, I1, a, b, c, d);
    }
  }
}

//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_2 */
void computeKijAnalytic(double* I0, double* I1, const Eigen::Vector2d& a,
                        const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                        const Eigen::Vector2d& d) {
  double hi = (b - a).squaredNorm();     // hi = norm(b-a) squared
  double hj = (d - c).squaredNorm();     // hj = norm(d-c) squared
  Eigen::Vector2d n = unitNormal(c, d);  // normal vector

  // TODO: complete the function to compute the values for the integrals I0 and I1
  // as mentioned in the lecture document.
}
/* SAM_LISTING_END_2 */

//-----------------------------------------------------------------------------
void computeKijSwappedAnalytic(double* I0, double* I1, const Eigen::Vector2d& a,
                               const Eigen::Vector2d& b,
                               const Eigen::Vector2d& c,
                               const Eigen::Vector2d& d) {
  double hi = (b - a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d - c).squaredNorm(); /* hj = norm(d-c)^2 */

  Eigen::Vector2d n = unitNormal(c, d); /* normal vector */
  Eigen::Vector2d u = d - c;
  Eigen::Vector2d v = a - b;
  Eigen::Vector2d w = c + d - a - b;
  Eigen::Vector2d wpu = w + u;
  Eigen::Vector2d wmu = w - u;
  Eigen::Vector2d wpv = w + v;
  Eigen::Vector2d wmv = w - v;

  double dot_v_n = v.dot(n); /* dot_v_n=<v,n> */
  double dot_w_n = w.dot(n);
  double dot_wpv_n = (wpv).dot(n); /* dot_wpv_n=<w+v,n> */
  double dot_wmv_n = (wmv).dot(n);

  double det = CrossProd2d(u, v);

  double lambda = 0.0, mu = 0.0;
  if (fabs(det) <= EPS * sqrt(hi * hj)) /* u,v linearly dependent */
  {
    if (fabs(u[0]) > fabs(u[1]))
      mu = v[0] / u[0];
    else
      mu = v[1] / u[1];

    *I0 = dot_w_n * (dlp(0, u, wpv) + dlp(0, u, wmv) +
                     mu * (dlp(1, v, wmu) - dlp(1, v, wpu)));
    *I1 = dot_w_n * (dlp(1, u, wpv) + dlp(1, u, wmv) +
                     mu * (-dlp(1, v, wmu) - dlp(1, v, wpu) +
                           0.5 * (dlp(0, u, wpv) - dlp(0, u, wmv) +
                                  mu * (dlp(2, v, wmu) - dlp(2, v, wpu)))));

  } else /* u,v linearly independent */
  {
    if (a[0] == d[0] && a[1] == d[1]) {
      *I0 = 2 * (dot_v_n * dlp(1, v, wmu) + dot_w_n * dlp(0, v, wmu) +
                 dot_wpv_n * dlp(0, u, wpv));
      *I1 = dot_wpv_n * dlp(1, u, wpv) - dot_v_n * dlp(1, v, wmu) -
            dot_w_n * dlp(0, v, wmu) + 0.5 * (*I0);
    } else if (b[0] == c[0] && b[1] == c[1]) {
      *I0 = 2 * (dot_v_n * dlp(1, v, wpu) + dot_w_n * dlp(0, v, wpu) +
                 dot_wmv_n * dlp(0, u, wmv));
      *I1 = dot_v_n * dlp(1, v, wpu) + dot_w_n * dlp(0, v, wpu) +
            dot_wmv_n * dlp(1, u, wmv) - 0.5 * (*I0);
    } else {
      mu = CrossProd2d(w, v) / det;
      lambda = CrossProd2d(u, w) / det;

      *I0 = (mu + 1) * (dot_v_n * dlp(1, v, wpu) + dot_w_n * dlp(0, v, wpu)) -
            (mu - 1) * (dot_v_n * dlp(1, v, wmu) + dot_w_n * dlp(0, v, wmu)) +
            (lambda + 1) * dot_wpv_n * dlp(0, u, wpv) -
            (lambda - 1) * dot_wmv_n * dlp(0, u, wmv);
      *I1 = 0.5 *
            ((mu + 1) * (dot_v_n * dlp(1, v, wpu) + dot_w_n * dlp(0, v, wpu)) +
             (mu - 1) * (dot_v_n * dlp(1, v, wmu) + dot_w_n * dlp(0, v, wmu)) +
             (lambda + 1) * dot_wpv_n * dlp(1, u, wpv) -
             (lambda - 1) * dot_wmv_n * dlp(1, u, wmv) - mu * (*I0));
    }
  }
  *I0 *= -0.125 * sqrt(hi * hj) / M_PI;
  *I1 *= -0.125 * sqrt(hi * hj) / M_PI;
}

//-----------------------------------------------------------------------------
void computeKijSemianalytic(double* I0, double* I1, const Eigen::Vector2d& a,
                            const Eigen::Vector2d& b, const Eigen::Vector2d& c,
                            const Eigen::Vector2d& d) {
  double hi = (b - a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d - c).squaredNorm(); /* hj = norm(d-c)^2 */

  /* 16-point Gaussian quadrature on [-1,1] */
  const int order = 16;
  const double* gauss_weight = getGaussWeights(order);
  const double* gauss_point = getGaussPoints(order);

  Eigen::Vector2d n = unitNormal(c, d);

  Eigen::Vector2d u = a - b;
  Eigen::Vector2d v = d - c;
  Eigen::Vector2d w = c + d - a - b;

  double I0tmp = 0.0;
  double I1tmp = 0.0;
  for (int i = 0; i < order; ++i) {
    Eigen::Vector2d z = gauss_point[i] * u + w;
    double dot_z_n = z.dot(n);

    I0tmp += gauss_weight[i] * dot_z_n * dlp(0, v, z);
    I1tmp += gauss_weight[i] * dot_z_n * dlp(1, v, z);
  }

  *I0 = -0.125 * sqrt(hi * hj) * I0tmp / M_PI;
  *I1 = -0.125 * sqrt(hi * hj) * I1tmp / M_PI;
}

//-----------------------------------------------------------------------------
void computeKijSwappedSemianalytic(double* I0, double* I1,
                                   const Eigen::Vector2d& a,
                                   const Eigen::Vector2d& b,
                                   const Eigen::Vector2d& c,
                                   const Eigen::Vector2d& d) {
  double hi = (b - a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d - c).squaredNorm(); /* hj = norm(d-c)^2 */

  /* 16-point Gaussian quadrature on [-1,1] */
  const int order = 16;
  const double* gauss_weight = getGaussWeights(order);
  const double* gauss_point = getGaussPoints(order);

  Eigen::Vector2d n = unitNormal(c, d); /* normal vector */

  Eigen::Vector2d u = d - c;
  Eigen::Vector2d v = a - b;
  Eigen::Vector2d w = c + d - a - b;

  double dot_v_n = v.dot(n);
  double dot_w_n = w.dot(n);

  double I0tmp = 0.0;
  double I1tmp = 0.0;
  for (int i = 0; i < order; ++i) {
    Eigen::Vector2d z = gauss_point[i] * u + w;

    I0tmp +=
        gauss_weight[i] * (dot_v_n * dlp(1, v, z) + dot_w_n * dlp(0, v, z));
    I1tmp += gauss_weight[i] * gauss_point[i] *
             (dot_v_n * dlp(1, v, z) + dot_w_n * dlp(0, v, z));
  }

  *I0 = -0.125 * sqrt(hi * hj) * I0tmp / M_PI;
  *I1 = -0.125 * sqrt(hi * hj) * I1tmp / M_PI;
}
}  // namespace mysolution
