///////////////////////////////////////////////////////////////////////////////
/// \file singleLayerPotential.cpp
/// \brief This file provides functions to compute single layer potential (slp)
///  type integrals.
///
///  This file contains only the implementation. For detailed documentation see
///  singleLayerPotential.hpp
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////

#include "singleLayerPotential.hpp"

#include <cstdlib>

#include "constants.hpp"
#include "gaussQuadrature.hpp"

namespace mysolution {
//-----------------------------------------------------------------------------
double slp(int k, const Eigen::Vector2d& u, const Eigen::Vector2d& v) {
  Eigen::VectorXd tmp = slpIterative(k, u, v);
  return tmp[k];
}

//------------------------------------------------------------------------------
// Computes the integrals I^(0)_k from Def. 2.2 in MAI08
// using a recursion formula given in Lemma 2.2
// returns all values for powers up to k in a vector
//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_1 */
Eigen::VectorXd slpIterative(int k, const Eigen::Vector2d& u,
                             const Eigen::Vector2d& v) {
  double a = u.squaredNorm();  // \cob{$\alpha = \N{\Bu}^2$}
  double b = 2. * u.dot(v);    // \cob{$\beta = 2\Bu\cdot\Bv$}
  double c = v.squaredNorm();  // \cob{$\gamma = \N{\Bv}^2$}
  double D = 0.;               // discriminant
  Eigen::VectorXd val(k + 1);  // return values

  // Ensure one non-zero argument vector
  double tmp = 4 * a * c - b * b;
  assert(fabs(u[0]) > EPS || fabs(u[1]) > EPS || fabs(v[0]) > EPS ||
         fabs(v[1]) > EPS);
  // By Cauchy-Schwarz inequality tmp >= 0
  assert(tmp >= -fabs(EPS * 4 * a * c));

  // Numerically sound way of testing if discriminant = 0
  if (tmp > EPS * 4 * a * c)
    D = sqrt(tmp);
  else
    D = 0.;

  // TODO: populate the output vecctor val, taking care of various possible values for k, D etc.

  return val;
}

//------------------------------------------------------------------------------
// Section 3 of MAI08,
// Used only for evaluation of Newton potential
//-----------------------------------------------------------------------------
double doubleSlp(int k, int l, const Eigen::Vector2d& u,
                 const Eigen::Vector2d& v, const Eigen::Vector2d& w) {
  double output = 0.;
  Eigen::VectorXd memTableUWpv, memTableUWmv, memTableVWpu, memTableVWmu;
  double normUSq = u.squaredNorm();
  double normVSq = v.squaredNorm();
  double normWSq = w.squaredNorm();
  double detUV = CrossProd2d(u, v);

  if (normUSq < EPS && normVSq < EPS) {
    if (k % 2 == 0) {
      if (l % 2 == 0) {
        if (normWSq < EPS)
          return 0;
        else
          return (double)2. / (double)((k + 1) * (l + 1)) * log(normWSq);
      }
    }
    return 0.;
  } else if (normUSq < EPS) {
    if (k % 2 == 0)
      return (double)1. / (k + 1) * slp(l, v, w);
    else
      return 0.;
  } else if (normVSq < EPS) {
    if (l % 2 == 0)
      return (double)1. / (l + 1) * slp(k, u, w);
    else
      return 0.;
  }

  if (fabs(detUV) < EPS * sqrt(normUSq * normVSq)) { /* u,v parallel */
    double mu = 0.;

    memTableUWpv = slpIterative(k, u, w + v);
    memTableUWmv = slpIterative(k, u, w - v);
    memTableVWpu = slpIterative(k + l + 1, v, w + u);
    memTableVWmu = slpIterative(k + l + 1, v, w - u);

    if (fabs(u[0]) < fabs(u[1]))
      mu = v[1] / u[1];
    else
      mu = v[0] / u[0];

    output = memTableUWpv[0] - mu * memTableVWpu[k + l + 1] +
             mu * memTableVWmu[k + l + 1];
    if ((k + l) % 2 == 0)
      output += memTableUWmv[0];
    else
      output -= memTableUWmv[0];
    output /= (2 * (k + l + 1));

    for (int i = 1; i <= k; ++i) {
      output *= 2 * i * mu;
      output += memTableUWpv[i] - mu * memTableVWpu[l + k + 1 - i];

      if ((k + l - i) % 2 == 0)
        output += memTableUWmv[i];
      else
        output -= memTableUWmv[i];

      if (i % 2 == 0)
        output += mu * memTableVWmu[l + k + 1 - i];
      else
        output -= mu * memTableVWmu[l + k + 1 - i];

      output /= (2 * (l + k + 1 - i));
    }
  } else {
    double mu1 = 0., mu2 = 0.;
    Eigen::VectorXd tmp(l + 1);

    if (fabs(w[0] + v[0] - u[0]) < fabs(w[0]) * EPS &&
        fabs(w[1] + v[1] - u[1]) < fabs(w[1]) * EPS) {
      memTableVWmu.resize(l + 1);
      memTableVWmu.setZero();
      memTableUWpv.resize(k + 1);
      memTableUWpv.setZero();
    } else {
      memTableVWmu = slpIterative(l, v, w - u);
      memTableUWpv = slpIterative(k, u, w + v);
    }

    if (fabs(w[0] + u[0] - v[0]) < fabs(w[0]) * EPS &&
        fabs(w[1] + u[1] - v[1]) < fabs(w[1]) * EPS) {
      memTableVWpu.resize(l + 1);
      memTableVWpu.setZero();
      memTableUWmv.resize(k + 1);
      memTableUWmv.setZero();
    } else {
      memTableVWpu = slpIterative(l, v, w + u);
      memTableUWmv = slpIterative(k, u, w - v);
    }

    mu1 = CrossProd2d(w, v) / detUV;
    mu2 = CrossProd2d(u, w) / detUV;

    tmp[0] = -2 + ((mu1 + 1) * memTableVWpu[0] - (mu1 - 1) * memTableVWmu[0] +
                   (mu2 + 1) * memTableUWpv[0] - (mu2 - 1) * memTableUWmv[0]) *
                      0.25;

    for (int i = 1; i <= l; ++i) {
      tmp[i] =
          0.5 * ((mu1 + 1) * memTableVWpu[i] - (mu1 - 1) * memTableVWmu[i] +
                 (mu2 + 1) * memTableUWpv[0]) -
          i * mu2 * tmp[i - 1];
      if (i % 2 == 0) {
        tmp[i] -= (double)4. / (double)(i + 1);
        tmp[i] -= 0.5 * (mu2 - 1) * memTableUWmv[0];
      } else
        tmp[i] += 0.5 * (mu2 - 1) * memTableUWmv[0];

      tmp[i] /= (i + 2);
    }

    for (int i = 1; i <= k; ++i) {
      tmp[0] =
          0.5 * ((mu1 + 1) * memTableVWpu[0] + (mu2 + 1) * memTableUWpv[i] -
                 (mu2 - 1) * memTableUWmv[i]) -
          i * mu1 * tmp[0];
      if (i % 2 == 0) {
        tmp[0] -= (double)4. / (double)(i + 1);
        tmp[0] -= 0.5 * (mu1 - 1) * memTableVWmu[0];
      } else {
        tmp[0] += 0.5 * (mu1 - 1) * memTableVWmu[0];
      }

      tmp[0] /= (i + 2);

      for (int j = 1; j <= l; ++j) {
        tmp[j] *= -i * mu1;
        tmp[j] -= j * mu2 * tmp[j - 1];
        tmp[j] +=
            0.5 * ((mu1 + 1) * memTableVWpu[j] + (mu2 + 1) * memTableUWpv[i]);
        if (i % 2 == 0) {
          if (j % 2 == 0) {
            tmp[j] -= (double)4. / (double)((i + 1) * (j + 1));
          }
          tmp[j] -= 0.5 * (mu1 - 1) * memTableVWmu[j];
        } else {
          tmp[j] += 0.5 * (mu1 - 1) * memTableVWmu[j];
        }

        if (j % 2 == 0) {
          tmp[j] -= 0.5 * (mu2 - 1) * memTableUWmv[i];
        } else {
          tmp[j] += 0.5 * (mu2 - 1) * memTableUWmv[i];
        }

        tmp[j] /= (i + j + 2);
      }
    }

    output = tmp[l];
  }

  return output;
}

//------------------------------------------------------------------------------
// Computation of entries of single layer Galerkin matrix for piecewise
// constantb trial and test functions.
//-----------------------------------------------------------------------------
double computeVij(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                  const Eigen::Vector2d& c, const Eigen::Vector2d& d,
                  double eta) {
  double hi = (b - a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d - c).squaredNorm(); /* hj = norm(d-c)^2 */

  return sqrt(hi * hj) * computeWij(a, b, c, d, eta);
}

//-----------------------------------------------------------------------------
double computeWij(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                  const Eigen::Vector2d& c, const Eigen::Vector2d& d,
                  double eta) {
  int swap = 0;
  double hi = (b - a).squaredNorm(); /* hi = norm(b-a)^2 */
  double hj = (d - c).squaredNorm(); /* hj = norm(d-c)^2 */
  double tmp = 0.;

  /* For stability reasons, we guarantee   hj <= hi   to ensure that *
   * outer integration is over smaller domain. This is done by       *
   * swapping Ej and Ei if necessary.                                */
  if (hj > hi) {
    swap = 1;
    std::swap(hi, hj);
  }

  if (eta == 0) { /* compute all matrix entries analytically */
    if (swap == 1)
      return computeWijAnalytic(c, d, a, b);
    else
      return computeWijAnalytic(a, b, c, d);
  } else { /* compute admissible matrix entries semi-analytically */
    if (distanceSegmentToSegment(a, b, c, d) > eta * sqrt(hj)) {
      if (swap == 1)
        return computeWijSemianalytic(c, d, a, b);
      else
        return computeWijSemianalytic(a, b, c, d);
    } else {
      if (swap == 1)
        return computeWijAnalytic(c, d, a, b);
      else
        return computeWijAnalytic(a, b, c, d);
    }
  }
}

//------------------------------------------------------------------------------
// Analytic integration of logarithmic kernel over two straight panels
//-----------------------------------------------------------------------------
/* SAM_LISTING_BEGIN_2 */
double computeWijAnalytic(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                          const Eigen::Vector2d& c, const Eigen::Vector2d& d) {
  double hi =
      (b - a).squaredNorm();  // length$^2$ of first panel \cob{$[\Ba,\Bb]$}
  double hj =
      (d - c).squaredNorm();  // lendth$^2$ of second panel \cob{$[\Bc,\Bd]$}
  double val = 0.;

  // There hold different recursion formulae when the panels
  // are parallel (det = 0) or not

  // TODO: implement the evaluation of det
  double det = 0;

  if (fabs(det) <=
      EPS * sqrt(hi * hj)) {  // parallel panels, Case II \Label[line]{Wij:2}

    // TODO:
    // Use function slp() to evaluate the required integral

  }  // \Label[line]{Wij:2a}
  else {  // \cob{$\Bx$} and \cob{$\By$} linearly independent, Case I \Label[line]{Wij:1}

    // TODO:
    // Use function slp() to evaluate the required integral
  }
  return val;
}
/* SAM_LISTING_END_2 */

//-----------------------------------------------------------------------------
double computeWijSemianalytic(const Eigen::Vector2d& a,
                              const Eigen::Vector2d& b,
                              const Eigen::Vector2d& c,
                              const Eigen::Vector2d& d) {
  double val = 0;
  const double* gauss_point = getGaussPoints(GAUSS_ORDER);
  const double* gauss_wht = getGaussWeights(GAUSS_ORDER);

  Eigen::Vector2d u = 0.5 * (a - b);

  for (int k = 0; k < GAUSS_ORDER; ++k) {
    /* transformation of quadrature nodes from [-1,1] to [a,b] */
    Eigen::Vector2d sx =
        ((1 - gauss_point[k]) * c + (1 + gauss_point[k]) * d) * 0.5;
    Eigen::Vector2d v = sx - 0.5 * (a + b);

    /* inner product wht*func(sx) */
    val += gauss_wht[k] * slp(0, u, v);
  }

  return -0.0625 * val / M_PI; /* = - 1/(16*M_PI) * int(log |.|^2) */
}
}  // namespace mysolution
