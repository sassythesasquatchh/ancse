///////////////////////////////////////////////////////////////////////////////
/// \file newtonPotential.hpp
/// \brief This file provides functions to compute the newton potential.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _NEWTONPOTENTIAL_HPPmysolution
#define _NEWTONPOTENTIAL_HPPmysolution

#include <cmath>

#include "geometry.hpp"
#include "singleLayerPotential.hpp"

namespace mysolution {
/**
 *  The function computeNkj() is called for the computation of the
 *  corresponding matrix entry.
 *
 *  @param[out] N  (nE x 3*nT) matrix. nE is the number of elements of the
 *                 boundary \f$\Gamma\f$ and nT the number of triangles of the
 *                 domain mesh of \f$\Omega\f$.
 *  @param[in] coordinates  (nC x 2) matrix containing the coordinates of the
 *                          vertices of the boundary mesh.
 *  @param[in] elements  (nE x 2) matrix containing the indices of the vertices
 *                       corresponding to each element of the boundary mesh.
 *  @param[in] vertices  (nV x 2) matrix containing the coordinates of the
 *                       vertices of the domain mesh.
 *  @param[in] triangles  (nT x 3) matrix containing the indices of the vertices
 *                        conforming the triangles of the domain mesh.
 *  @param[in] eta  Admissibility constant. It is greater or equal than 0.
 */
void computeN(Eigen::MatrixXd& N,
              const Eigen::Matrix<double, Eigen::Dynamic, 2>& coordinates,
              const Eigen::Matrix<int, Eigen::Dynamic, 2>& elements,
              const Eigen::Matrix<double, Eigen::Dynamic, 2>& vertices,
              const Eigen::Matrix<int, Eigen::Dynamic, 3>& triangles,
              double eta);

/**
 *  First, the line segment Ek = [a,b] and the triangle Tj are tested for
 *  admissibility. If they are admissible, the entry N_kj is computed
 *  semi-analytically. Otherwise, the entry is computed analytically.
 *
 *  @param[in] a,b  End points of the line segment Ek.
 *  @param[in] nodes  Corner points of the domain mesh element Tj.
 *  @param[in] eta  Admissibility parameter.
 *  @return The function returns the value of the discrete Newton potential,
 *   \f[ -\frac{1}{2\pi} \int_{Ek} \int_{Tj} \log{ \vert x-y \vert} dy ds_x.\f]
 */
double computeNkj(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                  const Eigen::Matrix<double, 3, 2>& nodes, double eta);

/**
 *  The entry N_kj is computed SEMI-ANALYTIC. This means the function uses the
 *  function newtonPotential() to compute the inner integral analytically. A
 *  Gauss quadrature rule is used to perform the outer integration over the
 *  line segment.
 *
 *  @param[in] a,b  End points of the boundary mesh element Ek
 *  @param[in] nodes  Corner points of the domain mesh element Tj.
 *  @return  The function returns the value
 *   \f[ -\frac{1}{2\pi} \int_{Ek} \int_{Tj} \log{ \vert x-y \vert} dy ds_x.\f]
 */
double computeNkjSemiAnalyticSegment(const Eigen::Vector2d& a,
                                     const Eigen::Vector2d& b,
                                     const Eigen::Matrix<double, 3, 2>& nodes);

/**
 *  The entry N_kj is computed "ANALYTIC". This means that the two integrals
 *  are computed analytically up to an integral of a smooth function. This
 *  remaining integral is evaluated by use of a Gaussian quadrature rule. This,
 *  clearly is not a fully analytic computation. Besides the above mentioned
 *  integral over a smooth function, we have to compute double-integrals of the
 *  type \f$ \int_{-1}^1\int_{-1}^1 s^k t^l \log{\vert su+tv+w \vert} dt ds,\f$
 *  which is done by calling the function doubleSlp().
 *
 *  @param[in] a, b are the points A,B of the line segment Ek.
 *  @param[in] nodes  Corner points of the domain mesh element Tj.
 *  @return  The function returns the entry N_{k,j} of the discrete Newton
 *  potential \f[ N_{kj} = -\frac{1}{2\pi} \int_{Ek} \int_{Tj} \log{ \vert x-y
 *  \vert} dy ds_x. \f]
 */
double computeNkjAnalytic(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                          const Eigen::Matrix<double, 3, 2>& nodes);

/**
 * The result is computed analytically up to a quadrature of a smooth function.
 *
 *  @param[in] nodes  Corner points of the domain mesh element T, given in
 *                    counter clockwise order.
 *  @param[in] x  Point in R^2.
 *  @return  The function returns the newton-potential \f$ \int_{T} \log{
 *           \vert x-t \vert} dt \f$
 */
double newtonPotential(const Eigen::Matrix<double, 3, 2>& nodes,
                       const Eigen::Vector2d& x);

/**
 *  We use a Gauss rule for the outer integral whereas the inner integral is
 *  computed analytically using the function innerAtanInt().
 *
 *  @param[in] nodes  Corner points of the domain mesh element T, given in
 *                    counter clockwise order.
 *  @param[in] x  Point in R^2.
 *  @return  The function returns the value of the integral
 *           \f[ \int_0^1 \int_0^{1-\xi} \frac{\Delta}{c + \eta b + \eta^2}
 *            d \eta d\xi dx \f]
 *  where \f$ \Delta, a, b, c \f$ are defined as in the documentation for the
 *  function innerAtanInt(). This integral arises as part of the computation of
 *  integrateAtanInt(), but it is also used from within the function
 *  newtonPotential().
 */
double evalAtanInt(const Eigen::Matrix<double, 3, 2>& nodes,
                   const Eigen::Vector2d& x);

/**
 *  Computes the integral and returns its value.
 *
 *  @param[in] xi  Value between 0 and 1.
 *  @param[in] a  Defined by \f$ a := \vert u \vert^2 \f$, where
 *                \f$ u:=-n_2+n_1, v:=-n_3+n_1, w := x-n_1. \f$.
 *  @param[in] b  Defined by \f$ b := 2u \cdot (w+\xi v) \f$
 *  @param[in] c  Defined by \f$ c := \vert w+\xi v \vert^2 \f$
 *
 *  @return  The function returns the value of the integral
 *           \f[ \int_0^\xi \Delta / (c+\eta b+\eta^2 a) d \eta. \f]
 *           where \f$\Delta := 4*a*c-b^2.\f$
 */
double innerAtanInt(double a, double b, double c, double xi);

/**
 *  Computes the integral by transforming it to an integral over the reference
 *  interval (-1,1) and the usage of a gauss quadrature formula. The inner
 *  integrals are evaluated using the evalAtanInt() function.
 *
 *  @param[in] a,b  End points of the boundary mesh element Ej
 *  @param[in] nodes  Corner points of the domain mesh element Tk.
 *  @return  The function returns the value of the integral
 */
double integrateAtanInt(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                        const Eigen::Matrix<double, 3, 2>& nodes);
}  // namespace mysolution

#endif
