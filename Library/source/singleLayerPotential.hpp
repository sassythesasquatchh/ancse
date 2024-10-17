///////////////////////////////////////////////////////////////////////////////
/// \file singleLayerPotential.hpp
/// \brief This file provides functions to compute single layer potential (slp)
///  type integrals.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////
#ifndef _SINGLELAYERPOTENTIAL_HPP
#define _SINGLELAYERPOTENTIAL_HPP

#include <cmath>
#include <cassert>
#include "geometry.hpp"
extern "C" {
#include "gaussQuadrature.h"
}


/**
 *  This function calculates the integral as stated above. Internally, it uses 
 *  the slpIterative-function and returns the last element of the returned 
 *  vector   
 *  
 *  @param[in] k   Non-negative integer. Power of s. 
 *  @param[in] u,v 2d vectors. 
 *  @return Value of the integral \f$ \int_{-1}^{+1} s^k \log{\vert s u+v 
 *          \vert^2} ds. \f$  
 */
double slp(int k, const Eigen::Vector2d& u, const Eigen::Vector2d& v);


/**
 *  @param[in] k   Non-negative integer. Power of s (and size of output vector). 
 *  @param[in] u,v 2d vectors. 
 *  @return (k+1) vector r that is given by \f$ r_i = \int_{-1}^1 s^i \log{\vert 
 *          s u+v \vert^2} ds. \f$ 
 */
Eigen::VectorXd slpIterative(int k, const Eigen::Vector2d& u,
			     const Eigen::Vector2d& v);


/**
 *  This function calculates the integral as stated above in a purely analytical 
 *  way.
 *
 *  @param[in] k,l   Non-negative integers.
 *  @param[in] u,v,w 2d vectors. 
 *  @return Value of the integral \f$ \int_{-1}^1 \int_{-1}^1 s^k t^l \log{\vert
 *           su+tv+w \vert} dt ds \f$ 
 */
double doubleSlp(int k, int l, const Eigen::Vector2d& u,
		 const Eigen::Vector2d& v, const Eigen::Vector2d& w);


/**
 *  This function uses computeWij() and multiplies the result  with |Ei|*|Ej| 
 *  to obtain the result.        
 *
 *  @param[in] a,b,c,d  2d vectors with the coordinates of these points.
 *  @param[in] eta      Admissible constant. It is a number between 0 and 1.
 *  @return Value of the Galerkin integral \f$ -\frac{1}{2\pi} \int_{Ej} 
 *          \int_{Ei} \log{ \vert x-y \vert } dsy dsx \f$, 
 *          where Ei = [a,b] and Ej = [c,d].
 */
double computeVij(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
		  const Eigen::Vector2d& c, const Eigen::Vector2d& d,
		  double eta);


/**
 *  This function checks whether the boundary elements Ei and Ej are admissible,
 *  this means that \f$ dist(Ei, Ej) > \eta \min(diam(Ei), diam(Ej)) \f$.
 *  In case, Ei and Ej are admissible computeWijSemianalytic() is called. 
 *  Otherwise computeWijAnalytic() is called.
 *
 *  @param[in] a,b,c,d  2d vectors with the coordinates of these points.
 *  @param[in] eta      Admissible constant. It is a number between 0 and 1.            
 *  @return Value of the Galerkin integral \f$ -\frac{1}{2 \pi} \frac{1}{|Ei|} 
 *           \frac{1}{|Ej|} \int_{Ej} \int_{Ei} \log{\vert x-y \vert} dsy dsx \f$
 *          where Ei = [a,b] and Ej = [c,d].                                   
 */
double computeWij(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
		  const Eigen::Vector2d& c, const Eigen::Vector2d& d, double eta);


/**
 *  This function calculates the double integral in a purely analytical way. It 
 *  is called when Ei and Ej are not admissible. 
 *
 *  @param[in] a,b,c,d  2d vectors with the coordinates of these points.
 *  @return Value of the Galerkin integral \f$ -\frac{1}{2 \pi} \frac{1}{|Ei|}  
 *          \frac{1}{|Ej|} \int_{Ej} \int_{Ei} \log{\vert x-y \vert} dsy dsx \f$
 *          whereas Ei = [a,b] and Ej = [c,d].      
 */
double computeWijAnalytic(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
			  const Eigen::Vector2d& c, const Eigen::Vector2d& d);


/**
 *  This function replaces the outer integral with gauss quadrature to avoid 
 *  cancellation effects. It must not get called whenever the segments Ei and Ej
 *   are not admissible.  
 *
 *  @param[in] a,b,c,d  2d vectors with the coordinates of these points. The  
 *                      segments Ei=[a,b], Ej=[c,d] are required to be admissible.
 *  @return Value of the Galerkin integral \f$ -\frac{1}{2 \pi} \frac{1}{|Ei|}  
 *          \frac{1}{|Ej|} \int_{Ej} \int_{Ei} \log{\vert x-y \vert} dsy dsx \f$
 *          whereas Ei = [a,b] and Ej = [c,d].      
 */
double computeWijSemianalytic(const Eigen::Vector2d& a,
			      const Eigen::Vector2d& b,
			      const Eigen::Vector2d& c,
			      const Eigen::Vector2d& d);

#endif
