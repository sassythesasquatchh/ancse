///////////////////////////////////////////////////////////////////////////////
/// \file evaluateW.cpp
/// \brief This file contains the function evaluateW that evaluates the
///        hypersingular integral operator W on any number of evaluation points
///        on the boundary \f$Gamma\f$.
///
///  This file contains only the implementation. For extensive documentation
///  consult the corresponding header-file.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>
#include <iostream>
#include "constants.hpp"
#include "geometry.hpp"
#include "doubleLayerPotential.hpp"
extern "C" {
#include "gaussQuadrature.h"
}
#include "evaluateW.hpp"


void evaluateW(Eigen::VectorXd& Wx, const BoundaryMesh& mesh, const Eigen::VectorXd &gh,
               const Eigen::MatrixXd &x, const Eigen::MatrixXd &n_x,
               double eta)
{
  int nX = x.rows();
  int nE = mesh.numElements();
  int nC = mesh.numVertices();
  // Initialize output vector
  Wx.resize(nX);  Wx.setZero();

  // Get quadrature points and weights
  const double* qp = getGaussPoints(GAUSS_ORDER);
  const double* qw = getGaussWeights(GAUSS_ORDER);

  for(int j = 0; j < nX; ++j){
    // save current point for readibility
    const Eigen::Vector2d& xj = x.row(j);

    /* Find the element which contains the current evaluation point.
     * Warn the user if
     *   - No such element is found.
     *   - More than one such element is found.
     *   - The evaluation point is very close to an end point of that
     *      element.
     * Either of these cases results in unreliable results and should
     * be avoided.
     */
    
    int found_element = 0; /* Number of elements containing the evaluation point. */
    int j_element = 0;     /* Id of the last found element. */
    double s = 0;          /* Ratio dist(x,a)/length of element.*/
    double gh_xj = 0.;     /* Function value of gh in evaluation point. */
    double g_0m1, g_1m1, g_0m2, g_1m2, g_2m2;

    for (int i = 0; i < nE; ++i){
      // get vertices indices and coordinates for Ei=[a,b]
      int aidx = mesh.getElementVertex(i,0);
      int bidx = mesh.getElementVertex(i,1);
      const Eigen::Vector2d& a = mesh.getVertex(aidx);
      const Eigen::Vector2d& b = mesh.getVertex(bidx);
      double lengthEi = (b-a).norm();

      double dist_x_Ei = distancePointToSegment(xj, a, b);

      if(dist_x_Ei < 5e-2 * EPS){
        found_element += 1;
        j_element = i;

        if ( (a-xj).norm() < EPS || (b-xj).norm() < EPS){
          fprintf(stderr, "Warning (evaluateW): The evaluation point "
            "(%g,%g) is too close to the end point of an element. The "
            "%d-th value of the return vector is probably wrong!\n",
            xj[0], xj[1], j);
        }

        double dist_xa = (xj-a).norm();
        /* Length of vector from the evaluation point to the start point
         * of the element containing the evaluation point. */
        assert(dist_xa >= 0. && dist_xa <= lengthEi);

        s = (2. * dist_xa/lengthEi) - 1.;
        assert(s < 1. && s > -1.);

        assert(aidx >= 0 && aidx < nC);
        assert(bidx >= 0 && bidx < nC);

        double gh_a = gh(aidx); /* Value of gh in start point of element. */
        double gh_b = gh(bidx); /* Value of gh in end point of element. */

        gh_xj = gh_a + (dist_xa / lengthEi) * (gh_b-gh_a);

        break;
      }
    } //end for loop over elements

    if (found_element == 0){
      fprintf(stderr, "Warning (evaluateW): The evaluation point "
        "(%g,%g) is not even close to an element! The %d-th entry of "
        "the return vector is set to NaN!\n", xj[0], xj[1], j);
      Wx[j] = NAN;
      continue;
    }

    if (found_element > 1){
      fprintf(stderr, "Warning (evaluateW): The evaluation point "
        "(%g,%g) is very close to %d elements. This either means that "
        "the evaluation point is very close to the end point of an element "
        "or that the input mesh is invalid. The %d-th entry of the return "
        "vector is set to NaN!\n", xj[0], xj[1], found_element, j);
      Wx[j] = NAN;
      continue;
    }

    // Computation of Wx element-wise
    double W_0=0, W_1=0;
    for(int i = 0; i < nE; ++i){
      // get vertices indices and coordinates for Ei=[a,b]
      int aidx = mesh.getElementVertex(i,0);
      int bidx = mesh.getElementVertex(i,1);
      const Eigen::Vector2d& a = mesh.getVertex(aidx);
      const Eigen::Vector2d& b = mesh.getVertex(bidx);
      // auxiliary vectors
      const Eigen::Vector2d& mp = 0.5*(a+b);
      const Eigen::Vector2d& vc  = (b-a);
      const Eigen::Vector2d& u = -0.5*vc;
      const Eigen::Vector2d& v = xj - mp;
      // auxiliary constants
      double u_sqnorm = u.squaredNorm();
      double v_sqnorm = v.squaredNorm();
      double udotv    = u.dot(v);
      double delta    = u_sqnorm*v_sqnorm - udotv*udotv;
      double upvsq    = u_sqnorm + 2*udotv + v_sqnorm;
      double umvsq    = u_sqnorm - 2*udotv + v_sqnorm;
      double sqrdelta = 2 * sqrt(delta);

      // Determine the normal vector
      Eigen::Vector2d n_vc;
      n_vc[0] =  0.5*vc[1] / sqrt(u_sqnorm);
      n_vc[1] = -0.5*vc[0] / sqrt(u_sqnorm);

      // Compute distance xj to Ei
      double dist_x_Ei = distancePointToSegment(xj, a, b);
      // Check whether the evaluation point is on the current element.
      if (dist_x_Ei < 5e-2 * EPS){
       // Compute p.v. of a locally regularized integral
        double lengthEi = sqrt(4.0*u_sqnorm);
        double alpha = (gh(bidx)-gh(aidx))/2.;
        assert(j_element == i);

        Wx(j) += 1./(lengthEi*M_PI)*( 2.0*gh_xj/(1.0-s*s)
                                      -alpha*log(fabs((s-1.0)/(s+1.0))) );
      }
      else{
        // Check if pair (x,T_j) is admissible
        if ((sqrt(4*u_sqnorm) > eta*dist_x_Ei) || eta==-1){
          // Compute the integrals analytically
          /* g_0^(-1)*/
          g_0m1 = dlp(0,u,v);
          g_1m1 = dlp(1,u,v);
  
          /* g_0^(-2)*/
          if(sqrdelta*sqrdelta > 4.0*EPS * u_sqnorm * v_sqnorm){
            g_0m2 = (  2*(u_sqnorm  + udotv)/(upvsq)
                     - 2*(-u_sqnorm + udotv)/(umvsq)
                     + 2*u_sqnorm*g_0m1 )/(sqrdelta*sqrdelta);
          }
          else{
            g_0m2 = 2*(u_sqnorm + 3*v_sqnorm)/(3*std::pow(v_sqnorm-u_sqnorm,3));
          }
  
          /* g_1^(-2)*/
          g_1m2 = ( -1/upvsq + 1/umvsq - 2*udotv*g_0m2 )/(2.*u_sqnorm);

          /* g_2^(-2)*/
          g_2m2 = ( g_0m1 - 2*udotv*g_1m2 - v_sqnorm*g_0m2 )/u_sqnorm;
  
          // Compute W_0
          W_0 = ( n_vc.dot(n_x.row(j))*g_0m1
                -2*( u.dot(n_x.row(j))*n_vc.dot(u)*g_1m2
                     + (v.dot(n_x.row(j))*n_vc.dot(u)
                        +u.dot(n_x.row(j))*n_vc.dot(v))*g_1m2
                     + u.dot(n_x.row(j))*n_vc.dot(v)*g_0m2 ) )/M_PI;
  
          // Compute W_1
          W_1 = ( n_vc.dot(n_x.row(j))*g_1m1
                 -2*( u.dot(n_x.row(j))*n_vc.dot(u)*g_2m2
                     + (v.dot(n_x.row(j))*n_vc.dot(u)
                        +u.dot(n_x.row(j))*n_vc.dot(v))*g_2m2
                     + v.dot(n_x.row(j))*n_vc.dot(v)*g_1m2 ) )/M_PI;
        }
        else{
          // Compute the integrals via quadrature rule*/
          W_0=0;
          W_1=0;
          for(int k=0;k<GAUSS_ORDER;++k){
            double denom = u_sqnorm*qp[k]*qp[k] + 2*udotv*qp[k] + v_sqnorm;
            W_0 += qw[k]*( n_vc.dot(n_x.row(j))/denom
                           -2*(n_x.row(j).dot(u*qp[k]+v)*n_vc.dot(u*qp[k]+v))
                           /(denom*denom) )/M_PI;
  
            W_1 += qw[k]*qp[k]*( n_vc.dot(n_x.row(j))/denom
                             -2*(n_x.row(j).dot(u*qp[k]+v)*n_vc.dot(u*qp[k]+v))
                                 /(denom*denom) )/M_PI;
          } // end quadrature for loop
        }

        // Putting the integrals together to obtain W(phi)(x)
        Wx(j) += 0.25*sqrt(u_sqnorm)*W_1*( gh(aidx) - gh(bidx) )
                -0.25*sqrt(u_sqnorm)*W_0*( gh(bidx) + gh(aidx) );

      }

    } // end for loop elements
  } //end for loop evaluation points

}

