///////////////////////////////////////////////////////////////////////////////
/// \file newtonPotential.cpp
/// \brief This file provides functions to compute the newton potential.
///
///  This file is part of the HILBERT program package for the numerical
///  solution of the Laplace equation with mixed boundary conditions by use of
///  BEM in 2D.
///
///  C++ adaptation for ANCSE17 of HILBERT V3.1 TUWien 2009-2013
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include "constants.hpp"
#include "newtonPotential.hpp"


//-----------------------------------------------------------------------------
void computeN(Eigen::MatrixXd& N,
              const Eigen::Matrix<double, Eigen::Dynamic, 2>& coordinates,
              const Eigen::Matrix<int   , Eigen::Dynamic, 2>& elements,
              const Eigen::Matrix<double, Eigen::Dynamic, 2>& vertices,
              const Eigen::Matrix<int   , Eigen::Dynamic, 3>& triangles,
              double eta)
{
  int nT = triangles.rows();
  int nE = elements.rows();
  N.resize(nE,nT);
  Eigen::MatrixXd nodes(3,2);

  for( int j=0;j<nT; ++j) { /*  running over triangles */
    nodes << vertices.row(triangles(j,0)),
             vertices.row(triangles(j,1)),
             vertices.row(triangles(j,2));

    for( int k=0; k<nE; ++k) { /* running over boundary elements */
      const Eigen::Vector2d& a = coordinates.row(elements(k,0));
      const Eigen::Vector2d& b = coordinates.row(elements(k,1));

      N(k,j) = computeNkj(a,b,nodes,eta);
    }
  }
}


//-----------------------------------------------------------------------------
double computeNkj(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                  const Eigen::Matrix<double,3,2>& nodes, double eta)
{
  double distET = distanceSegmentToSegment(a, b, nodes.row(0), nodes.row(1));
  double tmp = distanceSegmentToSegment(a, b, nodes.row(0), nodes.row(2));

  if (tmp < distET)
    distET = tmp;

  tmp = distanceSegmentToSegment(a, b, nodes.row(1), nodes.row(2));

  if (tmp < distET)
    distET = tmp;

  double diam = (a-b).squaredNorm();

  if (distET*eta >= diam) 	/* Semi-analytic */
     return computeNkjSemiAnalyticSegment(a,b,nodes);

  else 				/* analytic */
      return computeNkjAnalytic(a,b,nodes);
}


//-----------------------------------------------------------------------------
double computeNkjSemiAnalyticSegment(const Eigen::Vector2d& a,
                                     const Eigen::Vector2d& b,
                                     const Eigen::Matrix<double,3,2> &nodes)
{
  const double *gauss_points = getGaussPoints(GAUSS_ORDER);
  const double *gauss_weights = getGaussWeights(GAUSS_ORDER);

  double val = 0.;
  for (int k=0; k<GAUSS_ORDER; ++k){
    Eigen::Vector2d tmp = 0.5*(a+b + gauss_points[k]*(b-a));

    val += gauss_weights[k]*newtonPotential(nodes,tmp);
  }

  return -val*(a-b).norm()*0.25/ M_PI;
}


//-----------------------------------------------------------------------------
double computeNkjSemiAnalyticTriangle(const Eigen::Vector2d& a,
                                      const Eigen::Vector2d& b,
                                      const Eigen::Matrix<double,3,2>& nodes)
{
  const int nr_gauss_points = 7;
  const double* gp_x = getGaussPointsT(nr_gauss_points,0);
  const double* gp_y = getGaussPointsT(nr_gauss_points,1);
  const double* gauss_weights = getGaussWeightsT(nr_gauss_points);

  const Eigen::Vector2d& n1 = nodes.row(0);
  const Eigen::Vector2d& n2 = nodes.row(1);
  const Eigen::Vector2d& n3 = nodes.row(2);

  double detFT = fabs((n2[0]-n1[0])*(n3[1]-n1[1])-(n3[0]-n1[0])*(n2[1]-n1[1]));
  double tmp=0.;

  for (int i=0; i<nr_gauss_points; i++){
    Eigen::Vector2d aux1 = 0.5*(b-a);
    Eigen::Vector2d aux2 = 0.5*(b+a)-(n2-n1)*gp_x[i] - (n3-n1)*gp_y[i] - n1;
    tmp += gauss_weights[i]*slp(0,aux1,aux2);
  }

  return -0.125*(a-b).norm()*detFT*tmp/M_PI;
}


//-----------------------------------------------------------------------------
double computeNkjAnalytic(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                          const Eigen::Matrix<double,3,2>& nodes)
{
  const Eigen::Vector2d& n1 = nodes.row(0);
  const Eigen::Vector2d& n2 = nodes.row(1);
  const Eigen::Vector2d& n3 = nodes.row(2);

  double detFK   = fabs((n2[0]-n1[0])*(n3[1]-n1[1])-(n3[0]-n1[0])*(n2[1]-n1[1]));
  double volK    = 0.5*detFK;
  double lengthT = (a-b).norm();

  /* aux holds (n1-n3)/2 -n1 + mj */
  const Eigen::Vector2d& aux  = (a+b)/2. - (n1+n3)/2.;
  /* aux1 holds mj - (n2+n3)/2 */
  const Eigen::Vector2d& aux1 = (a+b)/2. - (n2+n3)/2.;
  /* aux2 holds mj - (n1+n3)/2 */
  const Eigen::Vector2d& aux2 = (a+b)/2. - (n1+n3)/2.;

  Eigen::VectorXd coeff(6);
  const Eigen::Vector2d& u = n1-n2;
  const Eigen::Vector2d& naux = (n1-n2)/(2*u.squaredNorm());
  coeff(0) = naux.dot(aux) + 0.25 ;
  coeff(1) = naux.dot((n1-n3)/2.) - 0.25 ;
  coeff(2) = naux.dot((b-a)/2.);
  coeff(3) = naux.dot(aux);
  coeff(4) = coeff(1) + 0.25 ;
  coeff(5) = coeff(2);

  Eigen::VectorXd J(6);
  J(0) = doubleSlp(0, 0, (b-a)/2., (n2-n3)/2., aux1);
  J(1) = doubleSlp(0, 1, (b-a)/2., (n2-n3)/2., aux1);
  J(2) = doubleSlp(1, 0, (b-a)/2., (n2-n3)/2., aux1);
  J(3) = doubleSlp(0, 0, (b-a)/2., (n2-n3)/2., aux2);
  J(4) = doubleSlp(0, 1, (b-a)/2., (n2-n3)/2., aux2);
  J(5) = doubleSlp(1, 0, (b-a)/2., (n2-n3)/2., aux2);

  double atanInt = integrateAtanInt(a,b,nodes);
  double ret = lengthT*(coeff.dot(J) - 1) + 0.5*atanInt/u.squaredNorm();

  return -0.5*volK*ret/M_PI;
}


//-----------------------------------------------------------------------------
double newtonPotential(const Eigen::Matrix<double,3,2>& nodes,
                       const Eigen::Vector2d& x)
{
  const Eigen::Vector2d& n1 = nodes.row(0);
  const Eigen::Vector2d& n2 = nodes.row(1);
  const Eigen::Vector2d& n3 = nodes.row(2);

  double volT = 0.5*fabs((n2[0]-n1[0])*(n3[1]-n1[1])-(n3[0]-n1[0])*(n2[1]-n1[1]));
  const Eigen::Vector2d& u = n1 - n2;
  const Eigen::Vector2d& v = n1 - n3;
  const Eigen::Vector2d& w = x  - n1;

  Eigen::VectorXd slpA = slpIterative(1, (v-u)/2., w + (u+v)/2.);
  Eigen::VectorXd slpB = slpIterative(1, v/2. , w + v/2.);

  double integralI = 0.25*(slpA[0]-slpA[1]);

  double integralII = u.dot(w)*slpA[0] + 0.5*u.dot(v)*(slpA[0]+slpA[1]);

  double integralIII = u.dot(w + v/2.)*slpB[0] + 0.5*u.dot(v)*slpB[1];

  double atanint = evalAtanInt(nodes,x);

  double coeff = volT/(2.*u.squaredNorm());
  double result = volT*integralI + coeff*integralII -
                  coeff*integralIII - volT + coeff*atanint;

  return result;
}


//-----------------------------------------------------------------------------
double evalAtanInt(const Eigen::Matrix<double,3,2>& nodes,
                   const Eigen::Vector2d& x)
{
  /* compute u and v */
  const Eigen::Vector2d& u = nodes.row(0) - nodes.row(1);
  const Eigen::Vector2d& v = nodes.row(0) - nodes.row(2);
  /* compute w */
  const Eigen::Vector2d& w = x  - nodes.row(0).transpose();

  const double* gauss_point = getGaussPoints(GAUSS_ORDER);
  const double* gauss_wht   = getGaussWeights(GAUSS_ORDER);

  double val=0.;
  for (int j= 0; j<GAUSS_ORDER; ++j ){
    /* transformation of quadrature points from [-1,1] to [0,1] */
    double xi = (gauss_point[j]+1)*0.5;

    double a = u.squaredNorm();
    double b = 2*(u.dot(w)+xi*u.dot(v));
    double c = w.squaredNorm()+xi*xi*v.squaredNorm()+2*xi*w.dot(v);
    double valInner = innerAtanInt(a,b,c,xi);
    val += gauss_wht[j]*valInner;
  }

  return val;
}


//-----------------------------------------------------------------------------
double innerAtanInt(double a, double b, double c, double xi)
{
  double delta = 4*a*c-b*b;
  double deltaRoot = sqrt(delta);
  double s = (2*a*(1-xi)+b)/deltaRoot;
  double t = -b/deltaRoot;
  double st = s*t;
  double valInner = 0.;

  if(delta>EPS*1*a*c) {
    if ( st<1-EPS )
      valInner = deltaRoot*atan( (s+t)/(1-st));
    else {
      if ( st>1+EPS) {
        if (s < 0)
      valInner = deltaRoot*(atan( (s+t)/(1-st)) - M_PI);
        else
      valInner = deltaRoot*(atan( (s+t)/(1-st)) + M_PI);
      }
      else {
        valInner = deltaRoot*(atan(s) + atan(t));
      }
    }
  }
  return valInner;
}


//-----------------------------------------------------------------------------
double integrateAtanInt(const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                        const Eigen::Matrix<double,3,2>& nodes)
{
  const double* gauss_point = getGaussPoints(GAUSS_ORDER);
  const double* gauss_wht   = getGaussWeights(GAUSS_ORDER);

  double val = 0.;
  for (int k=0; k<GAUSS_ORDER; ++k ){
    /* transformation of quadrature points from [-1,1] to [a,b] */
    const Eigen::Vector2d& sx = ((1-gauss_point[k])*a+(1+gauss_point[k])*b)*0.5;
    double valOuter = evalAtanInt(nodes,sx);

    val += gauss_wht[k]*valOuter;
  }

  val = val*0.5*(a-b).norm();
  return val;
}

