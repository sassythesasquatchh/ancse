#include <Eigen/Sparse>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

#include "meshReader.hpp"
#include "source/geometry.hpp"
#include "source/evaluateV.hpp"
#include "source/evaluateW.hpp"
#include "source/evaluateK.hpp"
#include "source/evaluateKadj.hpp"

//! Sparse Matrix type. Makes using this type easier.
typedef Eigen::SparseMatrix<double> SparseMatrix;

//! Vector type
typedef Eigen::VectorXd Vector;





int main(int, char**) {

  Eigen::MatrixXd coords;
  Eigen::MatrixXi els;

  readMesh("Lshape", coords, els);

  std::cout << "x" << std::endl;
  Eigen::MatrixXd x(3,2);
  Eigen::MatrixXd temp(4,2);
  temp << coords.row(0), coords.row(1), coords.row(2), coords.row(5);
  Eigen::MatrixXd nx(3,2);
  for(int i=0; i<3; i++){
    x.row(i) = (temp.row(i) + temp.row(i+1))/2.;
    nx(i,0) = x(i,1)/(x.row(i)).norm(), -x(i,0)/(x.row(i)).norm();
  }
  std::cout << x << std::endl;
  
  std::cout << "phi" << std::endl;
  Eigen::VectorXd phi(coords.rows());
  phi.setOnes();
  
  Eigen::VectorXd Vdom;
  evaluateV(Vdom, coords, els, phi, x, 0.001);

  Eigen::VectorXd Wdom;
  evaluateW(Wdom, coords, els, phi, x, nx, 0.001);

  Eigen::VectorXd Kdom;
  evaluateK(Kdom, coords, els, phi, x,  0.001);

  Eigen::VectorXd Kadom;
  evaluateKadj(Kadom, coords, els, phi, x, nx,  0.001);

  const std::string fnameV =  "Vdom_Lshape_Hilbert.dat";
  std::ofstream outV( fnameV.c_str() );
  outV << std::setprecision(18) << Vdom; 
  outV.close( );

  const std::string fnameW =  "Wdom_Lshape_Hilbert.dat";
  std::ofstream outW( fnameW.c_str() );
  outW << std::setprecision(18) << Wdom; 
  outW.close( );

  const std::string fnameK =  "Kdom_Lshape_Hilbert.dat";
  std::ofstream outK( fnameK.c_str() );
  outK << std::setprecision(18) << Kdom; 
  outK.close( );

    const std::string fnameKa =  "Kadom_Lshape_Hilbert.dat";
  std::ofstream outKa( fnameKa.c_str() );
  outKa << std::setprecision(18) << Kadom; 
  outKa.close( );
  
  std::cout << " Done " << std::endl;

  // loaded them in matlab and compared with Hilbert original implementation.
  // the norm of the difference is zero :)
  
  
  
}




