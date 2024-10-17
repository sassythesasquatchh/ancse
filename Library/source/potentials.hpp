#ifndef POTENTIALSHPP
#define POTENTIALSHPP

#include <Eigen/Dense>

//////////////////////////////////////////////////////////////////////////////
// lambda function for the potential u1
auto u1 = [](Eigen::Vector2d x) { return x(0); };

// Lambda function for the gradient of u1
auto gradu1 = [](Eigen::Vector2d x) {
  Eigen::Vector2d out(1, 0);
  return out;
};

//////////////////////////////////////////////////////////////////////////////
// Lambda function for potential u2
auto u2 = [](Eigen::Vector2d x) { return std::cosh(x(0)) * std::cos(x(1)); };

// Lambda function for gradient of potential u2
auto gradu2 = [](Eigen::Vector2d x) {
  Eigen::Vector2d out;
  out << std::sinh(x(0)) * std::cos(x(1)), -std::cosh(x(0)) * std::sin(x(1));
  return out;
};

//////////////////////////////////////////////////////////////////////////////
// Lambda function for potential u(x,y) = Re{Z^{2/3}}
auto u = [](Eigen::Vector2d x) {
  double alpha = 2. / 3.;
  // Getting the polar coordinates
  double r = x.norm();
  double phi = std::atan2(x(1), x(0));
  return std::pow(r, alpha) * std::cos(alpha * phi);
};

// Gradient of the above potential.
auto gradu = [](Eigen::Vector2d x) {
  Eigen::Vector2d out, rhat, thetahat;
  double alpha = 2. / 3.;

  // Polar coordinates
  double r = x.norm();
  double phi = std::atan2(x(1), x(0));

  // Basis in polar coordinates
  rhat << std::cos(phi), std::sin(phi);
  thetahat << -std::sin(phi), std::cos(phi);

  // Gradient in polar coordinates
  out = alpha * std::pow(r, alpha - 1) *
        (std::cos(alpha * phi) * rhat - std::sin(alpha * phi) * thetahat);
  return out;
};

#endif