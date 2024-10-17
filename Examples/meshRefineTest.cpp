#include <Eigen/Sparse>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

#include "BoundaryMesh.hpp"

int main(int, char **) {
  BoundaryMesh mesh("Lshape");

  Eigen::MatrixXd coords = mesh.getMeshVertices();
  coords /=4;
  Eigen::MatrixXi elems = mesh.getMeshElements();

  auto pair_output = mesh.uniformRefine();

  BoundaryMesh refinedMesh = pair_output.first;
  Eigen::MatrixXi father2son = pair_output.second;

  Eigen::MatrixXd coords_ref = refinedMesh.getMeshVertices();
  coords_ref /= 4;
  Eigen::MatrixXi elems_ref = refinedMesh.getMeshElements();

  std::cout << "Mesh coordinated = \n" << coords << std::endl;
  std::cout << "Mesh elements = \n" << elems << std::endl;

  std::cout << "Refined mesh coordinates = \n" << coords_ref << std::endl;
  std::cout << "Refined mesh elements = \n" << elems_ref << std::endl;

  std::cout << "Father2son = \n" << father2son << std::endl;

  return 0;
}
