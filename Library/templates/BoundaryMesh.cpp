///////////////////////////////////////////////////////////////////////////////
/// \file BoundaryMesh.cpp
/// \brief Implementation of BoundaryMesh class
///////////////////////////////////////////////////////////////////////////////

#include "BoundaryMesh.hpp"

#include <iomanip>
#include <utility>

namespace mysolution {
//------------------------------------------------------------------------------
std::pair<BoundaryMesh, Eigen::MatrixXi> BoundaryMesh::uniformRefine() const {
  Eigen::MatrixXd new_coordinates =
      Eigen::MatrixXd::Zero(coordinates_.rows() * 2, coordinates_.cols());
  Eigen::MatrixXi new_elements =
      Eigen::MatrixXi::Zero(elements_.rows() * 2, elements_.cols());

  Eigen::MatrixXi father2son(elements_.rows(), 2);

  // TODO: Populate the matrices for new coordinates, elements and the father to son relation

  BoundaryMesh newmesh(new_coordinates, new_elements);

  return std::make_pair(newmesh, father2son);
};

//------------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh(const coord_matrix_t& coords,
                           const elem_matrix_t& elems) {
  // TODO: implement the constructor
}

//------------------------------------------------------------------------------
int BoundaryMesh::numVertices() const {  // TODO: Implement the function
  return -1;
};

//------------------------------------------------------------------------------
int BoundaryMesh::numElements() const {  // TODO: Implement the function
  return -1;
};

//------------------------------------------------------------------------------
BoundaryMesh::coord_matrix_t BoundaryMesh::getMeshVertices() const {
  assert(isInitialized_);
  // TODO: implement the function

  BoundaryMesh::coord_matrix_t out;
  return out;
};

//------------------------------------------------------------------------------
BoundaryMesh::elem_matrix_t BoundaryMesh::getMeshElements() const {
  assert(isInitialized_);
  // TODO: implement the function

  BoundaryMesh::elem_matrix_t out;
  return out;
};

//------------------------------------------------------------------------------
Eigen::Vector2d BoundaryMesh::getVertex(int i) const {
  assert(isInitialized_);
  // TODO: implement the function
  return Eigen::Vector2d(0, 0);
};

//------------------------------------------------------------------------------
std::pair<Eigen::Vector2d, Eigen::Vector2d> BoundaryMesh::getElementVertices(
    int i) const {
  assert(isInitialized_);
  // TODO: implement the function

  return std::make_pair(Eigen::Vector2d(0, 0), Eigen::Vector2d(0, 0));
};

//------------------------------------------------------------------------------
int BoundaryMesh::getElementVertex(int i, int j) const {
  assert(isInitialized_);
  // TODO: implement the function

  return -1;
};

//------------------------------------------------------------------------------
void BoundaryMesh::loadMeshFromFile(const std::string& filename) {
  readData<coord_matrix_t>(filename + "_coordinates.dat", coordinates_);
  readData<elem_matrix_t>(filename + "_elements.dat", elements_);
  // elements file has indexing starting from 1. Fix it!
  elements_ =
      elements_ - Eigen::MatrixXi::Ones(elements_.rows(), elements_.cols());

  // Print mesh information
  std::cout << std::string(80, '=') << std::endl
            << std::string(27, ' ') << " READING MESH FROM FILE \n"
            << std::string(80, '=') << std::endl;
  std::cout << "Input file : " << filename << std::endl;
  std::cout << "Created " << coordinates_.rows() << " vertices "
            << "with coordinates :\n"
            << coordinates_ << std::endl
            << std::endl;
  std::cout << "Created " << elements_.rows() << " elements " << ": \n"
            << elements_ << std::endl
            << std::endl;
  std::cout << std::string(80, '=') << std::endl;

  isInitialized_ = 1;
};

//------------------------------------------------------------------------------
void BoundaryMesh::writeMeshToFile(const std::string& filename) {
  std::ofstream out_coords(filename + "_coordinates.dat");
  out_coords << std::setprecision(18) << coordinates_;
  out_coords.close();

  std::ofstream out_els(filename + "_elements.dat");
  out_els << std::setprecision(18) << elements_;
  out_els.close();
};

//------------------------------------------------------------------------------
template <typename T>
void BoundaryMesh::readData(const std::string& filename, T& data) {
  std::ifstream indata(filename);
  // check state
  if (!indata) {
    std::cout << "Could not open file '" << filename << " \n"
              << "File does not exist!" << std::endl;
    exit(-1);
  }

  std::vector<typename T::Scalar> values;
  std::string line;
  int rows = 0;
  // read every line from the stream
  while (std::getline(indata, line)) {
    std::stringstream dataStream(line);
    std::string dataCell;
    // read every cell from the line that is seperated by space
    // and put it into the vector or strings
    while (std::getline(dataStream, dataCell, '	')) {
      values.push_back(std::stod(dataCell));
    }
    rows++;
  }
  int cols = values.size() / rows;

  data.resize(rows, cols);
  data =
      Eigen::Map<const Eigen::Matrix<typename T::Scalar, T::RowsAtCompileTime,
                                     T::ColsAtCompileTime, Eigen::RowMajor>>(
          values.data(), rows, cols);
};
}  // namespace mysolution
