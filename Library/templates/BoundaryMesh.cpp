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
  coordinates_ = coords;
  elements_ = elems;
  isInitialized_ = 1;
}

//------------------------------------------------------------------------------
int BoundaryMesh::numVertices() const {  // TODO: Implement the function
  return coordinates_.rows();
};

//------------------------------------------------------------------------------
int BoundaryMesh::numElements() const {  // TODO: Implement the function
  return elements_.rows();
};

//------------------------------------------------------------------------------
BoundaryMesh::coord_matrix_t BoundaryMesh::getMeshVertices() const {
  assert(isInitialized_);
  // TODO: implement the function

  return coordinates_;
};

//------------------------------------------------------------------------------
BoundaryMesh::elem_matrix_t BoundaryMesh::getMeshElements() const {
  assert(isInitialized_);
  // TODO: implement the function

  return elements_;
};

//------------------------------------------------------------------------------
Eigen::Vector2d BoundaryMesh::getVertex(int i) const {
  assert(isInitialized_);
  // Ensure that the vertex index is not out of bounds
  assert(i < coordinates_.rows());
  // TODO: implement the function
  // Return the i-th row of the coordinates matrix
  return coordinates_.row(i);
};

//------------------------------------------------------------------------------
std::pair<Eigen::Vector2d, Eigen::Vector2d> BoundaryMesh::getElementVertices(
    int i) const {
  assert(isInitialized_);
  // TODO: implement the function
  // Ensure that the element index is not out of bounds
  assert(i < elements_.rows());

  // Get the first and second vertex idxs (row idxs) of the i-th element
  auto first_vertex_idx = elements_(i, 0);
  auto second_vertex_idx = elements_(i, 1);
  // Get the corresponding vertices
  return std::make_pair(coordinates_.row(first_vertex_idx), coordinates_.row(second_vertex_idx));
};

//------------------------------------------------------------------------------
int BoundaryMesh::getElementVertex(int i, int j) const {
  assert(isInitialized_);
  // TODO: implement the function
  // Ensure that the element index is not out of bounds
  assert(i < elements_.rows());
  // Each element has only 2 vertices
  assert(j < 2);
  // Return the corresponding vertex
  return _elements(i, j);
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
