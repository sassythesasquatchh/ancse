///////////////////////////////////////////////////////////////////////////////
/// \file BoundaryMesh.cpp
/// \brief Implementation of BoundaryMesh class
///////////////////////////////////////////////////////////////////////////////

#include "BoundaryMesh.hpp"

#include <iomanip>
#include <utility>

//------------------------------------------------------------------------------
std::pair<BoundaryMesh, Eigen::MatrixXi> BoundaryMesh::uniformRefine() const {
  Eigen::MatrixXd new_coordinates =
      Eigen::MatrixXd::Zero(coordinates_.rows() * 2, coordinates_.cols());
  Eigen::MatrixXi new_elements =
      Eigen::MatrixXi::Zero(elements_.rows() * 2, elements_.cols());

  // populate coordinates with the coordinates of the parent mesh
  // Preserves the index of the original coordinates
  new_coordinates.block(0, 0, coordinates_.rows(), coordinates_.cols()) =
      coordinates_;

  Eigen::MatrixXi father2son(elements_.rows(), 2);

  // Loop over parent elements
  for (int i = 0; i < elements_.rows(); ++i) {
    // Get the vertex indices for the ith element
    int aidx = getElementVertex(i, 0);
    int bidx = getElementVertex(i, 1);
    const Eigen::Vector2d& a = getVertex(aidx);
    const Eigen::Vector2d& b = getVertex(bidx);

    Eigen::Vector2d mid = 0.5 * (a + b);

    // Adding the mid-point to the list of coordinates
    new_coordinates.row(elements_.rows() + i) = mid.transpose();

    // Adding the new elements
    new_elements(2 * i, 0) = aidx;
    new_elements(2 * i, 1) = elements_.rows() + i;
    new_elements(2 * i + 1, 0) = elements_.rows() + i;
    new_elements(2 * i + 1, 1) = bidx;

    // father to son relation
    father2son(i, 0) = 2 * i;
    father2son(i, 1) = 2 * i + 1;
  }

  BoundaryMesh newmesh(new_coordinates, new_elements);

  return std::make_pair(newmesh, father2son);
};

//------------------------------------------------------------------------------

BoundaryMesh::BoundaryMesh(const coord_matrix_t& coords,
                           const elem_matrix_t& elems) {
  coordinates_ = coords;
  elements_ = elems;
  isInitialized_ = 1;
}

//------------------------------------------------------------------------------
int BoundaryMesh::numVertices() const { return coordinates_.rows(); };

//------------------------------------------------------------------------------
int BoundaryMesh::numElements() const { return elements_.rows(); };

//------------------------------------------------------------------------------
BoundaryMesh::coord_matrix_t BoundaryMesh::getMeshVertices() const {
  assert(isInitialized_);

  return coordinates_;
};

//------------------------------------------------------------------------------
BoundaryMesh::elem_matrix_t BoundaryMesh::getMeshElements() const {
  assert(isInitialized_);

  return elements_;
};

//------------------------------------------------------------------------------
Eigen::Vector2d BoundaryMesh::getVertex(int i) const {
  //std::cout << "function called with i = " << i << std::endl;
  assert(isInitialized_);
  assert(i < coordinates_.rows());

  return coordinates_.row(i);
};

//------------------------------------------------------------------------------
std::pair<Eigen::Vector2d, Eigen::Vector2d> BoundaryMesh::getElementVertices(
    int i) const {
  assert(isInitialized_);
  assert(i < elements_.rows());

  return std::make_pair<Eigen::Vector2d, Eigen::Vector2d>(
      coordinates_.row(elements_(i, 0)), coordinates_.row(elements_(i, 1)));
};

//------------------------------------------------------------------------------
int BoundaryMesh::getElementVertex(int i, int j) const {
  assert(isInitialized_);
  assert(i < elements_.rows());
  assert(j < 2);

  return elements_(i, j);
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
