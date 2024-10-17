///////////////////////////////////////////////////////////////////////////////
/// \file BoundaryMesh.hpp
/// \brief This class provides functions to read the meshes used in HILBERT and
///        handles the boundary mesh data.
///
/// The data structure follows the description in Section 3.3 of the HILBERT
/// documentation:
///
/// The set of nodes \f$ \{z_1,...,z_N \}\f$ of the mesh M_h with elements \f$
/// \{E_1,..,E_N \} \f$ is represented by the (Nx2)-Matrix coordinates. The j-th
/// row of __coordinates__ stores the coordinates of the j-th node \f$z_j=(x_j,
/// y_j)\f$.
///
/// If the boundary \f$\Gamma\f$ is not split into several parts, its mesh M_h
/// is represented by the (Nx2)-Matrix __elements__. The i-th boundary element
/// \f$E_i=[z_j,z_k]\f$ is stored as
///
///                            elements.row(i)=[j, k],
///
/// where the nodes are given in *counterclockwise* order,i.e. the parametrization
/// of the boundary element \f$ E_i \subset \Gamma \f$ is mathematically positive.
/// In other words, the outer normal vector n_i in R^2 of Gamma on a boundary
/// element \f$E_i = [z_j, z_k]\f$ reads
///
/// \f[ n_i = \frac{1}{\Vert z_k - z_j \Vert} \left( \begin{array}{c}
///              y_k - y_j \\ x_j - x_k \end{array} \right). \f]
///
///
/// If \f$\Gamma\f$ is split into Dirichlet boundary \f$\Gamma_D\f$ and Neumann
/// boundary \f$\Gamma_N\f$ which describe elements completely contained in one
/// part or the other, the mesh is represented by a submesh containing the
/// Dirichlet elements and another one containing the Neumann elements.
///
///////////////////////////////////////////////////////////////////////////////
#ifndef BOUNDARY_MESH_HPPmysolution
#define BOUNDARY_MESH_HPPmysolution

#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <istream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace mysolution {
class BoundaryMesh {
 private:
  /// type for mesh vertices (2d points)
  typedef Eigen::Matrix<double, Eigen::Dynamic, 2> coord_matrix_t;
  /// type for mesh elements (storing indices of the vertices conforming the given element)
  typedef Eigen::Matrix<int, Eigen::Dynamic, 2> elem_matrix_t;

  /// Class data
  coord_matrix_t coordinates_;
  elem_matrix_t elements_;
  bool isInitialized_;

 public:
  /**
   * Construct from matrices of coordinates and elements
   */
  BoundaryMesh(const coord_matrix_t& coords, const elem_matrix_t& elems);

  /**
   * Construct from data files
   */
  BoundaryMesh(const std::string& filename) { loadMeshFromFile(filename); }

  /**
   *  This function refines the mesh by splitting elements into half.
   */
  std::pair<BoundaryMesh, Eigen::MatrixXi> uniformRefine() const;

  /**
   *  This function returns the number of vertices on the mesh.
   */
  int numVertices() const;

  /**
   *  This function returns the number of elements on the mesh
   */
  int numElements() const;

  /**
   *  This function returns the matrix containing the mesh points.
   *
   *  @return Matrix containing the coordinates of each vertex of the 
   *  boundary mesh.
   */
  coord_matrix_t getMeshVertices() const;

  /**
   *  This function returns the matrix containing the indices of the vertices 
   *  corresponding to each element of the boundary mesh.
   */
  elem_matrix_t getMeshElements() const;

  /**
   *  This function returns the 2dvector containing the coordinates of the 
   *  vertex i.
   */
  Eigen::Vector2d getVertex(int i) const;

  /**
   *  This function returns the coordinates of the 2 vertices on the 
   *  i-th element.
   */
  std::pair<Eigen::Vector2d, Eigen::Vector2d> getElementVertices(int i) const;

  /**
   *  This function returns the index of the j-th vertex on the element i.
   */
  int getElementVertex(int i, int j) const;

  /**
   *  This function reads the .dat-files containing the mesh data and fills the
   *  matrices corresponding to elements and coordinates.
   *
   *  @param[in] filename string with the name of the mesh to be read.
   */
  void loadMeshFromFile(const std::string& filename);

  /**
   *  This function writes two .dat-files containing the mesh data: elements 
   *  and coordinates.
   *
   *  @param[in] filename string with the name of the files to be written:
   *             filename_coordinates.dat and filename_elements.dat
   */
  void writeMeshToFile(const std::string& filename);

 private:
  /**
   *  This function reads a .dat-file and uses it contents to fill the
   *  corresponding matrix data.
   *
   *  @tparam T  Type of Eigen matrix corresponding to data.
   *  @param[in] filename  String with the name of the mesh to be read.
   *  @param[out] data Matrix containing the data read from file.
   */
  template <typename T>
  void readData(const std::string& filename, T& data);

};  //end class
}  // namespace mysolution

#endif
