// Copyright (c) 2024 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Roberto Dyke, Sven Oesau
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_REGISTER_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_REGISTER_MESH_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>

#include <CGAL/config.h>

//#ifdef CGAL_EIGEN3_ENABLED

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Sparse>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Aff_transformation_3.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

typedef double ScalarType;
typedef Eigen::Vector<ScalarType, 3> Vertex;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, 3> Vertices;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3> Faces;
typedef Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Matrix;
typedef Eigen::Vector<ScalarType, Eigen::Dynamic> Vector;

typedef Eigen::SparseMatrix<ScalarType> SparseMat;
typedef Eigen::Triplet<ScalarType> SparseTriplet;

typedef Eigen::Index Index;

using Point_3 = Simple_cartesian<double>::Point_3;
class Eigen_matrix_to_point_map {
  const Vertices& points;
public:
  typedef Point_3 value_type;
  typedef const value_type reference;
  typedef std::size_t key_type;
  typedef boost::lvalue_property_map_tag category;
  Eigen_matrix_to_point_map(const Vertices& points) :points(points) {}
  reference operator[](key_type k) const { return Point_3(points(k, 0), points(k, 1), points(k, 2)); }
};

Eigen_matrix_to_point_map::reference get(const Eigen_matrix_to_point_map& ppmap, Eigen_matrix_to_point_map::key_type i) {
  return ppmap[i];
}

template<typename Matrix>
void dump(const Matrix& m, const std::string &filename) {
  std::ofstream file(filename);

  file << m.rows() << " " << m.cols() << std::endl;

  for (std::size_t r = 0; r < m.rows(); r++) {
    file << m(r, 0);
    for (std::size_t c = 1; c < m.cols(); c++)
      file << " " << m(r, c);
    file << std::endl;
  }

  file.close();
}

std::pair<Eigen::MatrixXi, Eigen::MatrixXf> nearest_neighbor(Vertices& points, Vertices& query, const Index k = 1) {
  using Search_traits = CGAL::Search_traits_adapter<std::size_t, Eigen_matrix_to_point_map, CGAL::Search_traits_3<Simple_cartesian<double>>>;
  using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Search_traits>;
  using KDTree = typename Neighbor_search::Tree;
  const int leaf_max_size = 10;
  const int ndim = 3;

  Eigen_matrix_to_point_map a(points);

  KDTree kdtree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.rows()), KDTree::Splitter(), Search_traits(Eigen_matrix_to_point_map(points)));
  kdtree.build();

  Eigen::MatrixXi idz(query.rows(), k);
  Eigen::MatrixXf dist(query.rows(), k);

  for (internal::Index i = 0; i < query.rows(); ++i) {
    Point_3 query_pt = { query(i, 0), query(i, 1), query(i, 2) };
    Neighbor_search search(kdtree, query_pt, k, 0, true, Neighbor_search::Distance(Eigen_matrix_to_point_map(points)));
    std::size_t j = 0;
    for (auto it = search.begin(); it != search.end() && j < k; ++it) {
      idz(i, j) = it->first;
      dist(i, j) = it->second;
    }
  }

  return std::make_pair(idz, dist);
}/*

Vertices calc_normals(Vertices& points, Faces& faces) {
  Vertices face_normals(faces.rows(), 3);
  Vertices normals = Vertices::Zero(points.rows(), 3);

  for (Index i = 0; i < faces.rows(); ++i) {
    Vertex v0 = points.row(faces(i, 0));
    Vertex v1 = points.row(faces(i, 1));
    Vertex v2 = points.row(faces(i, 2));

    Vertex n0 = (v1 - v0).cross(v2 - v1);
    face_normals.row(i) = n0.normalized();
    normals.row(faces(i, 0)) += n0; // unnormalized would respect facet areas
    normals.row(faces(i, 1)) += n0;
    normals.row(faces(i, 2)) += n0;
  }

  for (Index i = 0; i < points.rows(); ++i) {
    normals.row(i) = normals.row(i).normalized();
  }

  return normals;
}*/

template <typename T>
int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

void insertSparseMatrix(const SparseMat& mat, std::vector<SparseTriplet>& coefficients, size_t start_i = 0, size_t start_j = 0) {
  for (int k = 0; k < mat.outerSize(); ++k)
    for (SparseMat::InnerIterator it(mat, k); it; ++it)
      coefficients.push_back(SparseTriplet(start_i + it.row(), start_j + it.col(), it.value()));
}

/*
std::pair<Matrix, Vertex> point2point(Vertices X, Vertices Y) { // Transform step for rigid ICP
  assert(X.rows() == Y.rows());

  Eigen::Vector<ScalarType, 3> x_mean = X.colwise().mean();
  X.rowwise() -= x_mean.transpose();

  Eigen::Vector<ScalarType, 3> y_mean = Y.colwise().mean();
  Y.rowwise() -= y_mean.transpose();

  Matrix C = X.transpose() * Y;
  Eigen::JacobiSVD<Matrix> svd(C, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Eigen::Vector<ScalarType, 3> sgnVec;
  sgnVec << 1, 1, sign(svd.matrixU().determinant() * svd.matrixV().determinant());
  Matrix R = svd.matrixV() * sgnVec.asDiagonal() * svd.matrixU().transpose();
  Vertex t = y_mean - R * x_mean;

  return std::make_pair(R, t);
}*/

template <typename T>
Eigen::Matrix<T, 3, 3> rot(T a, T b, T c) {
  T ca = cos(a), cb = cos(b), cc = cos(c);
  T sa = sin(a), sb = sin(b), sc = sin(c);
  Eigen::Matrix<T, 3, 3> R;
  R << cb * cc, cc* sa* sb - ca * sc, ca* cc* sb + sa * sc,
    cb* sc, ca* cc + sa * sb * sc, ca* sb* sc - cc * sa,
    -sb, cb* sa, ca* cb;
  return R;
}

} // namespace internal

/*!
* \ingroup PMP_registration_grp
*
* \brief calculates non-rigid transformation of a mesh onto a set of oriented points.
*
* A non-rigid ICP, iterative closest point, method based on
* <A HREF="https://vgl.ict.usc.edu/Research/NonRigidRegistration/MODERN%20TECHNIQUES%20AND%20APPLICATIONS%20FOR%20REAL-TIME%20NON-RIGID%20REGISTRATION.pdf">a SIGGRAPH'16 Tutorial</A>.
* The method uses a few correspondences between the source and the target for the rough alignment. The iterative closest point method
* iteratively approaches the target by minimizing the distance between vertices of the source and points of the target.
*
* @tparam TriangleMesh a model of `MutableFaceGraph`.
* @tparam PointRange a model of the `concept ForwardRange` whose value type is the point type.
* @tparam VertexTranslationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
*   as key type and a \cgal Kernel `Vector_3` as value type.
* @tparam VertexRotationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
*   as key type and a \cgal Kernel `Aff_transformation_3` as value type.
* @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param source the triangle mesh that should be mapped onto target.
* @param target the target triangle mesh.
* @param vtm a writeable vertex property map of source to store the translation vector of the registration.
* @param vrm a writeable vertex property map of source to store the rotation part of the registration.
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" providing a point_map and normal_map for the PointRange
* @param correspondences a vector given matching points between the source and the target
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of registration iterations using ICP, iterative closest point}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`50`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_plane_energy}
*     \cgalParamDescription{the weight of the point to plane distance in the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_point_energy}
*     \cgalParamDescription{the weight of the point to matching point distance in the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{as_rigid_as_possible_energy}
*     \cgalParamDescription{uses the topology of the mesh to determine how a vertex it deforming with respect to its neighbors.}
*     \cgalParamType{double}
*     \cgalParamDefault{`20`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{max_matching_dist}
*     \cgalParamDescription{the maximum distance for a vertex in source to match with a point in target}
*     \cgalParamType{double}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `source`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`get_const_property_map(CGAL::vertex_point, source)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/

template <typename TriangleMesh,
          typename PointRange,
          typename VertexTranslationMap,
          typename VertexRotationMap,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
void non_rigid_mesh_to_points_registration(TriangleMesh& source,
  const PointRange& target,
  VertexTranslationMap& vtm,
  VertexRotationMap& vrm,
  const std::vector<std::pair<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, std::size_t>>& correspondences = std::vector<std::pair<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, std::size_t>>(),
  const NamedParameters1& np1 = parameters::default_values(),
  const NamedParameters2& np2 = parameters::default_values())
{
  const size_t iter = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::number_of_iterations), 50);
  const double w2 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::point_to_plane_energy), 2);
  const double w1 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::point_to_point_energy), 0.1);
  const double w3 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::as_rigid_as_possible_energy), 20);
  const double max_matching_dist = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::max_matching_dist), 0);

  internal::Vertices X(num_vertices(source), 3), Y(target.size(), 3);
  internal::Faces XF(num_faces(source), 3);

  using NP_helper = Point_set_processing_3_np_helper<PointRange, NamedParameters2>;
  using Point_map = typename NP_helper::Point_map;
  using Normal_map = typename NP_helper::Normal_map;

  Point_map point_map = NP_helper::get_const_point_map(target, np2);
  Normal_map normal_map = NP_helper::get_normal_map(target, np2);

  using Gt = typename GetGeomTraits<TriangleMesh, NamedParameters1>::type;
  using Vertex_point_map = typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type;
  using Point = typename Gt::Point_3;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, source));

  std::size_t idx = 0;
  for (auto v : vertices(source)) {
    X(idx, 0) = get(vpm, v).x();
    X(idx, 1) = get(vpm, v).y();
    X(idx, 2) = get(vpm, v).z();
    idx++;
  }

  idx = 0;
  for (auto f : faces(source)) {
    Vertex_around_face_circulator vit = Vertex_around_face_circulator(halfedge(f, source), source);
    XF(idx, 0) = *vit++;
    XF(idx, 1) = *vit++;
    XF(idx, 2) = *vit++;
    idx++;
  }
  std::cout << std::endl;

  Eigen::MatrixXi corr(correspondences.size(), 2);
  for (size_t i = 0; i < correspondences.size(); ++i) {
    corr.row(i) << correspondences[i].first, correspondences[i].second;
  }

  if (corr.rows() > 0)
    std::cout << "# correspondences = " << corr.rows() << std::endl;

  internal::Vertices NY(target.size(), 3);

  idx = 0;
  for (auto p : target) {
    Y(idx, 0) = get(point_map, p).x();
    Y(idx, 1) = get(point_map, p).y();
    Y(idx, 2) = get(point_map, p).z();
    NY(idx, 0) = get(normal_map, p).x();
    NY(idx, 1) = get(normal_map, p).y();
    NY(idx, 2) = get(normal_map, p).z();
    idx++;
  }
  std::cout << std::endl;

  std::vector<std::set<int>> neighbors(X.rows());
  for (internal::Index i = 0; i < XF.rows(); ++i) {
    int a = XF(i, 0);
    int b = XF(i, 1);
    int c = XF(i, 2);
    neighbors[a].insert(b);
    neighbors[a].insert(c);
    neighbors[b].insert(c);
    neighbors[b].insert(a);
    neighbors[c].insert(a);
    neighbors[c].insert(b);
  }

  // Non-rigid ICP
  internal::Vertices Z(X);
  internal::Index dim = Z.rows() * Z.cols();

  std::vector<internal::SparseTriplet> edge_coefficients;

  // build Ni
  Eigen::MatrixXi Ni(XF.rows() * XF.cols(), 1);
  idx = 0;
  for (internal::Index i = 0; i < XF.rows(); ++i) {
    int a = XF(i, 0);
    int b = XF(i, 1);
    int c = XF(i, 2);
    Ni(idx++) = b;
    Ni(idx++) = c;
    Ni(idx++) = a;
  }

  // build Nr
  edge_coefficients.clear();
  edge_coefficients.reserve(XF.rows() * XF.cols());
  idx = 0;
  for (internal::Index i = 0; i < XF.rows(); ++i) {
    int a = XF(i, 0);
    int b = XF(i, 1);
    int c = XF(i, 2);
    edge_coefficients.push_back(internal::SparseTriplet(idx, b, 1));
    idx++;
    edge_coefficients.push_back(internal::SparseTriplet(idx, c, 1));
    idx++;
    edge_coefficients.push_back(internal::SparseTriplet(idx, a, 1));
    idx++;
  }
  internal::SparseMat Nr(XF.rows() * XF.cols(), X.rows());
  Nr.setFromTriplets(edge_coefficients.begin(), edge_coefficients.end());

  // build MX
  edge_coefficients.clear();
  edge_coefficients.reserve(XF.rows() * XF.cols() * 2);
  idx = 0;
  for (internal::Index i = 0; i < XF.rows(); ++i) {
    int a = XF(i, 0);
    int b = XF(i, 1);
    int c = XF(i, 2);
    edge_coefficients.push_back(internal::SparseTriplet(idx, a, 1));
    edge_coefficients.push_back(internal::SparseTriplet(idx, b, -1));
    idx++;
    edge_coefficients.push_back(internal::SparseTriplet(idx, b, 1));
    edge_coefficients.push_back(internal::SparseTriplet(idx, c, -1));
    idx++;
    edge_coefficients.push_back(internal::SparseTriplet(idx, c, 1));
    edge_coefficients.push_back(internal::SparseTriplet(idx, a, -1));
    idx++;
  }
  internal::SparseMat B(XF.rows() * XF.cols(), X.rows());
  B.setFromTriplets(edge_coefficients.begin(), edge_coefficients.end());
  internal::Index edim = B.rows();

  internal::Vertices BX = B * X;
  internal::Vertices BX_original(BX);

  std::vector<internal::SparseTriplet> coefficients;
  internal::SparseMat A(Z.rows() + 2 * dim + 3 * edim, dim + dim);
  for (internal::Index i = 0; i < dim; ++i) {
    coefficients.push_back(internal::SparseTriplet(Z.rows() + i, i, 1));
  }

  internal::insertSparseMatrix(B, coefficients, Z.rows() + dim, 0);
  internal::insertSparseMatrix(B, coefficients, Z.rows() + dim + edim, Z.rows());
  internal::insertSparseMatrix(B, coefficients, Z.rows() + dim + 2 * edim, 2 * Z.rows());

  internal::Matrix b(Z.rows() + 2 * dim + 3 * edim, 1);

  std::vector<internal::SparseTriplet> weight_coefficients; // system regularizer
  internal::SparseMat W(dim + dim, dim + dim);
  for (internal::Index i = dim; i < dim + dim; ++i) {
    weight_coefficients.push_back(internal::SparseTriplet(i, i, 1));
  }
  W.setFromTriplets(weight_coefficients.begin(), weight_coefficients.end());

  weight_coefficients.clear();
  internal::SparseMat D(Z.rows() + 2 * dim + 3 * edim, Z.rows() + 2 * dim + 3 * edim);
  for (internal::Index i = 0; i < Z.rows(); ++i) {
    weight_coefficients.push_back(internal::SparseTriplet(i, i, w2));
  }
  for (internal::Index i = Z.rows(); i < Z.rows() + dim; ++i) {
    weight_coefficients.push_back(internal::SparseTriplet(i, i, w1));
  }
  for (internal::Index i = Z.rows() + dim; i < Z.rows() + dim + 3 * edim; ++i) {
    weight_coefficients.push_back(internal::SparseTriplet(i, i, w3));
  }
  D.setFromTriplets(weight_coefficients.begin(), weight_coefficients.end());

  // Setting very high weight for given correspondences.
  for (internal::Index i = 0; i < corr.rows(); ++i) {
    D.coeffRef(Z.rows() + corr(i, 0), Z.rows() + corr(i, 0)) = 1e15f;
    D.coeffRef(Z.rows() * 2 + corr(i, 0), Z.rows() * 2 + corr(i, 0)) = 1e15f;
    D.coeffRef(Z.rows() * 3 + corr(i, 0), Z.rows() * 3 + corr(i, 0)) = 1e15f;
  }

  // Solver
  Eigen::ConjugateGradient<internal::SparseMat, Eigen::Lower | Eigen::Upper> cg;
  //cg.setMaxIterations(1000);
  //cg.setTolerance(1e-6);
  Eigen::JacobiSVD<Eigen::Matrix<internal::ScalarType, 3, 3>> svd;
  std::vector<Eigen::Matrix<internal::ScalarType, 3, 3>> Rotations(Z.rows());

  size_t coefficients_size = coefficients.size();

  for (size_t it = 0; it < iter; ++it) {
    if (it > 0) {
      coefficients.erase(coefficients.begin() + coefficients_size, coefficients.end());
    }

    // Compute correspondence
    //Eigen::VectorXi idz = nearest_neighbor(V1, Z).first.col(0);
    std::pair<Eigen::MatrixXi, Eigen::MatrixXf> nn_result = internal::nearest_neighbor(Y, Z);
    Eigen::VectorXi idz = nn_result.first.col(0);
    Eigen::VectorXf dist = nn_result.second.col(0);

    if (max_matching_dist > 0) {
      D.setFromTriplets(weight_coefficients.begin(), weight_coefficients.end()); // reset weights
      // prune correspondences that are too distant
      int count = 0;
      for (internal::Index i = 0; i < Z.rows(); ++i) {
        if (dist[i] > max_matching_dist) {
          count++;
          D.coeffRef(i, i) = 0;
          D.coeffRef(Z.rows() + i, Z.rows() + i) = 0;
          D.coeffRef(Z.rows() * 2 + i, Z.rows() * 2 + i) = 0;
          D.coeffRef(Z.rows() * 3 + i, Z.rows() * 3 + i) = 0;
        }
      }
      std::cout << "active = " << Z.rows() - count << " / " << Z.rows() << std::endl;
      // ensure hard correspondences have not been removed
      for (internal::Index i = 0; i < corr.rows(); ++i) {
        D.coeffRef(Z.rows() + corr(i, 0), Z.rows() + corr(i, 0)) = 1e15f;
        D.coeffRef(Z.rows() * 2 + corr(i, 0), Z.rows() * 2 + corr(i, 0)) = 1e15f;
        D.coeffRef(Z.rows() * 3 + corr(i, 0), Z.rows() * 3 + corr(i, 0)) = 1e15f;
      }
    }

    for (internal::Index i = 0; i < corr.rows(); ++i) {
      idz(corr(i, 0)) = corr(i, 1);
    }

    internal::Vertices P(Y(idz, Eigen::indexing::all)); // target points
    internal::Vertices NP(NY(idz, Eigen::indexing::all)); // target normals

    for (internal::Index i = 0; i < Z.rows(); ++i) {
      coefficients.push_back(internal::SparseTriplet(i, i, NP(i, 0)));
      coefficients.push_back(internal::SparseTriplet(i, Z.rows() + i, NP(i, 1)));
      coefficients.push_back(internal::SparseTriplet(i, 2 * Z.rows() + i, NP(i, 2)));
    }

    for (internal::Index i = 0; i < edim; ++i) {
      size_t ni = Ni(i);
      coefficients.push_back(internal::SparseTriplet(Z.rows() + dim + i, ni + dim, BX(i, 1)));
      coefficients.push_back(internal::SparseTriplet(Z.rows() + dim + i, Z.rows() + ni + dim, -BX(i, 2)));
      coefficients.push_back(internal::SparseTriplet(Z.rows() + dim + edim + i, ni + dim, -BX(i, 0)));
      coefficients.push_back(internal::SparseTriplet(Z.rows() + dim + edim + i, 2 * Z.rows() + ni + dim, BX(i, 2)));
      coefficients.push_back(internal::SparseTriplet(Z.rows() + dim + 2 * edim + i, Z.rows() + ni + dim, BX(i, 0)));
      coefficients.push_back(internal::SparseTriplet(Z.rows() + dim + 2 * edim + i, 2 * Z.rows() + ni + dim, -BX(i, 1)));
    }

    for (internal::Index i = 0; i < Z.rows(); ++i) {
      b(i) = NP(i, 0) * P(i, 0) + NP(i, 1) * P(i, 1) + NP(i, 2) * P(i, 2);
      b(i + Z.rows()) = P(i, 0);
      b(i + Z.rows() * 2) = P(i, 1);
      b(i + Z.rows() * 3) = P(i, 2);
    }
    for (internal::Index i = 0; i < edim; ++i) {
      b(i + Z.rows() + dim) = BX(i, 0);
      b(i + Z.rows() + dim + edim) = BX(i, 1);
      b(i + Z.rows() + dim + edim * 2) = BX(i, 2);
    }

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // There does not seem to be an advantage in performance
    /*if (it == 0) {
      cg.analyzePattern(A.transpose() * D * A + W);
    }*/
    cg.factorize(A.transpose() * D * A + W);
    internal::Vector x = cg.solve(A.transpose() * D * b);

    // Solution
    Z(Eigen::indexing::all, 0) = x(Eigen::seq(0, Z.rows() - 1));
    Z(Eigen::indexing::all, 1) = x(Eigen::seq(Z.rows(), 2 * Z.rows() - 1));
    Z(Eigen::indexing::all, 2) = x(Eigen::seq(2 * Z.rows(), 3 * Z.rows() - 1));


    // Update edge neighborhoods by new local rotation
    if (it == (iter - 1)) {
      // Should replace with CGAL's ARAP implementation https://github.com/CGAL/cgal/blob/master/Surface_mesh_deformation/include/CGAL/Surface_mesh_deformation.h
      // See also: compute_close_rotation https://github.com/CGAL/cgal/blob/master/Surface_mesh_deformation/include/CGAL/Deformation_Eigen_closest_rotation_traits_3.h
      for (internal::Index i = 0; i < Z.rows(); ++i) {
        std::vector<int> nbrs(neighbors[i].begin(), neighbors[i].end());
        internal::Matrix A = X(nbrs, Eigen::indexing::all).rowwise() - X.row(i);
        internal::Matrix B_ = Z(nbrs, Eigen::indexing::all).rowwise() - Z.row(i);

        svd.compute(A.transpose() * B_, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Rotations[i] = svd.matrixV() * svd.matrixU().transpose();
        if (Rotations[i].determinant() < 0) {
          Eigen::Matrix<internal::ScalarType, 3, 3> M = Eigen::Matrix3d::Identity();
          M(2, 2) = -1;
          Rotations[i] = svd.matrixV() * M * svd.matrixU().transpose();
        }
      }
      for (internal::Index i = 0; i < edim; ++i) {
        internal::Matrix R = Rotations[Ni(i)];
        BX.row(i) = BX.row(i) * R.transpose();
      }
    }
    else {
      // Regular matrix update is happening here (taking linearized coefficients, recreating matrix and updating rotations)
      for (internal::Index i = 0; i < edim; ++i) {
        int ni = Ni(i);
        internal::Matrix R = internal::rot(x(dim + 2 * Z.rows() + ni),
          x(dim + Z.rows() + ni),
          x(dim + ni));
        BX.row(i) = BX.row(i) * R.transpose();
      }
    }
  }

  idx = 0;
  for (auto v : vertices(source)) {
    Point z(Z(idx, 0), Z(idx, 1), Z(idx, 2));
    Point x(X(idx, 0), X(idx, 1), X(idx, 2));
    Aff_transformation_3<Gt> t1(CGAL::TRANSLATION, CGAL::ORIGIN - x);
    Aff_transformation_3<Gt> t2(CGAL::TRANSLATION, -(CGAL::ORIGIN - z));
    const auto& r = Rotations[idx];
    Aff_transformation_3<Gt> rota(r(0, 0), r(0, 1), r(0, 2), 0,
      r(1, 0), r(1, 1), r(1, 2), 0,
      r(2, 0), r(2, 1), r(2, 2), 0,
      1);
    put(vtm, v, z - x);
    put(vrm, v, rota);
    idx++;
  }
}

/*!
* \ingroup PMP_registration_grp
*
* \brief calculates non-rigid transformation of a mesh onto another mesh.
*
* A non-rigid ICP, iterative closest point, method based on
* <A HREF="https://vgl.ict.usc.edu/Research/NonRigidRegistration/MODERN%20TECHNIQUES%20AND%20APPLICATIONS%20FOR%20REAL-TIME%20NON-RIGID%20REGISTRATION.pdf">a SIGGRAPH'16 Tutorial</A>.
* The method uses optional correspondences between the source and the target for the rough alignment. The iterative closest point method
* iteratively approaches the target by minimizing the distance between vertices of the source and points of the target.
*
* @tparam TriangleMesh1 a model of `MutableFaceGraph`.
* @tparam TriangleMesh2 a const model of the `MutableFaceGraph`.
* @tparam VertexTranslationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Vector_3` as value type.
* @tparam VertexRotationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Aff_transformation_3` as value type.
* @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters1"
* @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters2"
*
* @param source the triangle mesh that should be mapped onto target.
* @param target the target point range with oriented normals.
* @param vtm a writeable vertex property map of source to store the translation vector of the registration.
* @param vrm a writeable vertex property map of source to store the rotation part of the registration.
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters1" of the source and the method among the ones listed below
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters2" of the target providing a vertex point map and a vertex normal map.
* @param correspondences a vector given matching points between the source and the target
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of registration iterations using ICP, iterative closest point}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`50`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_plane_energy}
*     \cgalParamDescription{the weight of the point to plane distance in the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_point_energy}
*     \cgalParamDescription{the weight of the point to matching point distance in the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{as_rigid_as_possible_energy}
*     \cgalParamDescription{uses the topology of the mesh to determine how a vertex it deforming with respect to its neighbors.}
*     \cgalParamType{double}
*     \cgalParamDefault{`50`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{max_matching_dist}
*     \cgalParamDescription{the maximum distance for a vertex in source to match with a point in target. 0 means that there is no max distance.}
*     \cgalParamType{double}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `source`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`get_const_property_map(CGAL::vertex_point, source)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh1`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/

template <typename TriangleMesh1, typename TriangleMesh2,
  typename VertexTranslationMap,
  typename VertexRotationMap,
  typename NamedParameters1 = parameters::Default_named_parameters,
  typename NamedParameters2 = parameters::Default_named_parameters>
void non_rigid_mesh_to_mesh_registration(TriangleMesh1& source,
  const TriangleMesh2& target,
  VertexTranslationMap& vtm,
  VertexRotationMap& vrm,
  const std::vector<std::pair<typename boost::graph_traits<TriangleMesh1>::vertex_descriptor, typename boost::graph_traits<TriangleMesh2>::vertex_descriptor>>& correspondences = std::vector<std::pair<typename boost::graph_traits<TriangleMesh1>::vertex_descriptor, typename boost::graph_traits<TriangleMesh2>::vertex_descriptor>>(),
  const NamedParameters1& np1 = parameters::default_values(),
  const NamedParameters2& np2 = parameters::default_values())
{
  using Gt2 = typename GetGeomTraits<TriangleMesh2, NamedParameters2>::type;
  //using Point = typename GeomTraits2::Point_3;
  using Vertex_point_map = typename GetVertexPointMap<TriangleMesh2, NamedParameters2>::type;
  using Vector_map_tag = dynamic_vertex_property_t<typename Gt2::Vector_3>;
  using Default_vector_map = typename boost::property_map<TriangleMesh2, Vector_map_tag>::const_type;
  using Vertex_normal_map = typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
    NamedParameters2,
    Default_vector_map>::type;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np2, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, target));
  Vertex_normal_map vnm = parameters::choose_parameter(parameters::get_parameter(np2, internal_np::vertex_normal_map), get(Vector_map_tag(), target));

//   if constexpr (!parameters::is_default_parameter<NamedParameters2, internal_np::vertex_normal_map_t>::value)
//     vnm = get(Vector_map_tag(), target).first;
//   else
//     vnm = parameters::get_parameter(np2, internal_np::vertex_normal_map);

  // if the normal map is not provided, compute it
   if (parameters::is_default_parameter<NamedParameters2, internal_np::vertex_normal_map_t>::value)
     compute_vertex_normals(target, vnm, np2);

  using Gt1 = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;

  typedef std::pair<typename Gt1::Point_3, typename Gt1::Vector_3> Point_with_normal;
  typedef std::vector<Point_with_normal> Pwn_vector;
  typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
  typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
  Pwn_vector points;
  points.reserve(target.num_vertices());

  for (auto v : target.vertices())
    points.push_back(std::make_pair(get(vpm, v), get(vnm, v)));

  std::vector<std::pair<typename boost::graph_traits<TriangleMesh1>::vertex_descriptor, std::size_t>> correspondences_pts;
  correspondences_pts.reserve(correspondences.size());
  for (auto p : correspondences)
    correspondences_pts.push_back(std::make_pair(p.first, static_cast<std::size_t>(p.second)));

  non_rigid_mesh_to_points_registration(source, points, vtm, vrm, correspondences_pts, np1, parameters::point_map(Point_map()).normal_map(Normal_map()));
}

/*!
* \ingroup PMP_registration_grp
*
* \brief applies a non-rigid transformation to the vertices of the mesh. Face and vertex normal vectors are invalid after transformation.
*
* @tparam TriangleMesh a model of `MutableFaceGraph`.
* @tparam VertexTranslationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Vector_3` as value type.
* @tparam VertexRotationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Aff_transformation_3` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param mesh the triangle mesh that should be transformed.
* @param vtm a readable vertex property map of source to store the translation vector of the registration.
* @param vrm a readable vertex property map of source to store the rotation part of the registration.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `source`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`get_const_property_map(CGAL::vertex_point, source)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*/
template <typename TriangleMesh,
  typename VertexTranslationMap,
  typename VertexRotationMap,
  typename NamedParameters = parameters::Default_named_parameters>
void apply_non_rigid_transformation(TriangleMesh& mesh,
                                    const VertexTranslationMap& vtm,
                                    const VertexRotationMap& vrm,
                                    const NamedParameters& np = parameters::default_values()) {
  using Gt =  typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using Vertex_point_map = typename GetVertexPointMap<TriangleMesh, NamedParameters>::type;
  using Vector_map_tag = dynamic_vertex_property_t<typename Gt::Vector_3>;
  using Vector_map_tag = dynamic_vertex_property_t<typename Gt::Vector_3>;
  using Default_vector_map = typename boost::property_map<TriangleMesh, Vector_map_tag>::type;
  using Vertex_normal_map = typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
    NamedParameters,
    Default_vector_map>::type;
  using Point = typename Gt::Point_3;
  using Vector = typename Gt::Vector_3;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point), get_property_map(CGAL::vertex_point, mesh));
  Vertex_normal_map vnm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_normal_map), get(Vector_map_tag(), mesh));


  for (auto v : vertices(mesh)) {
    Point p = get(vpm, v);
    p += get(vtm, v);
    put(vpm, v, p);
/*
    if (!parameters::is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value) {
      Vector n = get(vnm, v);
      auto rotation = get(vrm, v);
      put(vnm, v, rotation.inverse().transform(n));
    }*/
  }
}
} // namespace Polygon_mesh_processing
} // namespace CGAL

//#endif // CGAL_EIGEN3_ENABLED

#endif