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

#include <CGAL/license/Polygon_mesh_processing/registration.h>

#include <CGAL/config.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#endif

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Surface_mesh_deformation.h>
#include <CGAL/Deformation_Eigen_closest_rotation_traits_3.h>
#include <CGAL/Weights/cotangent_weights.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Splitters.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Aff_transformation_3.h>

bool new_arap = true;

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {
namespace registration {

#ifdef CGAL_EIGEN3_ENABLED
using ScalarType = double;
using Vertices = Eigen::Matrix<ScalarType, Eigen::Dynamic, 3>;
using Faces = Eigen::Matrix<int, Eigen::Dynamic, 3>;
using Matrix = Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
using Vector = Eigen::VectorXd;

using SparseMat = Eigen::SparseMatrix<ScalarType>;
using SparseTriplet = Eigen::Triplet<ScalarType>;

using Index = Eigen::Index;

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

  Eigen_matrix_to_point_map a(points);

  KDTree kdtree(boost::counting_iterator<std::size_t>(0), boost::counting_iterator<std::size_t>(points.rows()), KDTree::Splitter(), Search_traits(Eigen_matrix_to_point_map(points)));
  kdtree.build();

  Eigen::MatrixXi idz(query.rows(), k);
  Eigen::MatrixXf dist(query.rows(), k);

  for (Index i = 0; i < query.rows(); ++i) {
    Point_3 query_pt = { query(i, 0), query(i, 1), query(i, 2) };
    Neighbor_search search(kdtree, query_pt, k, 0, true, Neighbor_search::Distance(Eigen_matrix_to_point_map(points)));
    Index j = 0;
    for (auto it = search.begin(); it != search.end() && j < k; ++it) {
      idz(i, j) = it->first;
      dist(i, j) = it->second;
    }
  }

  return std::make_pair(idz, dist);
}

template <typename T>
int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

void insertSparseMatrix(const SparseMat& mat, std::vector<SparseTriplet>& coefficients, size_t start_i = 0, size_t start_j = 0) {
  for (int k = 0; k < mat.outerSize(); ++k)
    for (SparseMat::InnerIterator it(mat, k); it; ++it)
      coefficients.push_back(SparseTriplet(start_i + it.row(), start_j + it.col(), it.value()));
}

template<typename T, typename Visitor, typename TriangleMesh>
void rotation(std::size_t index, Visitor &v, const TriangleMesh &source, const Vertices& X, const Vertices& Z, const std::set<int>& neighbors, const std::vector<T> &weights, std::vector<Eigen::Matrix<ScalarType, 3, 3>> &rotations) {
  using Vertex_index = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  Deformation_Eigen_closest_rotation_traits_3 cr;

  Eigen::Matrix3d cov;
  cov.setZero();

  Eigen::Vector3d x_src = X.row(index), z_src = Z.row(index);
  auto n = neighbors.begin();
  auto w = weights.begin();

  v.rotation_matrix_pre(Vertex_index(index), source);

  for (; (n != neighbors.end()) && (w != weights.end()); n++, w++) {
    Eigen::Vector3d x_dst = X.row(*n);
    Eigen::Vector3d z_dst = Z.row(*n);
    Eigen::Vector3d source_edge = x_src - x_dst;
    Eigen::Vector3d moving_edge = z_src - z_dst;

    cr.add_scalar_t_vector_t_vector_transpose(cov, *w, source_edge, moving_edge);
    v.update_covariance_matrix(cov, rotations[index]);
  }

  cr.compute_close_rotation(cov, rotations[index]);
}

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

#endif
} // namespace registration
} // namespace internal

/*!
* \ingroup PMP_registration_grp
*
* \brief computes non-rigid transformation of a mesh onto a set of oriented points.
*
* A non-rigid ICP, iterative closest point, method based on
* <A HREF="https://vgl.ict.usc.edu/Research/NonRigidRegistration/MODERN%20TECHNIQUES%20AND%20APPLICATIONS%20FOR%20REAL-TIME%20NON-RIGID%20REGISTRATION.pdf">a SIGGRAPH'16 Tutorial</A>.
* The method uses a few `correspondences` between the `source` and the `target` for the rough alignment. The iterative closest point method
* iteratively approaches the `target` by minimizing the distance between vertices of the `source` and points of the `target`.
*
* @note This function requires the \ref thirdpartyEigen library.
*
* @tparam TriangleMesh a const model of `FaceGraph`.
* @tparam PointRange is a model of `ConstRange`. The value type of
*   its iterator is the key type of the named parameter `point_map`.
* @tparam VertexTranslationMap is a property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*   as key type and a \cgal Kernel `Vector_3` as value type.
* @tparam VertexRotationMap is a property map with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*   as key type and a \cgal Kernel `Aff_transformation_3` as value type.
* @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters".
* @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param source the triangle mesh to be mapped onto `target`.
* @param target the target point set.
* @param vtm a writable vertex property map of `source` to store the translation vector of the registration.
* @param vrm a writable vertex property map of `source` to store the rotation part of the registration.
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters 1" among the ones listed below.
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters 2" providing a point_map and normal_map for the `PointRange` as listed below.
*
* \cgalNamedParamsBegin{Named Parameters 1}
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of registration iterations using ICP, iterative closest point}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`50`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_plane_weight}
*     \cgalParamDescription{the weight \f$w_2\f$ of the point to plane energy in the registration. }
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*     \cgalParamExtra{\f$w_2\f$ needs to be 0 or positive. See \ref PMPNonRigidRegistrationParameters.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_point_weight}
*     \cgalParamDescription{the weight \f$w_1\f$ of the point to matching point energy in the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*     \cgalParamExtra{\f$w_1\f$ needs to be 0 or positive. See \ref PMPNonRigidRegistrationParameters.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{as_rigid_as_possible_weight}
*     \cgalParamDescription{defines the rigidity of the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`50`}
*     \cgalParamExtra{The weight \f$w_3\f$ needs to be 0 or positive. See \ref PMPNonRigidRegistrationParameters.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{maximum_matching_distance}
*     \cgalParamDescription{the maximum distance for a vertex in `source` to match with a point in `target`. The default value 0 means no maximum matching distance.}
*     \cgalParamType{double}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{correspondences}
*     \cgalParamDescription{a range of matching vertex-point pairs between the `source` and the `target`.}
*     \cgalParamType{`ConstRange` whose value type is a pair of `boost::graph_traits<TriangleMesh>::%vertex_descriptor` and the value type of `PointRange`.}
*     \cgalParamDefault{empty}
*     \cgalParamExtra{to avoid copies, this parameter can be passed using `std::cref`.}
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
*
* \cgalNamedParamsBegin{Named Parameters 2}
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the point set `target`.}
*     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
*                    of the iterator of `PointRange` and whose value type is `geom_traits::Point_3`.}
*     \cgalParamDefault{`CGAL::Identity_property_map<geom_traits::Point_3>`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{normal_map}
*     \cgalParamDescription{a property map associating normals to the elements of the point set `target`.}
*     \cgalParamType{a model of `ReadablePropertyMap` whose key type is the value type
*                    of the iterator of `PointRange` and whose value type is `geom_traits::Vector_3`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/

template <typename TriangleMesh,
          typename PointRange,
          typename VertexTranslationMap,
          typename VertexRotationMap,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
void non_rigid_mesh_to_points_registration(const TriangleMesh& source,
  const PointRange& target,
  VertexTranslationMap& vtm,
  VertexRotationMap& vrm,
  const NamedParameters1& np1 = parameters::default_values(),
  const NamedParameters2& np2 = parameters::default_values())
{
#ifdef CGAL_EIGEN3_ENABLED
  const size_t iter = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::number_of_iterations), 50);
  const double w2 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::point_to_plane_weight), 2);
  const double w1 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::point_to_point_weight), 0.1);
  const double w3 = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::as_rigid_as_possible_weight), 20);
  const double max_matching_dist = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::maximum_matching_distance), 0);

  using namespace internal::registration;
  Vertices X(num_vertices(source), 3), Y(target.size(), 3);
  Faces XF(num_faces(source), 3);

  using Vertex_index = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  using NP_helper = Point_set_processing_3_np_helper<PointRange, NamedParameters2>;
  using Point_map = typename NP_helper::Point_map;
  using Normal_map = typename NP_helper::Normal_map;

  Point_map point_map = NP_helper::get_const_point_map(target, np2);
  Normal_map normal_map = NP_helper::get_normal_map(target, np2);

  using Gt = typename GetGeomTraits<TriangleMesh, NamedParameters1>::type;
  using Vertex_point_map = typename GetVertexPointMap<TriangleMesh, NamedParameters1>::type;
  using Point = typename Gt::Point_3;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np1, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, source));

  typedef typename CGAL::internal_np::Lookup_named_param_def<internal_np::correspondences_t, NamedParameters1, std::vector<std::pair<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, std::size_t>>>::reference Correspondences_type;
  Correspondences_type correspondences = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter_reference(np1, CGAL::internal_np::correspondences), std::vector<std::pair<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, std::size_t>>());

  using ARAP_visitor = typename CGAL::internal::Types_selectors < TriangleMesh, Vertex_point_map, SRE_ARAP>::ARAP_visitor;
  using WC = typename CGAL::internal::Types_selectors < TriangleMesh, Vertex_point_map, SRE_ARAP>::Weight_calculator;
  WC wc;

  std::size_t idx = 0;
  std::map<typename boost::graph_traits<TriangleMesh>::vertex_descriptor, std::size_t> vi;
  std::vector<typename boost::graph_traits<TriangleMesh>::vertex_descriptor> iv;
  std::vector<std::vector<double> > he_weights;

  iv.resize(num_vertices(source));
  for (auto v : vertices(source)) {
    X(idx, 0) = get(vpm, v).x();
    X(idx, 1) = get(vpm, v).y();
    X(idx, 2) = get(vpm, v).z();
    iv[idx] = v;
    vi[v] = idx;
    idx++;
  }

  idx = 0;
  for (auto f : faces(source)) {
    auto h = halfedge(f, source);
    XF(idx, 0) = CGAL::target(h, source);
    h = next(h, source);
    XF(idx, 1) = CGAL::target(h, source);
    h = next(h, source);
    XF(idx, 2) = CGAL::target(h, source);
    idx++;
  }
  std::cout << std::endl;

  Eigen::MatrixXi corr(correspondences.size(), 2);
  for (size_t i = 0; i < correspondences.size(); ++i) {
    corr.row(i) << correspondences[i].first, correspondences[i].second;
  }

  if (corr.rows() > 0)
    std::cout << "# correspondences = " << corr.rows() << std::endl;

  Vertices NY(target.size(), 3);

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
  for (Index i = 0; i < XF.rows(); ++i) {
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

  // Calculate edge weights
  he_weights.resize(X.rows());
  for (Index i = 0; i < X.rows(); ++i) {
    std::size_t vi = XF(i, 0);
    he_weights[vi].resize(neighbors[vi].size());
    std::size_t idx = 0;
    for (int vj : neighbors[vi]) {
      auto h = halfedge(Vertex_index(vi), Vertex_index(vj), source);
      assert(h.second);
      he_weights[vi][idx++] = wc(h.first, source, vpm);
    }
  }

  ARAP_visitor visitor;
  visitor.init(source, vpm);

  // Non-rigid ICP
  Vertices Z(X);
  Index dim = Z.rows() * Z.cols();

  std::vector<SparseTriplet> edge_coefficients;

  // build Ni
  Eigen::MatrixXi Ni(XF.rows() * XF.cols(), 1);
  idx = 0;
  for (Index i = 0; i < XF.rows(); ++i) {
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
  for (Index i = 0; i < XF.rows(); ++i) {
    int a = XF(i, 0);
    int b = XF(i, 1);
    int c = XF(i, 2);
    edge_coefficients.push_back(SparseTriplet(idx, b, 1));
    idx++;
    edge_coefficients.push_back(SparseTriplet(idx, c, 1));
    idx++;
    edge_coefficients.push_back(SparseTriplet(idx, a, 1));
    idx++;
  }
  SparseMat Nr(XF.rows() * XF.cols(), X.rows());
  Nr.setFromTriplets(edge_coefficients.begin(), edge_coefficients.end());

  // build MX
  edge_coefficients.clear();
  edge_coefficients.reserve(XF.rows() * XF.cols() * 2);
  idx = 0;
  for (Index i = 0; i < XF.rows(); ++i) {
    int a = XF(i, 0);
    int b = XF(i, 1);
    int c = XF(i, 2);
    if (new_arap) {
      auto ha = halfedge(Vertex_index(a), Vertex_index(b), source);
      double w_a = wc(ha.first, source, vpm);
      auto hb = halfedge(Vertex_index(b), Vertex_index(c), source);
      double w_b = wc(hb.first, source, vpm);
      auto hc = halfedge(Vertex_index(c), Vertex_index(a), source);
      double w_c = wc(hc.first, source, vpm);
      edge_coefficients.push_back(SparseTriplet(idx, a, w_a));
      edge_coefficients.push_back(SparseTriplet(idx, b, -w_a));
      idx++;
      edge_coefficients.push_back(SparseTriplet(idx, b, w_b));
      edge_coefficients.push_back(SparseTriplet(idx, c, -w_b));
      idx++;
      edge_coefficients.push_back(SparseTriplet(idx, c, w_c));
      edge_coefficients.push_back(SparseTriplet(idx, a, -w_c));
      idx++;
    }
    else {
      edge_coefficients.push_back(SparseTriplet(idx, a, 1));
      edge_coefficients.push_back(SparseTriplet(idx, b, -1));
      idx++;
      edge_coefficients.push_back(SparseTriplet(idx, b, 1));
      edge_coefficients.push_back(SparseTriplet(idx, c, -1));
      idx++;
      edge_coefficients.push_back(SparseTriplet(idx, c, 1));
      edge_coefficients.push_back(SparseTriplet(idx, a, -1));
      idx++;
    }
  }
  SparseMat B(XF.rows() * XF.cols(), X.rows());
  B.setFromTriplets(edge_coefficients.begin(), edge_coefficients.end());
  Index edim = B.rows();

  Vertices BX = B * X;// edges in order: b to a, c to b, a to c.
  Vertices BX_original(BX);

  std::vector<SparseTriplet> coefficients;
  SparseMat A(Z.rows() + 2 * dim + 3 * edim, dim + dim);
  for (Index i = 0; i < dim; ++i) {
    coefficients.push_back(SparseTriplet(Z.rows() + i, i, 1));
  }

  insertSparseMatrix(B, coefficients, Z.rows() + dim, 0);
  insertSparseMatrix(B, coefficients, Z.rows() + dim + edim, Z.rows());
  insertSparseMatrix(B, coefficients, Z.rows() + dim + 2 * edim, 2 * Z.rows());

  Matrix b(Z.rows() + 2 * dim + 3 * edim, 1);

  std::vector<SparseTriplet> weight_coefficients; // system regularizer
  SparseMat W(dim + dim, dim + dim);
  for (Index i = dim; i < dim + dim; ++i) {
    weight_coefficients.push_back(SparseTriplet(i, i, 1));
  }
  W.setFromTriplets(weight_coefficients.begin(), weight_coefficients.end());

  weight_coefficients.clear();
  SparseMat D(Z.rows() + 2 * dim + 3 * edim, Z.rows() + 2 * dim + 3 * edim);
  for (Index i = 0; i < Z.rows(); ++i) {
    weight_coefficients.push_back(SparseTriplet(i, i, w2)); // point to plane energy
  }
  for (Index i = Z.rows(); i < Z.rows() + dim; ++i) {
    weight_coefficients.push_back(SparseTriplet(i, i, w1)); // point to point energy
  }
  for (Index i = Z.rows() + dim; i < Z.rows() + dim + 3 * edim; ++i) {
    weight_coefficients.push_back(SparseTriplet(i, i, w3)); // arap energy
  }
  D.setFromTriplets(weight_coefficients.begin(), weight_coefficients.end());

  // Setting very high weight for given correspondences.
  for (Index i = 0; i < corr.rows(); ++i) {
    D.coeffRef(Z.rows() + corr(i, 0), Z.rows() + corr(i, 0)) = 1e15f;
    D.coeffRef(Z.rows() * 2 + corr(i, 0), Z.rows() * 2 + corr(i, 0)) = 1e15f;
    D.coeffRef(Z.rows() * 3 + corr(i, 0), Z.rows() * 3 + corr(i, 0)) = 1e15f;
  }

  // Solver
  Eigen::ConjugateGradient<SparseMat, Eigen::Lower | Eigen::Upper> cg;
  //cg.setMaxIterations(1000);
  //cg.setTolerance(1e-6);

#if EIGEN_VERSION_AT_LEAST(3,4,90)
  Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeFullU | Eigen::ComputeFullV> svd;
#else
  Eigen::JacobiSVD<Eigen::Matrix<ScalarType, 3, 3>> svd;
#endif

  std::vector<Eigen::Matrix<ScalarType, 3, 3>> rotations(Z.rows());

  for (std::size_t i = 0; i < std::size_t(Z.rows()); i++)
    rotations[i].setIdentity();

  size_t coefficients_size = coefficients.size();

  for (size_t it = 0; it < iter; ++it) {
    std::cout << "." << std::flush;
    // Reset coefficients (removes point pairs, e.g.)
    if (it > 0) {
      coefficients.erase(coefficients.begin() + coefficients_size, coefficients.end());
    }

    // Compute correspondence
    //Eigen::VectorXi idz = nearest_neighbor(V1, Z).first.col(0);
    std::pair<Eigen::MatrixXi, Eigen::MatrixXf> nn_result = nearest_neighbor(Y, Z);
    Eigen::VectorXi idz = nn_result.first.col(0);
    Eigen::VectorXf dist = nn_result.second.col(0);

    if (max_matching_dist > 0) {
      D.setFromTriplets(weight_coefficients.begin(), weight_coefficients.end()); // reset weights
      // prune correspondences that are too distant
      int count = 0;
      for (Index i = 0; i < Z.rows(); ++i) {
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
      for (Index i = 0; i < corr.rows(); ++i) {
        D.coeffRef(Z.rows() + corr(i, 0), Z.rows() + corr(i, 0)) = 1e15f;
        D.coeffRef(Z.rows() * 2 + corr(i, 0), Z.rows() * 2 + corr(i, 0)) = 1e15f;
        D.coeffRef(Z.rows() * 3 + corr(i, 0), Z.rows() * 3 + corr(i, 0)) = 1e15f;
      }
    }

    // Restore hard correspondences
    for (Index i = 0; i < corr.rows(); ++i) {
      idz(corr(i, 0)) = corr(i, 1);
    }

#if EIGEN_VERSION_AT_LEAST(3,4,0)
    Vertices P(Y(idz, Eigen::indexing::all)); // target points
    Vertices NP(NY(idz, Eigen::indexing::all)); // target normals
#else
    Vertices P(idz.size(), 3); // target points
    Vertices NP(idz.size(), 3); // target normals
    for (Index i = 0; i < idz.rows(); i++) {
      P.row(i) = Y.row(idz(i));
      NP.row(i) = NY.row(idz(i));
    }
#endif

    // Filling coefficients for first part of A -> point_to_plane_energy
    for (Index i = 0; i < Z.rows(); ++i) {
      coefficients.push_back(SparseTriplet(i, i, NP(i, 0)));
      coefficients.push_back(SparseTriplet(i, Z.rows() + i, NP(i, 1)));
      coefficients.push_back(SparseTriplet(i, 2 * Z.rows() + i, NP(i, 2)));
    }

    // Filling coefficients from edges for last part of A -> as_rigid_as_possible_energy
    for (Index i = 0; i < edim; ++i) {
      size_t ni = Ni(i);
      coefficients.push_back(SparseTriplet(Z.rows() + dim + i, ni + dim, BX(i, 1)));
      coefficients.push_back(SparseTriplet(Z.rows() + dim + i, Z.rows() + ni + dim, -BX(i, 2)));
      coefficients.push_back(SparseTriplet(Z.rows() + dim + edim + i, ni + dim, -BX(i, 0)));
      coefficients.push_back(SparseTriplet(Z.rows() + dim + edim + i, 2 * Z.rows() + ni + dim, BX(i, 2)));
      coefficients.push_back(SparseTriplet(Z.rows() + dim + 2 * edim + i, Z.rows() + ni + dim, BX(i, 0)));
      coefficients.push_back(SparseTriplet(Z.rows() + dim + 2 * edim + i, 2 * Z.rows() + ni + dim, -BX(i, 1)));
    }

    for (Index i = 0; i < Z.rows(); ++i) {
      b(i) = NP(i, 0) * P(i, 0) + NP(i, 1) * P(i, 1) + NP(i, 2) * P(i, 2); // distance of matched point to origin for point to plane energy
      // Coordinates of matched point for point to point energy
      b(i + Z.rows()) = P(i, 0);
      b(i + Z.rows() * 2) = P(i, 1);
      b(i + Z.rows() * 3) = P(i, 2);
    }
    for (Index i = 0; i < edim; ++i) {
      b(i + Z.rows() + dim) = BX(i, 0);
      b(i + Z.rows() + dim + edim) = BX(i, 1);
      b(i + Z.rows() + dim + edim * 2) = BX(i, 2);
    }

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    if (it == 0) {
      cg.analyzePattern(A.transpose() * D * A + W);
    }

    cg.factorize(A.transpose() * D * A + W);
    Vector x = cg.solve(A.transpose() * D * b);

    // Solution
#if EIGEN_VERSION_AT_LEAST(3,4,0)
    Z(Eigen::indexing::all, 0) = x(Eigen::seq(0, Z.rows() - 1));
    Z(Eigen::indexing::all, 1) = x(Eigen::seq(Z.rows(), 2 * Z.rows() - 1));
    Z(Eigen::indexing::all, 2) = x(Eigen::seq(2 * Z.rows(), 3 * Z.rows() - 1));
#else
    for (Index i = 0; i < Z.rows(); i++) {
      Z(i, 0) = x(i);
      Z(i, 1) = x(Z.rows() + i);
      Z(i, 2) = x(2 * Z.rows() + i);
    }
#endif

    // Update edge neighborhoods by new local rotation
    if (it == (iter - 1)) {
      if (new_arap) {
        for (Index i = 0; i < Z.rows(); ++i)
          rotation(std::size_t(i), visitor, source, X, Z, neighbors[i], he_weights[i], rotations);

        for (Index i = 0; i < edim; ++i) {
          BX.row(i) = BX_original.row(i) * rotations[Ni(i)].transpose();
        }
      }
      else {
        for (Index i = 0; i < Z.rows(); ++i) {
#if EIGEN_VERSION_AT_LEAST(3,4,0)
          std::vector<int> nbrs(neighbors[i].begin(), neighbors[i].end());
          Matrix A = X(nbrs, Eigen::indexing::all).rowwise() - X.row(i);
          Matrix B_ = Z(nbrs, Eigen::indexing::all).rowwise() - Z.row(i);
#else
          Matrix A(neighbors[i].size(), 3);
          Matrix B_(neighbors[i].size(), 3);
          auto xi = X.row(i);
          auto zi = Z.row(i);
          auto nit = neighbors[i].begin();
          for (std::size_t j = 0; j < neighbors[i].size(); j++, nit++) {
            A.row(j) = X.row(*nit) - xi;
            B_.row(j) = Z.row(*nit) - zi;
          }
#endif

#if EIGEN_VERSION_AT_LEAST(3,4,90)
          svd.compute(A.transpose() * B_);
#else
          svd.compute(A.transpose() * B_, Eigen::ComputeFullU | Eigen::ComputeFullV);
#endif

          rotations[i] = svd.matrixV() * svd.matrixU().transpose();
          if (rotations[i].determinant() < 0) {
            Eigen::Matrix<ScalarType, 3, 3> M = Eigen::Matrix3d::Identity();
            M(2, 2) = -1;
            rotations[i] = svd.matrixV() * M * svd.matrixU().transpose();
          }
        }
        for (Index i = 0; i < edim; ++i) {
          Matrix R = rotations[Ni(i)];
          BX.row(i) = BX.row(i) * R.transpose();
        }
      }
    }
    else
      if (new_arap) {
        for (Index i = 0; i < Z.rows(); ++i)
          rotation(std::size_t(i), visitor, source, X, Z, neighbors[i], he_weights[i], rotations);

        for (Index i = 0; i < edim; ++i)
          BX.row(i) = BX_original.row(i) * rotations[Ni(i)].transpose();
      }
      else {
        // Regular matrix update is happening here (taking linearized coefficients, recreating matrix and rotating edges)
        for (Index i = 0; i < edim; ++i) {
          int ni = Ni(i);
          Matrix R = rot(x(dim + 2 * Z.rows() + ni),
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
/*
    Aff_transformation_3<Gt> t1(CGAL::TRANSLATION, CGAL::ORIGIN - x);
    Aff_transformation_3<Gt> t2(CGAL::TRANSLATION, -(CGAL::ORIGIN - z));*/
    const auto& r = rotations[idx];
    Aff_transformation_3<Gt> rota(r(0, 0), r(0, 1), r(0, 2), 0,
      r(1, 0), r(1, 1), r(1, 2), 0,
      r(2, 0), r(2, 1), r(2, 2), 0,
      1);
    put(vtm, v, z - x);
    put(vrm, v, rota);
    idx++;
  }
#else
static_assert(false, "Eigen library is required for non-rigid mesh registration");
#endif
}

/*!
* \ingroup PMP_registration_grp
*
* \brief computes non-rigid transformation of a mesh onto another mesh.
*
* A non-rigid ICP, iterative closest point, method based on
* <A HREF="https://vgl.ict.usc.edu/Research/NonRigidRegistration/MODERN%20TECHNIQUES%20AND%20APPLICATIONS%20FOR%20REAL-TIME%20NON-RIGID%20REGISTRATION.pdf">a SIGGRAPH'16 Tutorial</A>.
* The method uses optional `correspondences` between the `source` and the `target` for the rough alignment. The iterative closest point method
* iteratively approaches the `target` by minimizing the distance between vertices of the `source` and vertices of the `target`.
*
* @note This function requires the \ref thirdpartyEigen library.
*
* @tparam TriangleMesh1 a const model of `FaceGraph`.
* @tparam TriangleMesh2 a const model of `FaceGraph`.
* @tparam VertexTranslationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
*   as key type and a \cgal Kernel `Vector_3` as value type.
* @tparam VertexRotationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
*   as key type and a \cgal Kernel `Aff_transformation_3` as value type.
* @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters1".
* @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters2".
*
* @param source the triangle mesh to be mapped onto `target`.
* @param target the target triangle mesh.
* @param vtm a writable vertex property map of `source` to store the translation vector of the registration.
* @param vrm a writable vertex property map of `source` to store the rotation part of the registration.
* @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters 1" of the `source` and the method among the ones listed below.
* @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters 2" of the `target` providing a vertex point map and a vertex normal map as listed below.
*
* \cgalNamedParamsBegin{Named Parameters 1}
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of registration iterations using ICP, iterative closest point}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`50`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_plane_weight}
*     \cgalParamDescription{the weight \f$w_2\f$ of the point to plane energy in the registration. }
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*     \cgalParamExtra{\f$w_2\f$ needs to be 0 or positive. See \ref PMPNonRigidRegistrationParameters.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{point_to_point_weight}
*     \cgalParamDescription{the weight \f$w_1\f$ of the point to matching point energy in the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`1`}
*     \cgalParamExtra{\f$w_1\f$ needs to be 0 or positive. See \ref PMPNonRigidRegistrationParameters.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{as_rigid_as_possible_weight}
*     \cgalParamDescription{defines the rigidity of the registration}
*     \cgalParamType{double}
*     \cgalParamDefault{`50`}
*     \cgalParamExtra{The weight \f$w_3\f$ needs to be 0 or positive. See \ref PMPNonRigidRegistrationParameters.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{maximum_matching_distance}
*     \cgalParamDescription{the maximum distance for a vertex in `source` to match with a point in `target`. The default value 0 means no maximum matching distance.}
*     \cgalParamType{double}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{correspondences}
*     \cgalParamDescription{a range of matching vertex pairs between the `source` and the `target`.}
*     \cgalParamType{`ConstRange` whose value type is a pair of `boost::graph_traits<TriangleMesh1>::%vertex_descriptor` and `boost::graph_traits<TriangleMesh2>::%vertex_descriptor`.}
*     \cgalParamDefault{empty}
*     \cgalParamExtra{to avoid copies, this parameter can be passed using `std::cref`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `source`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
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
* \cgalNamedParamsBegin{Named Parameters 2}
*   \cgalParamNBegin{vertex_normal_map}
*     \cgalParamDescription{a property map associating normals to the vertices of `target`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Vector_3` as value type}
*     \cgalParamDefault{`get(dynamic_vertex_property_t<Vector_3>(), target)`}
*     \cgalParamExtra{If this parameter is omitted, vertex normals will be computed using `compute_vertex_normals()`.}
*   \cgalParamNEnd
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `target`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh2>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`get_const_property_map(CGAL::vertex_point, target)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh2`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/

template <typename TriangleMesh1, typename TriangleMesh2,
  typename VertexTranslationMap,
  typename VertexRotationMap,
  typename NamedParameters1 = parameters::Default_named_parameters,
  typename NamedParameters2 = parameters::Default_named_parameters>
void non_rigid_mesh_to_mesh_registration(const TriangleMesh1& source,
  const TriangleMesh2& target,
  VertexTranslationMap& vtm,
  VertexRotationMap& vrm,
  const NamedParameters1& np1 = parameters::default_values(),
  const NamedParameters2& np2 = parameters::default_values())
{
#ifdef CGAL_EIGEN3_ENABLED
  using Gt2 = typename GetGeomTraits<TriangleMesh2, NamedParameters2>::type;
  using Vertex_point_map = typename GetVertexPointMap<TriangleMesh2, NamedParameters2>::type;
  using Vector_map_tag = dynamic_vertex_property_t<typename Gt2::Vector_3>;
  using Default_vector_map = typename boost::property_map<TriangleMesh2, Vector_map_tag>::const_type;
  using Vertex_normal_map = typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
    NamedParameters2,
    Default_vector_map>::type;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np2, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, target));
  Vertex_normal_map vnm = parameters::choose_parameter(parameters::get_parameter(np2, internal_np::vertex_normal_map), get(Vector_map_tag(), target));

  // if the normal map is not provided, compute it
  if (parameters::is_default_parameter<NamedParameters2, internal_np::vertex_normal_map_t>::value)
    compute_vertex_normals(target, vnm, np2);

  typedef typename CGAL::internal_np::Lookup_named_param_def<internal_np::correspondences_t, NamedParameters1, std::vector<std::pair<typename boost::graph_traits<TriangleMesh1>::vertex_descriptor, typename boost::graph_traits<TriangleMesh2>::vertex_descriptor>>>::reference Correspondences_type;
  Correspondences_type correspondences = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter_reference(np1, CGAL::internal_np::correspondences), std::vector<std::pair<typename boost::graph_traits<TriangleMesh1>::vertex_descriptor, typename boost::graph_traits<TriangleMesh2>::vertex_descriptor>>());

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

  non_rigid_mesh_to_points_registration(source, points, vtm, vrm, np1.combine(parameters::correspondences(correspondences_pts)), parameters::point_map(Point_map()).normal_map(Normal_map()));

#else
  static_assert(false, "Eigen library is required for non-rigid mesh registration");
#endif
}

/*!
* \ingroup PMP_registration_grp
*
* \brief applies a non-rigid transformation to the vertices of the mesh. Face and vertex normal vectors are not updated after transformation.
*
* @tparam TriangleMesh a model of `FaceGraph`.
* @tparam VertexTranslationMap is a property map with `boost::graph_traits<TriangleMesh1>::%vertex_descriptor`
 *   as key type and a \cgal Kernel `Vector_3` as value type.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters".
*
* @param mesh the triangle mesh to be transformed.
* @param vtm a readable vertex property map of `mesh` to store the translation vector of the registration.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.
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
  typename NamedParameters = parameters::Default_named_parameters>
void apply_non_rigid_transformation(const TriangleMesh& mesh,
                                    const VertexTranslationMap& vtm,
                                    const NamedParameters& np = parameters::default_values()) {
  using Gt = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using Vertex_point_map = typename GetVertexPointMap<TriangleMesh, NamedParameters>::type;

  using Point = typename Gt::Point_3;

  Vertex_point_map vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point), get_const_property_map(CGAL::vertex_point, mesh));

  for (auto v : vertices(mesh)) {
    Point p = get(vpm, v);
    p += get(vtm, v);
    put(vpm, v, p);
  }
}
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif
