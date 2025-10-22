// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Qijia Huang

#ifndef CGAL_FAST_WINDING_NUMBER_H
#define CGAL_FAST_WINDING_NUMBER_H

#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/license/Variational_medial_axis.h>
namespace CGAL {

#ifndef DOXYGEN_RUNNING
template <typename GT>
struct Tensor3
{
  using FT = typename GT::FT;
  using Vector_3 = typename GT::Vector_3;
  using Vec3 = Eigen::Matrix<FT, 3, 1>;

  using Mat3 = Eigen::Matrix<FT, 3, 3>;
  std::array<Mat3, 3> M;

  Tensor3() {
    for(auto& A : M)
      A.setZero();
  }

  inline void add(const Vector_3& v1, const Vector_3& v2, const Vector_3& v3, const Vector_3& normal, const FT& area) {
    Vec3 n(normal.x(), normal.y(), normal.z());
    Vec3 v1_eigen(v1.x(), v1.y(), v1.z());
    Vec3 v2_eigen(v2.x(), v2.y(), v2.z());
    Vec3 v3_eigen(v3.x(), v3.y(), v3.z());

    Mat3 C = (1.0 / 3.0) *
             (v1_eigen * v1_eigen.transpose() + v2_eigen * v2_eigen.transpose() + v3_eigen * v3_eigen.transpose());
    M[0].noalias() += (area * n[0]) * C;
    M[1].noalias() += (area * n[1]) * C;
    M[2].noalias() += (area * n[2]) * C;
  }

  Tensor3& operator+=(const Tensor3& o) {
    for(int k = 0; k < 3; ++k)
      M[k] += o.M[k];
    return *this;
  }

  FT operator*(const Tensor3& T) const {
    FT sum = 0;
    for(int k = 0; k < 3; ++k)
      sum += M[k].cwiseProduct(T.M[k]).sum();
    return sum;
  }
};

// Fast winding number coefficients for different orders of Taylor expansion
// For more information about the derivation of the coefficients,
// please refer to the original paper: Fast Winding Numbers for Soups and Clouds of Points
// by Barill et al., Siggraph 2018
template <class TriangleMesh, class GeomTraits = Default, int ORDER = 3> struct Fast_winding_number_Coeff;

template <class TriangleMesh, class GeomTraits>
struct Fast_winding_number_Coeff<TriangleMesh, GeomTraits, 1>
{
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  inline void add(const Vector_3& /* centered_centroid */, const Vector_3& /* normal */, const FT& /* area */) {}
  FT sum_area = 0;
  Point_3 weighted_centroid = Point_3(0, 0, 0);
  Vector_3 weighted_normal = Vector_3(0, 0, 0);
};

template <class TriangleMesh, class GeomTraits>
struct Fast_winding_number_Coeff<TriangleMesh, GeomTraits, 2>
{
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;
  using FT = typename GT::FT;
  using Vec3 = Eigen::Matrix<FT, 3, 1>;
  using Mat3 = Eigen::Matrix<FT, 3, 3>;

  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  inline void add(const Vector_3& centered_centroid, const Vector_3& normal, const FT& area) {
    Vec3 centroid_eigen(centered_centroid.x(), centered_centroid.y(), centered_centroid.z());
    Vec3 normal_eigen(normal.x(), normal.y(), normal.z());
    Q += area * (centroid_eigen * normal_eigen.transpose());
  }
  FT sum_area = 0;
  Point_3 weighted_centroid = Point_3(0, 0, 0);
  Vector_3 weighted_normal = Vector_3(0, 0, 0);
  Mat3 Q = Mat3::Zero();
};

template <class TriangleMesh, class GeomTraits>
struct Fast_winding_number_Coeff<TriangleMesh, GeomTraits, 3>
{
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;
  using FT = typename GT::FT;
  using Vec3 = Eigen::Matrix<FT, 3, 1>;
  using Mat3 = Eigen::Matrix<FT, 3, 3>;

  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;

  inline void add(const Vector_3& v1,
                  const Vector_3& v2,
                  const Vector_3& v3,
                  const Vector_3& centered_centroid,
                  const Vector_3& normal,
                  const FT& area) {
    Vec3 centroid_eigen(centered_centroid.x(), centered_centroid.y(), centered_centroid.z());
    Vec3 normal_eigen(normal.x(), normal.y(), normal.z());
    Q += area * (centroid_eigen * normal_eigen.transpose());
    T.add(v1, v2, v3, normal, area);
  }
  FT sum_area = 0;
  Point_3 weighted_centroid = Point_3(0, 0, 0);
  Vector_3 weighted_normal = Vector_3(0, 0, 0);
  Mat3 Q = Mat3::Zero();
  Tensor3<GT> T;
};

#endif // DOXYGEN_RUNNING

/**
 * \ingroup PkgVMASRef
 * \class CGAL::Fast_winding_number
 * \brief %Fast evaluation of the (normalized) winding number of a triangle mesh.
 *
 * This class implements a hierarchical (multipole / cluster–far field) approximation of
 * the generalized winding number for a triangle mesh.
 *
 * \tparam TriangleMesh
 *         a model of `FaceListGraph`
 * \tparam FaceNormalMap
 *         a readable property map contains the normal of each face: face_descriptor -> Vector_3
 * \tparam FaceAreaMap
 *         a readable property map contains the area of each face: face_descriptor -> FT
 * \tparam FaceCentroidMap
 *         a readable property map contains the centroid of each face: face_descriptor -> Point_3
 * \tparam Tree
 *         a built AABB tree over the faces of `TriangleMesh`
 * \tparam ORDER
 *         the order of the taylor expansion, must be 1, 2, or 3.<br>
 *        <b>%Default:</b> `3`
 * \tparam GeomTraits
 *        a model of `Kernel`<br>
 *       <b>%Default:</b>
 * \code
 *    CGAL::Kernel_traits<
 *     boost::property_traits<
 *       boost::property_map<TriangleMesh, CGAL::vertex_point_t>::type
 *    >::value_type
 *  >::Kernel
 * \endcode
 *         a model of `ReadWritePropertyMap`
 *       with `boost::graph_traits<TriangleMesh_>::%vertex_descriptor` as key and
 *       `GeomTraits_::Point_3` as value type.<br>
 *       <b>%Default:</b>
 * \code
 *   boost::property_map<TriangleMesh_, CGAL::vertex_point_t>::const_type.
 * \endcode
 *
 */
template <class TriangleMesh,
          class FaceNormalMap,
          class FaceAreaMap,
          class FaceCentroidMap,
          class Tree,
          int ORDER = 3,
          class GeomTraits = Default,
          class VertexPointMap = Default>
class Fast_winding_number
{
  using VPM = typename Default::Get<VertexPointMap,
                                    typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::type;
  using GT = typename Default::Get<
      GeomTraits,
      typename Kernel_traits<typename boost::property_traits<
          typename boost::property_map<TriangleMesh, vertex_point_t>::type>::value_type>::Kernel>::type;
  using FT = typename GT::FT;
  using Point_3 = typename GT::Point_3;
  using Vector_3 = typename GT::Vector_3;
  using Sphere_3 = typename GT::Sphere_3;
  using Coeff = Fast_winding_number_Coeff<TriangleMesh, GeomTraits, ORDER>;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using Node = CGAL::AABB_node<typename Tree::AABB_traits>;
  using Vec3 = Eigen::Matrix<FT, 3, 1>;
  using Mat3 = Eigen::Matrix<FT, 3, 3>;

public:
  /** \name Constructor
   *@{
   * The constructor of a vmas object.
   *
   * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
   *
   * @param tmesh
   *        The input triangle mesh without borders.
   * @param fnm
   *        A readable property map that maps each face of `tmesh` to its normal vector.
   * @param fam
   *        A readable property map that maps each face of `tmesh` to its area.
   * @param fcm
   *        A readable property map that maps each face of `tmesh` to its centroid.
   * @param tree
   *        A built AABB tree over the faces of `tmesh`.
   * @param np
   *        An optional sequence of \ref bgl_namedparameters "Named Parameters", listed below:
   *
   * \cgalNamedParamsBegin
   *   \cgalParamNBegin{vertex_point_map}
   *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
   *     \cgalParamType{a class model of `ReadablePropertyMap` with
   *     `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
   *                    as key type and `%Point_3` as value type}
   *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
   *   \cgalParamNEnd
   *   \cgalParamNBegin{fast_winding_number_beta}
   *     \cgalParamDescription{The parameter to control the accuracy/speed trade-off of the fast winding number
   *     evaluation. little values lead to more accurate results, but slower evaluation.}
   *     \cgalParamType{FT}
   *     \cgalParamDefault{FT(2.0)}
   *     \cgalParamExtra{The range of this parameter is (0,+inf).}
   *   \cgalParamNEnd
   * \cgalNamedParamsEnd
   */
  template <class NamedParameters = parameters::Default_named_parameters>
  Fast_winding_number(const TriangleMesh& tmesh,
                      const FaceNormalMap& fnm,
                      const FaceAreaMap& fam,
                      const FaceCentroidMap& fcm,
                      Tree& tree,
                      const NamedParameters& np = parameters::default_values()
                      )
      : tmesh_(tmesh)
      , fnm_(fnm)
      , fam_(fam)
      , fcm_(fcm)
      , tree_(tree)
      {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    vpm_ = choose_parameter(get_parameter(np, internal_np::vertex_point),
                            get_const_property_map(CGAL::vertex_point, tmesh));
    beta_ = choose_parameter(get_parameter(np, internal_np::fast_winding_number_beta), FT(2.0));
    std::size_t nb_nodes = tree_.nb_node();
    node_radius_.resize(nb_nodes, FT(0));
    coeffs_.resize(nb_nodes);
    precompute_coeffs();
  }
  ///@}

  /**
   *
   * \brief computes the fast winding number of a given query point.
   *
   * This function computes the fast winding number of a given query point `p` recursively.
   * For a given query point `p`, we first compute its distance `r` to the center of root node in the AABB tree.
   * If `r` is larger than `beta` times the radius of the root node, we consider the root node as a single dipole
   * and use the precomputed Taylor expansion coefficients to evaluate the winding number.
   * Otherwise, we traverse down the tree recursively and repeat the same process for each child node.
   *
   * @param p the query point
   * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters", listed below:
   * \cgalNamedParamsBegin
   *  \cgalParamNBegin{fast_winding_number_beta}
   *   \cgalParamDescription{The parameter to control the accuracy/speed trade-off of the fast winding number
   *  evaluation. little values lead to more accurate results, but slower evaluation.}
   *  \cgalParamType{FT}
   *  \cgalParamDefault{FT(2.0)}
   * \cgalParamExtra{The range of this parameter is (0,+inf).}
   *
   * @return the fast winding number of `p`
   */
  template <class NamedParameters = parameters::Default_named_parameters>
  FT evaluate_fast_winding_number(const Point_3& p, const NamedParameters& np = parameters::default_values())  {
    using parameters::choose_parameter;
    using parameters::get_parameter;
    beta_ = choose_parameter(get_parameter(np, internal_np::fast_winding_number_beta), FT(2.0));

    FT winding_number = FT(0);
    const Node* root = tree_.root_node();
    std::size_t nb_prim = tree_.size();
    std::vector<std::tuple<const Node*, typename Tree::Primitive_iterator, std::size_t>> traversal_queue;
    traversal_queue.emplace_back(root, tree_.primitives_begin(), nb_prim);
    while(!traversal_queue.empty()) {
      auto [node, prim_it, nb_prim] = traversal_queue.back();
      traversal_queue.pop_back();

      std::size_t n_id = node_id(node);
      Vector_3 R = coeffs_[n_id].weighted_centroid - p;
      FT r = CGAL::approximate_sqrt(R.squared_length());
      if(r > beta_ * node_radius_[n_id]) {
        // The query point is far enough from the node, so use single dipole approximation
        winding_number += direct_eval(node, p);
      } else {
        switch(nb_prim) {
        case 2: {
          // leaf node
          winding_number += solid_angle_leaf(node, p);
        } break;
        case 3: {
          // partial leaf node
          for(auto it = prim_it; it < prim_it + nb_prim; ++it) {
            face_descriptor f = it->id();
            winding_number += solid_angle(p, f) / (FT(4) * CGAL_PI);
          }

        } break;
        default: {
          std::size_t nb_left = nb_prim / 2;
          std::size_t nb_right = nb_prim - nb_left;

          traversal_queue.emplace_back(std::addressof(node->left_child()), prim_it, nb_left);
          traversal_queue.emplace_back(std::addressof(node->right_child()), prim_it + nb_left, nb_right);
        } break;
        }
      }
    }
    return winding_number;
  }

  /**
   * \brief computes the exact winding number of a given query point `p` and the input mesh.
   *
   * @param p the query point
   * @return the exact winding number of `p`
   */
  FT exact_winding_number(const Point_3& p) const {
    FT winding_number = FT(0);
    for(face_descriptor f : faces(tmesh_)) {
      winding_number += solid_angle(p, f);
    }
    return winding_number / (FT(4) * CGAL_PI);
  }

  /**
   * \brief tests whether a given query point is inside the input mesh.
   *
   * @param p the query point
   * @return `true` if `p` is inside the input mesh, `false` otherwise.
   */
  bool is_inside(const Point_3& p) { return evaluate_fast_winding_number(p) > FT(0.5); }

private:
  std::size_t node_id(const Node* node) const { return std::size_t(node - tree_.root_node()); }

  //This function precommutes the coefficients for all nodes in the AABB tree.
  //We search from top to bottom, and for each node, we compute the coefficients
  //by looping over all of the faces it contains.
  void precompute_coeffs() {

    const Node* root = tree_.root_node();

    std::size_t nb_primitives = tree_.size();
    std::vector<std::tuple<const Node*, typename Tree::Primitive_iterator, std::size_t>> traversal_queue;
    traversal_queue.emplace_back(root, tree_.primitives_begin(), nb_primitives);

    while(!traversal_queue.empty()) {
      auto [node, prim_it, nb_prim] = traversal_queue.back();
      traversal_queue.pop_back();

      Coeff& coeff = coeffs_[node_id(node)];
      Vector_3 sum_centroid(0, 0, 0);
      Vector_3 sum_normal(0, 0, 0);
      FT total_area = FT(0);

      //First loop:
      //for each node, we loop over all of the faces it contains and compute the weighted centroid, weighted normal, and total area of the each node
      for(auto it = prim_it; it < prim_it + nb_prim; ++it) {
        face_descriptor f = it->id();
        FT area = get(fam_, f);
        const Vector_3& normal = get(fnm_, f);
        const Point_3& centroid = get(fcm_, f);

        sum_centroid += (centroid - CGAL::ORIGIN) * area;
        sum_normal += area * normal;
        total_area += area;
      }
      CGAL_assertion(total_area > FT(0));
      coeff.weighted_centroid =
          Point_3(sum_centroid.x() / total_area, sum_centroid.y() / total_area, sum_centroid.z() / total_area);
      //Second loop:
      //we loop over all of the faces again to compute the higher order coefficients
      //The reason we do not compute the higher order coefficients in the first loop is that
      //we need the weighted centroid of the node (vector R is centralized)
      FT max_norm_sq = 0;
      for(auto it = prim_it; it < prim_it + nb_prim; ++it) {
        face_descriptor f = it->id();
        FT area = get(fam_, f);
        const Vector_3& normal = get(fnm_, f);
        const Point_3& centroid = get(fcm_, f);
        Vector_3 R = centroid - coeff.weighted_centroid;

        for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh_), tmesh_)) {
          const Point_3& pv = get(vpm_, v);
          max_norm_sq = (std::max)(max_norm_sq, (pv - coeff.weighted_centroid).squared_length());
        }
        if constexpr(ORDER < 3)
          coeff.add(R, normal, area);
        else {
          Vector_3 iv[3];
          int i = 0;
          for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh_), tmesh_)) {
            const Point_3& pv = get(vpm_, v);
            iv[i] = Vector_3(pv.x(), pv.y(), pv.z());
            ++i;
          }
          //See Appendix_B in the original paper for more details about this formula
          Vector_3 v1 = 0.5 * (iv[0] + iv[1]) - (coeff.weighted_centroid - CGAL::ORIGIN);
          Vector_3 v2 = 0.5 * (iv[1] + iv[2]) - (coeff.weighted_centroid - CGAL::ORIGIN);
          Vector_3 v3 = 0.5 * (iv[2] + iv[0]) - (coeff.weighted_centroid - CGAL::ORIGIN);
          coeff.add(v1, v2, v3, R, normal, area);
        }
      }
      node_radius_[node_id(node)] = CGAL::approximate_sqrt(max_norm_sq);
      coeff.sum_area = total_area;
      coeff.weighted_normal = sum_normal;

      switch(nb_prim) {
      case 2: {
        // leaf node
      } break;
      case 3: {
        // partial leaf node
        // left contains data, right is a node with the two remaining primitives
        traversal_queue.emplace_back(std::addressof(node->right_child()), prim_it + 1, std::size_t(2));
      } break;
      default: {
        std::size_t nb_left = nb_prim / 2;
        std::size_t nb_right = nb_prim - nb_left;

        traversal_queue.emplace_back(std::addressof(node->left_child()), prim_it, nb_left);
        traversal_queue.emplace_back(std::addressof(node->right_child()), prim_it + nb_left, nb_right);
      } break;
      }
    }
  }

  // This function computes the solid angle of a triangle face `f` seen from a point `p`.
  FT solid_angle(const Point_3& p, const face_descriptor& f) const {
    Vector_3 iv[3];
    int i = 0;
    for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh_), tmesh_)) {
      iv[i] = get(vpm_, v) - p;
      ++i;
    }
    FT l[3];
    for(std::size_t j = 0; j < 3; ++j)
      l[j] = CGAL::approximate_sqrt(iv[j].squared_length());
    if(l[0] < 1e-15 || l[1] < 1e-15 || l[2] < 1e-15) {
      return 0;
    }
    for(std::size_t j = 0; j < 3; ++j)
      iv[j] = iv[j] / l[j];
    FT numerator = CGAL::scalar_product(iv[0], CGAL::cross_product((iv[1] - iv[0]), iv[2] - iv[1]));
    if(numerator == 0)
      return 0;
    FT denominator = FT(1) + CGAL::scalar_product(iv[0], iv[1]) + CGAL::scalar_product(iv[1], iv[2]) +
                     CGAL::scalar_product(iv[2], iv[0]);
    return 2 * std::atan2(numerator, denominator);
  }

  //helper function to compute the solid angle of a leaf node
  FT solid_angle_leaf(const Node* node, const Point_3& p) const {
    FT w = FT(0);
    face_descriptor f1 = node->left_data().id();
    w += solid_angle(p, f1);
    face_descriptor f2 = node->right_data().id();
    w += solid_angle(p, f2);
    return w / (4 * CGAL_PI);
  }

  // helper function to evaluate the node as a single dipole.
  // This is the case when the query point is far enough from the node.
  FT direct_eval(const Node* node, const Point_3& p) const {
    FT w = FT(0);
    const Vector_3 R = coeffs_[node_id(node)].weighted_centroid - p;
    Vec3 R_eigen(R.x(), R.y(), R.z());
    FT r = (std::max)(CGAL::approximate_sqrt(R.squared_length()), FT(1e-15));
    if constexpr(ORDER >= 1) {
      Vector_3 G1 = R / (4 * CGAL_PI * r * r * r);
      w += CGAL::scalar_product(G1, coeffs_[node_id(node)].weighted_normal);
    }
    if constexpr(ORDER >= 2) {
      FT r2 = r * r;
      FT r3 = r2 * r;
      FT r5 = r2 * r3;
      Mat3 G2 = (r2 * Mat3::Identity() - 3 * R_eigen * R_eigen.transpose()) / (4 * CGAL_PI * r5);
      w += (G2.cwiseProduct(coeffs_[node_id(node)].Q)).sum();
    }
    if constexpr(ORDER >= 3) {
      Tensor3<GT> K;
      const FT c1 = 15.0 / (4 * CGAL_PI * std::pow(r, 7));
      const FT c2 = 3.0 / (4 * CGAL_PI * std::pow(r, 5));
      for(int k = 0; k < 3; ++k) {

        for(int i = 0; i < 3; ++i) {
          for(int j = 0; j < 3; ++j) {

            K.M[k](i, j) = c1 * R_eigen(i) * R_eigen(j) * R_eigen(k) -
                           c2 * (R_eigen(k) * (i == j) + R_eigen(j) * (i == k) + R_eigen(i) * (k == j));
          }
        }
      }
      w += K * coeffs_[node_id(node)].T;
    }
    return w;
  }

private:
  const TriangleMesh& tmesh_;
  const FaceNormalMap& fnm_;
  const FaceAreaMap& fam_;
  const FaceCentroidMap& fcm_;
  VPM vpm_;

  std::vector<FT> node_radius_;
  std::vector<Coeff> coeffs_;
  Tree& tree_;
  FT beta_;
};

} // end of namespace CGAL

#endif // CGAL_FAST_WINDING_NUMBER_H
