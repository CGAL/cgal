// Copyright (c) 2021  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_PLANE_MESH_SEGMENTATION_H
#define CGAL_PLANE_MESH_SEGMENTATION_H

#include <CGAL/license/Surface_mesh_segmentation.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Filtered_predicate.h>
#include <CGAL/Exact_kernel_selector.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

namespace CGAL {

/// @TODO these predicates should probably end-up in the Kernel
namespace Predicates {

// predicates
template <typename K>
struct Is_dihedral_angle_close_to_Pi_impl{
  typedef bool result_type;

  /// Compares the value of the dihedral angle between the triangles `pqr` and `pqs`
  /// with the angles \$f\theta\f$ and \f$\pi + \theta\f$.
  /// `false` is returned if the dihedral angle is less or equal to \f$\pi\f$.
  /// If the cosinus of the dihedral angle is smaller or equal to the cosinus of \f$\theta\f$
  /// `true` is returned and `false` otherwise. \f$theta\f$ is provided using the square of its
  /// cosinus.
  bool
  operator()(const typename K::Point_3& p,
             const typename K::Point_3& q,
             const typename K::Point_3& r,
             const typename K::Point_3& s,
             const typename K::FT min_cosinus_squared) const
  {
    const typename K::Vector_3 pq(p,q);
    const typename K::Vector_3 pr(p,r);
    const typename K::Vector_3 ps(p,s);

    typename K::Construct_cross_product_vector_3 cross_product =
      K().construct_cross_product_vector_3_object();
    const typename K::Vector_3 pqpr = cross_product(pq, pr);
    const typename K::Vector_3 pqps = cross_product(pq, ps);
    const typename K::FT sc_prod_1 = pqpr * pqps;

    typename K::Compute_squared_length_3 squared_length = K().compute_squared_length_3_object();
    if (sign(sc_prod_1) != NEGATIVE) return false;
    return compare( square(sc_prod_1),
                    min_cosinus_squared
                      * squared_length(pqpr)
                      * squared_length(pqps) ) != SMALLER;
  }
};

template <typename K, template <class Kernel> class Pred>
struct Get_filtered_predicate_RT
{
  typedef typename Exact_kernel_selector<K>::Exact_kernel_rt Exact_kernel_rt;
  typedef typename Exact_kernel_selector<K>::C2E_rt C2E_rt;
  typedef Simple_cartesian<Interval_nt_advanced>        Approximate_kernel;
  typedef Cartesian_converter<K, Approximate_kernel>   C2F;


  typedef Filtered_predicate<Pred<Exact_kernel_rt>,
                             Pred<Approximate_kernel>,
                             C2E_rt, C2F> type;
};

template<typename K, bool has_filtered_predicates = K::Has_filtered_predicates>
struct Is_dihedral_angle_close_to_Pi
  : public Is_dihedral_angle_close_to_Pi_impl<K>
{
  using Is_dihedral_angle_close_to_Pi_impl<K>::operator();
};

template<typename K>
struct Is_dihedral_angle_close_to_Pi<K, true>
  : public Get_filtered_predicate_RT<K, Is_dihedral_angle_close_to_Pi_impl>::type
{
  using Get_filtered_predicate_RT<K, Is_dihedral_angle_close_to_Pi_impl>::type::operator();
};

} // namespace Predicates

namespace Planar_segmentation {

template <typename TriangleMesh,
          typename VertexPointMap,
          typename edge_descriptor>
bool is_edge_between_coplanar_faces(edge_descriptor e,
                                    const TriangleMesh& tm,
                                    double min_cosinus_squared,
                                    const VertexPointMap& vpm)
{
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::type K;
  typedef typename boost::property_traits<VertexPointMap>::reference Point_ref_3;
  if (is_border(e, tm)) return false;
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
    h =  halfedge(e, tm);
  Point_ref_3 p = get(vpm, source(h, tm) );
  Point_ref_3 q = get(vpm, target(h, tm) );
  Point_ref_3 r = get(vpm, target(next(h, tm), tm) );
  Point_ref_3 s = get(vpm, target(next(opposite(h, tm), tm), tm) );

  if (min_cosinus_squared==1)
    return coplanar(p, q, r, s);
  else
  {
    Predicates::Is_dihedral_angle_close_to_Pi<K> pred;
    return pred(p, q, r, s, min_cosinus_squared);
  }
}

template <typename TriangleMesh,
          typename EdgeIsConstrainedMap,
          typename VertexPointMap>
void
mark_constrained_edges(
  TriangleMesh& tm,
  EdgeIsConstrainedMap edge_is_constrained,
  double min_cosinus_squared,
  const VertexPointMap& vpm)
{
  for(typename boost::graph_traits<TriangleMesh>::edge_descriptor e : edges(tm))
  {
    if (!get(edge_is_constrained,e))
      if (!is_edge_between_coplanar_faces(e, tm, min_cosinus_squared, vpm))
        put(edge_is_constrained, e, true);
  }
}

template <typename ECM, typename TriangleMesh>
ECM get_ecm(ECM ecm, const TriangleMesh& /*tm*/)
{
  return ecm;
}

template <typename ECM, typename TriangleMesh>
const ECM& get_ecm(const internal_np::Param_not_found&, const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  ECM ecm = get(CGAL::dynamic_edge_property_t<bool>(), tm);

  for(edge_descriptor e : edges(tm))
    put(ecm, e, false);
  return ecm;
}

#ifndef CGAL_DO_NOT_USE_PCA
template <typename TriangleMesh,
          typename FaceCCIdMap,
          typename VertexPointMap>
typename boost::property_traits<FaceCCIdMap>::value_type
coplanarity_segmentation_with_region_growing(const TriangleMesh& triangle_mesh,
                                             const double max_distance,
                                             const double cos_value_squared,
                                             FaceCCIdMap& face_cc_ids,
                                             const VertexPointMap& vpm)
{
  using Kernel = typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel;
  using CC_ID = typename boost::property_traits<FaceCCIdMap>::value_type;
  using FT = typename Kernel::FT;

  // Get face range and types.
  const auto face_range = faces(triangle_mesh);

  using Face_range          = decltype(face_range);
  using Triangle_mesh       = TriangleMesh;
  using Vertex_to_point_map = VertexPointMap;

  using Neighbor_query = CGAL::Shape_detection::Polygon_mesh::
    One_ring_neighbor_query<Triangle_mesh, Face_range>;
  using Region_type    = CGAL::Shape_detection::Polygon_mesh::
    Least_squares_plane_fit_region<Kernel, Triangle_mesh, Face_range, Vertex_to_point_map>;
  using Sorting        = CGAL::Shape_detection::Polygon_mesh::
    Least_squares_plane_fit_sorting<Kernel, Triangle_mesh, Neighbor_query, Face_range, Vertex_to_point_map>;

  using Region_growing = CGAL::Shape_detection::
    Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(triangle_mesh);

  Region_type region_type(
    triangle_mesh,
    CGAL::parameters::
    maximum_distance(static_cast<FT>(max_distance)).
    cosine_value(static_cast<FT>(CGAL::sqrt(cos_value_squared))).
    vertex_point_map(vpm));

  // Sort face indices.
  Sorting sorting(
    triangle_mesh, neighbor_query, CGAL::parameters::vertex_point_map(vpm));
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    face_range, neighbor_query, region_type, sorting.seed_map());

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  // std::cout << "- found regions: " << regions.size() << std::endl;

  for (std::size_t i = 0; i < regions.size(); ++i) {
    for (const std::size_t face_index : regions[i]) {
      const auto& face = *(face_range.begin() + face_index);
      put(face_cc_ids, face, static_cast<CC_ID>(i));
    }
  }
  return static_cast<CC_ID>(regions.size());
}

template <typename TriangleMesh,
          typename FaceCCIdMap,
          typename VertexPointMap>
typename boost::property_traits<FaceCCIdMap>::value_type
coplanarity_segmentation_with_pca(TriangleMesh& tm,
                                  double max_frechet_distance,
                                  double min_cosinus_squared,
                                  FaceCCIdMap& face_cc_ids,
                                  const VertexPointMap& vpm)
{
  typedef typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel IK; // input kernel
  typedef CGAL::Exact_predicates_inexact_constructions_kernel PCA_K;
  typedef boost::graph_traits<TriangleMesh> graph_traits;

  typedef typename boost::property_traits<FaceCCIdMap>::value_type CC_ID;

  std::size_t nb_faces = std::distance(faces(tm).first, faces(tm).second);
  std::size_t faces_tagged = 0;
  CC_ID cc_id(-1);
  typename graph_traits::face_iterator fit_seed = boost::begin(faces(tm));

  CGAL::Cartesian_converter<IK, PCA_K> to_pca_k;

  const double max_squared_frechet_distance = max_frechet_distance * max_frechet_distance;

  auto get_triangle = [&tm, &vpm, &to_pca_k](typename graph_traits::face_descriptor f)
  {
    typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
      h = halfedge(f, tm);
    return typename PCA_K::Triangle_3(to_pca_k(get(vpm, source(h, tm))),
                                      to_pca_k(get(vpm, target(h, tm))),
                                      to_pca_k(get(vpm, target(next(h, tm), tm))));
  };

  while(faces_tagged!=nb_faces)
  {
    CGAL_assertion(faces_tagged<=nb_faces);

    while( get(face_cc_ids, *fit_seed)!=CC_ID(-1) )
    {
      ++fit_seed;
      CGAL_assertion( fit_seed!=boost::end(faces(tm)) );
    }
    typename graph_traits::face_descriptor seed = *fit_seed;

    std::vector<typename graph_traits::halfedge_descriptor> queue; /// not sorted for now @TODO
    queue.push_back( halfedge(seed, tm) );
    queue.push_back( next(queue.back(), tm) );
    queue.push_back( next(queue.back(), tm) );

    std::vector<typename PCA_K::Triangle_3> current_selection; /// @TODO compare lls-fitting with only points
    current_selection.push_back(get_triangle(seed));
    put(face_cc_ids, seed, ++cc_id);
    ++faces_tagged;

    auto does_fitting_respect_distance_bound = [&to_pca_k, &vpm, max_squared_frechet_distance](
      std::unordered_set<typename graph_traits::vertex_descriptor>& vertices,
      const PCA_K::Plane_3& plane)
    {
      typename PCA_K::Compare_squared_distance_3 compare_squared_distance;
      for (typename graph_traits::vertex_descriptor v : vertices)
      {
        if (compare_squared_distance(to_pca_k(get(vpm, v)), plane, max_squared_frechet_distance) == LARGER)
          return false;
      }
      return true;
    };

    std::unordered_set<typename graph_traits::vertex_descriptor> vertex_selection;
    vertex_selection.insert(target(queue[0], tm));
    vertex_selection.insert(target(queue[1], tm));
    vertex_selection.insert(target(queue[2], tm));
    while(!queue.empty())
    {
      typename graph_traits::halfedge_descriptor h = queue.back(),
                                               opp = opposite(h, tm);
      queue.pop_back();
      if (is_border(opp, tm) || get(face_cc_ids, face(opp, tm))!=CC_ID(-1)) continue;
      if (!CGAL::Planar_segmentation::is_edge_between_coplanar_faces(edge(h, tm), tm, min_cosinus_squared, vpm) )
        continue;
      current_selection.push_back(get_triangle(face(opp, tm)));

      bool new_vertex_added = vertex_selection.insert(target(next(opp, tm), tm)).second;

      typename PCA_K::Plane_3 plane;
      typename PCA_K::Point_3 centroid;

      linear_least_squares_fitting_3(current_selection.begin(),
                                     current_selection.end(),
                                     plane,
                                     centroid,
                                     Dimension_tag<2>());

      if (!new_vertex_added || does_fitting_respect_distance_bound(vertex_selection, plane))
      {
        put(face_cc_ids, face(opp, tm), cc_id);
        ++faces_tagged;
        queue.push_back(next(opp, tm));
        queue.push_back(prev(opp, tm));
      }
      else{
        /// @TODO add an opti to avoid testing several time a face rejected
        current_selection.pop_back();
        vertex_selection.erase(target(next(opp, tm), tm));
      }
    }
  }

  return cc_id+1;
}
#endif

} // namespace Planar_segmentation

/*!
 * @TODO: Add doc with vpm, ecm, and plane per segment property map
 * ecm is both in and out
 * @TODO ecm as input is ignored if PCA is used...
 */
template <class TriangleMesh, class SegmentPropertyMap, class NamedParameters>
typename boost::property_traits<SegmentPropertyMap>::value_type
segment_via_plane_fitting(const TriangleMesh& tm, SegmentPropertyMap segment_ids, const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh,
                                       CGAL::dynamic_edge_property_t<bool> >::type Default_ecm;

  // User or default edge is-constrained map
  typedef typename internal_np::Lookup_named_param_def <
    internal_np::edge_is_constrained_t,
    NamedParameters,
    Default_ecm
  > ::type Ecm;

  Ecm ecm = Planar_segmentation::get_ecm<Ecm>(get_parameter(np, internal_np::edge_is_constrained), tm);

  // Vpm
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters>::const_type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  // input parameter for example (exact or approximate)
  const double min_cosinus_squared = choose_parameter<double>(get_parameter(np, internal_np::min_cosinus_squared), 1);
  const double max_frechet_distance = choose_parameter<double>(get_parameter(np, internal_np::max_Frechet_distance), 33);
  const bool use_region_growing = choose_parameter<bool>(get_parameter(np, internal_np::use_region_growing), true);

  if (max_frechet_distance == 33 && !use_region_growing) // default parameter value
  {

    // mark constrained edges
    /// @TODO: hardcoded path if min_cosinus_squared == 1
    Planar_segmentation::mark_constrained_edges(tm, ecm, min_cosinus_squared, vpm);

    // segment connected components (cc) delimited by constrained edges
    typename boost::property_traits<SegmentPropertyMap>::value_type
      res = Polygon_mesh_processing::connected_components(
        tm, segment_ids, parameters::edge_is_constrained_map(ecm));

    // for consistency in case approximate angle is used
    if (min_cosinus_squared!=1)
    {
      for(edge_descriptor e : edges(tm))
      {
        if (get(ecm, e) && !is_border(e, tm))
        {
          halfedge_descriptor h = halfedge(e, tm);
          if ( get(segment_ids, face(h, tm))==get(segment_ids, face(opposite(h, tm), tm)) )
            put(ecm, e, false);
        }
      }
    }
    return res;
  }
  else
  {
    typename boost::property_traits<SegmentPropertyMap>::value_type res;
    if (!use_region_growing) {
      // use PCA version
      res = Planar_segmentation::coplanarity_segmentation_with_pca(tm,
                                                                   max_frechet_distance,
                                                                   min_cosinus_squared,
                                                                   segment_ids,
                                                                   vpm);
    } else {
      // use Region Growing version
      res = Planar_segmentation::coplanarity_segmentation_with_region_growing(tm,
                                                                              max_frechet_distance,
                                                                              min_cosinus_squared,
                                                                              segment_ids,
                                                                              vpm);
    }

    // mark constrained edges
    for(edge_descriptor e : edges(tm))
    {
      halfedge_descriptor h = halfedge(e, tm);
      if (is_border(e, tm) || get(segment_ids, face(h, tm)) !=
                              get(segment_ids, face(opposite(h, tm), tm)))
      {
        put(ecm, e, true);
      }
    }
    return res;
  }
}

template <class TriangleMesh, class SegmentPropertyMap>
typename boost::property_traits<SegmentPropertyMap>::value_type
segment_via_plane_fitting(const TriangleMesh& tm, SegmentPropertyMap segment_ids)
{
  return segment_via_plane_fitting(tm, segment_ids, parameters::all_default());
}

} // namespace CGAL

#endif
