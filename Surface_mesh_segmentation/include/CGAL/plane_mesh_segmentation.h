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

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

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
ECM get_ecm(ECM ecm, const TriangleMesh& tm)
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

} // namespace Planar_segmentation

/*!
 * @TODO: Add doc with vpm, ecm, and plane per segment property map
 * ecm is both in and out
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
  double min_cosinus_squared = choose_parameter<double>(get_parameter(np, internal_np::min_cosinus_squared), 1);

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

template <class TriangleMesh, class SegmentPropertyMap>
std::size_t
segment_via_plane_fitting(const TriangleMesh& tm, SegmentPropertyMap segment_ids)
{
  return segment_via_plane_fitting(tm, segment_ids, parameters::all_default());
}

} // namespace CGAL

#endif
