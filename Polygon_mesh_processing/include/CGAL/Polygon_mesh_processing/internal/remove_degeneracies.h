// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMOVE_DEGENERACIES_H
#define CGAL_POLYGON_MESH_PROCESSING_REMOVE_DEGENERACIES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <set>
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
#include <sstream>
#include <fstream>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename TriangleMesh, typename VPM, typename ECM, typename Traits>
std::array<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor, 2>
is_badly_shaped(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                TriangleMesh& tmesh,
                const VPM& vpm,
                const ECM& ecm,
                const Traits& gt,
                const double cap_threshold, // angle over 160° ==> cap
                const double needle_threshold, // longest edge / shortest edge over this ratio ==> needle
                const double collapse_length_threshold) // max length of edges allowed to be collapsed
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

  const halfedge_descriptor null_h = boost::graph_traits<TriangleMesh>::null_halfedge();

  halfedge_descriptor res = PMP::is_needle_triangle_face(f, tmesh, needle_threshold,
                                                         parameters::vertex_point_map(vpm)
                                                                    .geom_traits(gt));
  if(res != null_h && !get(ecm, edge(res, tmesh)))
  {
    // don't want to collapse edges that are too large
    if(collapse_length_threshold == 0  ||
       edge_length(res, tmesh, parameters::vertex_point_map(vpm).geom_traits(gt)) <= collapse_length_threshold)
    {
      return make_array(res, null_h);
    }
  }
  else // let's not make it possible to have a face be both a cap and a needle (for now)
  {
    res = PMP::is_cap_triangle_face(f, tmesh, cap_threshold, parameters::vertex_point_map(vpm).geom_traits(gt));
    if(res != null_h && !get(ecm, edge(res, tmesh)))
      return make_array(null_h, res);
  }

  return make_array(null_h, null_h);
}

template <typename TriangleMesh, typename EdgeContainer,
          typename VPM, typename ECM, typename Traits>
void collect_badly_shaped_triangles(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                    TriangleMesh& tmesh,
                                    const VPM& vpm,
                                    const ECM& ecm,
                                    const Traits& gt,
                                    const double cap_threshold, // angle over this threshold (as a cosine) ==> cap
                                    const double needle_threshold, // longest edge / shortest edge over this ratio ==> needle
                                    const double collapse_length_threshold, // max length of edges allowed to be collapsed
                                    EdgeContainer& edges_to_collapse,
                                    EdgeContainer& edges_to_flip)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

  std::array<halfedge_descriptor, 2> res = is_badly_shaped(f, tmesh, vpm, ecm, gt, cap_threshold,
                                                           needle_threshold, collapse_length_threshold);

  if(res[0] != boost::graph_traits<TriangleMesh>::null_halfedge())
  {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
    std::cout << "add new needle: " << edge(res[0], tmesh) << std::endl;
#endif
    edges_to_collapse.insert(edge(res[0], tmesh));
  }
  else // let's not make it possible to have a face be both a cap and a needle (for now)
  {
    if(res[1] != boost::graph_traits<TriangleMesh>::null_halfedge())
    {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "add new cap: " << edge(res[1],tmesh) << std::endl;
#endif
      edges_to_flip.insert(edge(res[1], tmesh));
    }
  }
}

/*
// Following Ronfard et al. 96 we look at variation of the normal after the collapse
// the collapse must be topologically valid
template <class TriangleMesh, class NamedParameters>
bool is_collapse_geometrically_valid(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                                     const TriangleMesh& tmesh,
                                     const NamedParameters& np)
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type   VPM;
  typedef typename boost::property_traits<VPM>::reference                         Point_ref;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             Traits;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, tmesh));
  Traits gt = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

/// @todo handle boundary edges

  h = opposite(h, tmesh); // Euler::collapse edge keeps the target and removes the source

  // source is kept, target is removed
  CGAL_assertion(target(h, tmesh) == vertex_removed);
  Point_ref kept = get(vpm, source(h, tmesh));
  Point_ref removed= get(vpm, target(h, tmesh));

  // consider triangles incident to the vertex removed
  halfedge_descriptor stop = prev(opposite(h, tmesh), tmesh);
  halfedge_descriptor hi = opposite(next(h, tmesh), tmesh);

  std::vector<halfedge_descriptor> triangles;
  while(hi != stop)
  {
    if(!is_border(hi, tmesh))
    {
      Point_ref a = get(vpm, target(next(hi, tmesh), tmesh));
      Point_ref b = get(vpm, source(hi, tmesh));

      //ack a-b-point_remove and a-b-point_kept has a compatible orientation
      /// @todo use a predicate
      typename Traits::Vector_3 n1 = gt.construct_cross_product_vector_3_object()(removed-a, b-a);
      typename Traits::Vector_3 n2 = gt.construct_cross_product_vector_3_object()(kept-a, b-a);
      if(gt.compute_scalar_product_3_object()(n1, n2) <= 0)
        return false;
    }

    hi = opposite(next(hi, tmesh), tmesh);
  }

  return true;
}
*/

template <class TriangleMesh, typename VPM, typename Traits>
boost::optional<double> get_collapse_volume(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                                            const TriangleMesh& tmesh,
                                            const VPM& vpm,
                                            const Traits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                         Point_ref;
  typedef typename Traits::Vector_3                                               Vector_3;

  const typename Traits::Point_3 origin(ORIGIN);

/// @todo handle boundary edges

  h = opposite(h, tmesh); // Euler::collapse edge keeps the target and removes the source

  // source is kept, target is removed
  Point_ref kept = get(vpm, source(h, tmesh));
  Point_ref removed= get(vpm, target(h, tmesh));

  // init volume with incident triangles (reversed orientation
  double delta_vol = volume(removed, kept, get(vpm, target(next(h, tmesh), tmesh)), origin) +
                     volume(kept, removed, get(vpm, target(next(opposite(h, tmesh), tmesh), tmesh)), origin);

  // consider triangles incident to the vertex removed
  halfedge_descriptor stop = prev(opposite(h, tmesh), tmesh);
  halfedge_descriptor hi = opposite(next(h, tmesh), tmesh);

  std::vector<halfedge_descriptor> triangles;
  while(hi != stop)
  {
    if(!is_border(hi, tmesh))
    {
      Point_ref a = get(vpm, target(next(hi, tmesh), tmesh));
      Point_ref b = get(vpm, source(hi, tmesh));

      //ack a-b-point_remove and a-b-point_kept has a compatible orientation
      /// @todo use a predicate
      Vector_3 v_ab = gt.construct_vector_3_object()(a, b);
      Vector_3 v_ar = gt.construct_vector_3_object()(a, removed);
      Vector_3 v_ak = gt.construct_vector_3_object()(a, kept);

      Vector_3 n1 = gt.construct_cross_product_vector_3_object()(v_ar, v_ab);
      Vector_3 n2 = gt.construct_cross_product_vector_3_object()(v_ak, v_ab);
      if(gt.compute_scalar_product_3_object()(n1, n2) <= 0)
        return boost::none;

      delta_vol += volume(b, a, removed, origin) + volume(a, b, kept, origin); // opposite orientation
    }

    hi = opposite(next(hi, tmesh), tmesh);
  }

  return CGAL::abs(delta_vol);
}

template <typename TriangleMesh, typename VPM, typename VCM, typename Traits>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
get_best_edge_orientation(typename boost::graph_traits<TriangleMesh>::edge_descriptor e,
                          const TriangleMesh& tmesh,
                          const VPM& vpm,
                          const VCM& vcm,
                          const Traits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = halfedge(e, tmesh), ho = opposite(h, tmesh);

  CGAL_assertion(!get(vcm, source(h, tmesh)) || !get(vcm, target(h, tmesh)));

  boost::optional<double> dv1 = get_collapse_volume(h, tmesh, vpm, gt);
  boost::optional<double> dv2 = get_collapse_volume(ho, tmesh, vpm, gt);

  // the resulting point of the collapse of a halfedge is the target of the halfedge before collapse
  if(get(vcm, source(h, tmesh)))
     return dv2 != boost::none ? ho
                               : boost::graph_traits<TriangleMesh>::null_halfedge();

  if(get(vcm, target(h, tmesh)))
     return dv1 != boost::none ? h
                               : boost::graph_traits<TriangleMesh>::null_halfedge();

  if(dv1 != boost::none)
  {
    if(dv2 != boost::none)
      return (*dv1 < *dv2) ? h : ho;

    return h;
  }

  if(dv2 != boost::none)
    return ho;

  return boost::graph_traits<TriangleMesh>::null_halfedge();
}

// adapted from triangulate_faces
template <typename TriangleMesh, typename VPM, typename Traits>
bool should_flip(typename boost::graph_traits<TriangleMesh>::edge_descriptor e,
                 const TriangleMesh& tmesh,
                 const VPM& vpm,
                 const Traits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_traits<VPM>::reference                 Point_ref;
  typedef typename Traits::Vector_3                                       Vector_3;

  CGAL_precondition(!is_border(e, tmesh));

  halfedge_descriptor h = halfedge(e, tmesh);

  Point_ref p0 = get(vpm, target(h, tmesh));
  Point_ref p1 = get(vpm, target(next(h, tmesh), tmesh));
  Point_ref p2 = get(vpm, source(h, tmesh));
  Point_ref p3 = get(vpm, target(next(opposite(h, tmesh), tmesh), tmesh));

  /* Chooses the diagonal that will split the quad in two triangles that maximize
   * the scalar product of of the un-normalized normals of the two triangles.
   * The lengths of the un-normalized normals (computed using cross-products of two vectors)
   *  are proportional to the area of the triangles.
   * Maximize the scalar product of the two normals will avoid skinny triangles,
   * and will also taken into account the cosine of the angle between the two normals.
   * In particular, if the two triangles are oriented in different directions,
   * the scalar product will be negative.
   */

//  CGAL::cross_product(p2-p1, p3-p2) * CGAL::cross_product(p0-p3, p1-p0);
//  CGAL::cross_product(p1-p0, p1-p2) * CGAL::cross_product(p3-p2, p3-p0);

  const Vector_3 v01 = gt.construct_vector_3_object()(p0, p1);
  const Vector_3 v12 = gt.construct_vector_3_object()(p1, p2);
  const Vector_3 v23 = gt.construct_vector_3_object()(p2, p3);
  const Vector_3 v30 = gt.construct_vector_3_object()(p3, p0);

  const double p1p3 = gt.compute_scalar_product_3_object()(
                        gt.construct_cross_product_vector_3_object()(v12, v23),
                        gt.construct_cross_product_vector_3_object()(v30, v01));

  const Vector_3 v21 = gt.construct_opposite_vector_3_object()(v12);
  const Vector_3 v03 = gt.construct_opposite_vector_3_object()(v30);

  const double p0p2 = gt.compute_scalar_product_3_object()(
                        gt.construct_cross_product_vector_3_object()(v01, v21),
                        gt.construct_cross_product_vector_3_object()(v23, v03));

  return p0p2 <= p1p3;
}

} // namespace internal

namespace experimental {

// @todo check what to use as priority queue with removable elements, set might not be optimal
template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_almost_degenerate_faces(const FaceRange& face_range,
                                    TriangleMesh& tmesh,
                                    const double cap_threshold,
                                    const double needle_threshold,
                                    const double collapse_length_threshold,
                                    const NamedParameters& np)
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;

  typedef Constant_property_map<vertex_descriptor, bool>                        Default_VCM;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                       NamedParameters,
                                                       Default_VCM>::type       VCM;
  VCM vcm_np = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained), Default_VCM(false));

  typedef Constant_property_map<edge_descriptor, bool>                          Default_ECM;
  typedef typename internal_np::Lookup_named_param_def<internal_np::edge_is_constrained_t,
                                                       NamedParameters,
                                                       Default_ECM>::type       ECM;
  ECM ecm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained), Default_ECM(false));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type           Traits;
  Traits gt = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  // Vertex property map that combines the VCM and the fact that extremities of a constrained edge should be constrained
  typedef CGAL::dynamic_vertex_property_t<bool>                                 Vertex_property_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_property_tag>::type DVCM;
  DVCM vcm = get(Vertex_property_tag(), tmesh);

  for(face_descriptor f : face_range)
  {
    if(f == boost::graph_traits<TriangleMesh>::null_face())
      continue;

    for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
    {
      if(get(ecm, edge(h, tmesh)))
      {
        put(vcm, source(h, tmesh), true);
        put(vcm, target(h, tmesh), true);
      }
      else if(get(vcm_np, target(h, tmesh)))
      {
        put(vcm, target(h, tmesh), true);
      }
    }
  }

  // Start the process of removing bad elements
  std::set<edge_descriptor> edges_to_collapse;
  std::set<edge_descriptor> edges_to_flip;

  // @todo could probably do something a bit better by looping edges, consider the incident faces
  // f1 / f2 and look at f1 if f1<f2, and the edge is smaller than the two other edges...
  for(face_descriptor f : face_range)
    internal::collect_badly_shaped_triangles(f, tmesh, vpm, ecm, gt,
                                             cap_threshold, needle_threshold, collapse_length_threshold,
                                             edges_to_collapse, edges_to_flip);

#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
  int iter = 0;
#endif

  for(;;)
  {
    bool something_was_done = false;

#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
    std::cout << edges_to_collapse.size() << " needles and " << edges_to_flip.size() << " caps" << std::endl;
    std::ostringstream oss;
    oss << "degen_cleaning_iter_" << iter++ << ".off";
    std::ofstream out(oss.str().c_str());
    out << std::setprecision(17);
    out << tmesh;
    out.close();
#endif

    if(edges_to_collapse.empty() && edges_to_flip.empty())
      return true;

    /// @todo maybe using a priority queue handling the more almost degenerate elements should be used
    std::set<edge_descriptor> next_edges_to_collapse;
    std::set<edge_descriptor> next_edges_to_flip;

    // treat needles
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
    int kk=0;
    std::ofstream(std::string("tmp/n-00000.off")) << tmesh;
#endif
    while(!edges_to_collapse.empty())
    {
      edge_descriptor e = *edges_to_collapse.begin();
      edges_to_collapse.erase(edges_to_collapse.begin());

      CGAL_assertion(!get(ecm, e));

      if(get(vcm, source(e, tmesh)) && get(vcm, target(e, tmesh)))
        continue;

#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "  treat needle: " << e << " (" << tmesh.point(source (e, tmesh))
                                        << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
#endif
      if(CGAL::Euler::does_satisfy_link_condition(e, tmesh))
      {
        // the following edges are removed by the collapse
        halfedge_descriptor h = halfedge(e, tmesh);
        CGAL_assertion(!is_border(h, tmesh)); // because extracted from a face

        std::array<halfedge_descriptor, 2> nc =
          internal::is_badly_shaped(face(h, tmesh), tmesh, vpm, ecm, gt,
                                    cap_threshold, needle_threshold, collapse_length_threshold);

        if(nc[0] != h)
        {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
          std::cerr << "Warning: Needle criteria no longer verified " << tmesh.point(source(e, tmesh)) << " "
                                                                      << tmesh.point(target(e, tmesh)) << std::endl;
#endif
          // the opposite edge might also have been inserted in the set and might still be a needle
          h = opposite(h, tmesh);
          if(is_border(h, tmesh))
            continue;

          nc = internal::is_badly_shaped(face(h, tmesh), tmesh, vpm, ecm, gt,
                                         cap_threshold, needle_threshold,
                                         collapse_length_threshold);
          if(nc[0] != h)
            continue;
        }

        for(int i=0; i<2; ++i)
        {
          if(!is_border(h, tmesh))
          {
            edge_descriptor pe = edge(prev(h, tmesh), tmesh);
            edges_to_flip.erase(pe);
            next_edges_to_collapse.erase(pe);
            edges_to_collapse.erase(pe);
          }

          h = opposite(h, tmesh);
        }

        // pick the orientation of edge to keep the vertex minimizing the volume variation
        h = internal::get_best_edge_orientation(e, tmesh, vpm, vcm, gt);

        if(h == boost::graph_traits<TriangleMesh>::null_halfedge())
        {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
            std::cerr << "Warning: geometrically invalid edge collapse! "
                      << tmesh.point(source(e, tmesh)) << " "
                      << tmesh.point(target(e, tmesh)) << std::endl;
#endif
          next_edges_to_collapse.insert(e);
          continue;
        }

        edges_to_flip.erase(e);
        next_edges_to_collapse.erase(e); // for edges added in faces incident to a vertex kept after a collapse
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
        std::cerr << "  " << kk << " -- Collapsing " << tmesh.point(source(h, tmesh)) << "  "
                                                     << tmesh.point(target(h, tmesh)) << std::endl;
#endif
        // moving to the midpoint is not a good idea. On a circle for example you might endpoint with
        // a bad geometry because you iteratively move one point
        // auto mp = midpoint(tmesh.point(source(h, tmesh)), tmesh.point(target(h, tmesh)));

        vertex_descriptor v = Euler::collapse_edge(edge(h, tmesh), tmesh);

        //tmesh.point(v) = mp;
        // examine all faces incident to the vertex kept
        for(halfedge_descriptor hv : halfedges_around_target(v, tmesh))
        {
          if(!is_border(hv, tmesh))
          {
            internal::collect_badly_shaped_triangles(face(hv, tmesh), tmesh, vpm, ecm, gt,
                                                     cap_threshold, needle_threshold, collapse_length_threshold,
                                                     edges_to_collapse, edges_to_flip);
          }
        }

#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
        std::string nb = std::to_string(++kk);
        if(kk<10) nb = std::string("0")+nb;
        if(kk<100) nb = std::string("0")+nb;
        if(kk<1000) nb = std::string("0")+nb;
        if(kk<10000) nb = std::string("0")+nb;
        std::ofstream(std::string("tmp/n-")+nb+std::string(".off")) << tmesh;
#endif
        something_was_done = true;
      }
      else
      {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cerr << "Warning: uncollapsable edge! " << tmesh.point(source(e, tmesh)) << " "
                                                     << tmesh.point(target(e, tmesh)) << std::endl;
#endif
        next_edges_to_collapse.insert(e);
      }
    }

    // treat caps
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
    kk=0;
    std::ofstream(std::string("tmp/c-000.off")) << tmesh;
#endif
    while(!edges_to_flip.empty())
    {
      edge_descriptor e = *edges_to_flip.begin();
      edges_to_flip.erase(edges_to_flip.begin());

      CGAL_assertion(!get(ecm, e));

      if(get(vcm, source(e, tmesh)) && get(vcm, target(e, tmesh)))
        continue;

#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "treat cap: " << e << " (" << tmesh.point(source(e, tmesh))
                                   << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
#endif

      halfedge_descriptor h = halfedge(e, tmesh);
      std::array<halfedge_descriptor,2> nc = internal::is_badly_shaped(face(h, tmesh), tmesh, vpm, ecm, gt,
                                                                       cap_threshold, needle_threshold,
                                                                       collapse_length_threshold);
      // First check the triangle is still a cap
      if(nc[1] != h)
      {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cerr << "Warning: Cap criteria no longer verified " << tmesh.point(source(e, tmesh)) << " --- "
                                                                 << tmesh.point(target(e, tmesh)) << std::endl;
#endif
        // the opposite edge might also have been inserted in the set and might still be a cap
        h = opposite(h, tmesh);
        if(is_border(h, tmesh))
          continue;

        nc = internal::is_badly_shaped(face(h, tmesh), tmesh, vpm, ecm, gt,
                                       cap_threshold, needle_threshold, collapse_length_threshold);
        if(nc[1] != h)
          continue;
      }

      // special case on the border
      if(is_border(opposite(h, tmesh), tmesh))
      {
        // remove the triangle
        edges_to_flip.erase(edge(prev(h, tmesh), tmesh));
        edges_to_flip.erase(edge(next(h, tmesh), tmesh));
        next_edges_to_collapse.erase(edge(prev(h, tmesh), tmesh));
        next_edges_to_collapse.erase(edge(next(h, tmesh), tmesh));
        Euler::remove_face(h, tmesh);
        something_was_done = true;
        continue;
      }

      // condition for the flip to be valid (the edge to be created does not already exist)
      if(!halfedge(target(next(h, tmesh), tmesh),
                   target(next(opposite(h, tmesh), tmesh), tmesh), tmesh).second)
      {

        if(!internal::should_flip(e, tmesh, vpm, gt))
        {
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
          std::cout << "Flipping prevented: not the best diagonal" << std::endl;
#endif
          next_edges_to_flip.insert(e);
          continue;
        }

#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cout << "Flipping" << std::endl;
#endif
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
        std::cerr << "step " << kk << "\n";
        std::cerr << "  Flipping " << tmesh.point(source(h, tmesh)) << "  "
                                                   << tmesh.point(target(h, tmesh)) << std::endl;
#endif
        Euler::flip_edge(h, tmesh);
        CGAL_assertion(edge(h, tmesh) == e);

        // handle face updates
        for(int i=0; i<2; ++i)
        {
          CGAL_assertion(!is_border(h, tmesh));
          std::array<halfedge_descriptor, 2> nc =
            internal::is_badly_shaped(face(h, tmesh), tmesh, vpm, ecm, gt,
                                      cap_threshold, needle_threshold, collapse_length_threshold);

          if(nc[1] != boost::graph_traits<TriangleMesh>::null_halfedge())
          {
            if(edge(nc[1], tmesh) != e)
              next_edges_to_flip.insert(edge(nc[1], tmesh));
          }
          else
          {
            if(nc[0] != boost::graph_traits<TriangleMesh>::null_halfedge())
            {
              next_edges_to_collapse.insert(edge(nc[0], tmesh));
            }
          }
          h = opposite(h, tmesh);
        }
        something_was_done = true;
      }
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      else
      {
        std::cerr << "Warning: unflippable edge! " << tmesh.point(source(h, tmesh)) << " --- "
                                                   << tmesh.point(target(h, tmesh)) << std::endl;
        next_edges_to_flip.insert(e);
      }
#endif
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
      std::string nb = std::to_string(++kk);
      if(kk<10) nb = std::string("0")+nb;
      if(kk<100) nb = std::string("0")+nb;
      if(kk<1000) nb = std::string("0")+nb;
      if(kk<10000) nb = std::string("0")+nb;
      std::ofstream(std::string("tmp/c-")+nb+std::string(".off")) << tmesh;
#endif
    }

    std::swap(edges_to_collapse, next_edges_to_collapse);
    std::swap(edges_to_flip, next_edges_to_flip);

    if(!something_was_done)
      return false;
  }

  return false;
}

template <typename FaceRange, typename TriangleMesh>
bool remove_almost_degenerate_faces(const FaceRange& face_range,
                                    TriangleMesh& tmesh,
                                    const double cap_threshold,
                                    const double needle_threshold,
                                    const double collapse_length_threshold)
{
  return remove_almost_degenerate_faces(face_range, tmesh,
                                        cap_threshold, needle_threshold, collapse_length_threshold,
                                        CGAL::parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_almost_degenerate_faces(TriangleMesh& tmesh,
                                    const double cap_threshold,
                                    const double needle_threshold,
                                    const double collapse_length_threshold,
                                    const CGAL_PMP_NP_CLASS& np)
{
  return remove_almost_degenerate_faces(faces(tmesh), tmesh, cap_threshold, needle_threshold,
                                        collapse_length_threshold, np);
}

template<class TriangleMesh>
bool remove_almost_degenerate_faces(TriangleMesh& tmesh,
                                    const double cap_threshold,
                                    const double needle_threshold,
                                    const double collapse_length_threshold)
{
  return remove_almost_degenerate_faces(faces(tmesh), tmesh,
                                        cap_threshold, needle_threshold, collapse_length_threshold,
                                        CGAL::parameters::all_default());
}

} // namespace experimental
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REMOVE_DEGENERACIES_H
