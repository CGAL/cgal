// Copyright (c) 2015-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_DEGENERACIES_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_DEGENERACIES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/boost/graph/selection.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/property_map.h>
#include <CGAL/Union_find.h>

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/IO/OFF_reader.h>
#endif

#include <boost/algorithm/minmax_element.hpp>
#include <boost/utility/enable_if.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

// First part of the file: remove_ALMOST_degenerate_faces (needles/caps)
// Second part of the file: remove_degenerate_edges/faces

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

  // @todo handle boundary edges

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
      // @todo use a predicate
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
boost::optional<typename Traits::FT>
get_collapse_volume(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                    const TriangleMesh& tmesh,
                    const VPM& vpm,
                    const Traits& gt)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;

  typedef typename Traits::FT                                                     FT;
  typedef typename boost::property_traits<VPM>::reference                         Point_ref;
  typedef typename Traits::Vector_3                                               Vector_3;

  const typename Traits::Point_3 origin(ORIGIN);

// @todo handle boundary edges

  h = opposite(h, tmesh); // Euler::collapse edge keeps the target and removes the source

  // source is kept, target is removed
  Point_ref kept = get(vpm, source(h, tmesh));
  Point_ref removed= get(vpm, target(h, tmesh));

  // init volume with incident triangles (reversed orientation
  FT delta_vol = volume(removed, kept, get(vpm, target(next(h, tmesh), tmesh)), origin) +
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
      // @todo use a predicate
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
  typedef typename Traits::FT                                             FT;

  halfedge_descriptor h = halfedge(e, tmesh), ho = opposite(h, tmesh);

  CGAL_assertion(!get(vcm, source(h, tmesh)) || !get(vcm, target(h, tmesh)));

  boost::optional<FT> dv1 = get_collapse_volume(h, tmesh, vpm, gt);
  boost::optional<FT> dv2 = get_collapse_volume(ho, tmesh, vpm, gt);

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

  typedef typename Traits::FT                                             FT;
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

  const FT p1p3 = gt.compute_scalar_product_3_object()(
                    gt.construct_cross_product_vector_3_object()(v12, v23),
                    gt.construct_cross_product_vector_3_object()(v30, v01));

  const Vector_3 v21 = gt.construct_opposite_vector_3_object()(v12);
  const Vector_3 v03 = gt.construct_opposite_vector_3_object()(v30);

  const FT p0p2 = gt.compute_scalar_product_3_object()(
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

  typedef Static_boolean_property_map<vertex_descriptor, false>                 Default_VCM;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                       NamedParameters,
                                                       Default_VCM>::type       VCM;
  VCM vcm_np = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained), Default_VCM());

  typedef Static_boolean_property_map<edge_descriptor, false>                   Default_ECM;
  typedef typename internal_np::Lookup_named_param_def<internal_np::edge_is_constrained_t,
                                                       NamedParameters,
                                                       Default_ECM>::type       ECM;
  ECM ecm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained), Default_ECM());

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

    // @todo maybe using a priority queue handling the more almost degenerate elements should be used
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

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/// Remove degenerate_edges/faces

namespace internal {

template <typename HalfedgeGraph, typename VertexPointMap, typename Traits>
struct Less_vertex_point
{
  typedef typename boost::graph_traits<HalfedgeGraph>::vertex_descriptor    vertex_descriptor;

  Less_vertex_point(const Traits& traits, const VertexPointMap& vpmap)
    : m_traits(traits), m_vpmap(vpmap)
  {}

  bool operator()(vertex_descriptor v1, vertex_descriptor v2) const {
    return m_traits.less_xyz_3_object()(get(m_vpmap, v1), get(m_vpmap, v2));
  }

  const Traits& m_traits;
  const VertexPointMap& m_vpmap;
};

} // end namespace internal

// this function removes a border edge even if it does not satisfy the link condition.
// null_vertex() is returned if the removal changes the topology of the input
template <typename TriangleMesh, typename EdgeSet, typename FaceSet>
typename boost::graph_traits<TriangleMesh>::vertex_descriptor
remove_a_border_edge(typename boost::graph_traits<TriangleMesh>::edge_descriptor ed,
                     TriangleMesh& tm,
                     EdgeSet& input_range,
                     EdgeSet& edge_set,
                     FaceSet& face_set)
{
  typedef boost::graph_traits<TriangleMesh>                             GT;
  typedef typename GT::vertex_descriptor                                vertex_descriptor;
  typedef typename GT::halfedge_descriptor                              halfedge_descriptor;
  typedef typename GT::edge_descriptor                                  edge_descriptor;
  typedef typename GT::face_descriptor                                  face_descriptor;

  halfedge_descriptor h = halfedge(ed, tm);

  if(is_border(h, tm))
    h = opposite(h, tm);

  halfedge_descriptor opp_h = opposite(h, tm);
  CGAL_assertion(is_border(opp_h, tm));
  CGAL_assertion(!is_border(h, tm));

  CGAL_assertion(next(next(opp_h, tm), tm) != opp_h); // not working for a hole made of 2 edges
  CGAL_assertion(next(next(next(opp_h, tm), tm), tm) != opp_h); // not working for a hole make of 3 edges

  if(CGAL::Euler::does_satisfy_link_condition(edge(h, tm), tm))
  {
    edge_set.erase(ed);
    input_range.erase(ed);
    halfedge_descriptor h = halfedge(ed, tm);
    if(is_border(h, tm))
      h = opposite(h, tm);

    const edge_descriptor prev_e = edge(prev(h, tm), tm);
    edge_set.erase(prev_e);
    input_range.erase(prev_e);
    face_set.erase(face(h, tm));

    return CGAL::Euler::collapse_edge(ed, tm);
  }

  // collect edges that have one vertex in the link of
  // the vertices of h and one of the vertex of h as other vertex
  std::set<edge_descriptor> common_incident_edges;
  for(halfedge_descriptor hos : halfedges_around_source(h, tm))
  {
    for(halfedge_descriptor hot : halfedges_around_target(h, tm))
    {
      if(target(hos, tm) == source(hot, tm))
      {
        common_incident_edges.insert(edge(hot, tm));
        common_incident_edges.insert(edge(hos, tm));
      }
    }
  }

  CGAL_assertion(common_incident_edges.size() >= 2);

  // in the following loop, we visit define a connected component of
  // faces bounded by edges in common_incident_edges and h. We look
  // for the maximal one. This set of faces is the one that will
  // disappear while collapsing ed
  std::set<face_descriptor> marked_faces;

  std::vector<halfedge_descriptor> queue;
  queue.push_back(opposite(next(h, tm), tm));
  queue.push_back(opposite(prev(h, tm), tm));
  marked_faces.insert(face(h, tm));

  do
  {
    std::vector<halfedge_descriptor> boundary;
    while(!queue.empty())
    {
      halfedge_descriptor back=queue.back();
      queue.pop_back();
      face_descriptor fback=face(back, tm);
      if(common_incident_edges.count(edge(back, tm)))
      {
        boundary.push_back(back);
        continue;
      }

      if(fback==GT::null_face() || !marked_faces.insert(fback).second)
        continue;

      queue.push_back(opposite(next(back, tm), tm));
      if(is_border(queue.back(), tm))
      {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "Boundary reached during exploration, the region to be removed is not a topological disk, not handled for now.\n";
#endif
        return GT::null_vertex();
      }

      queue.push_back(opposite(prev(back, tm), tm));
      if(is_border(queue.back(), tm))
      {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "Boundary reached during exploration, the region to be removed is not a topological disk, not handled for now.\n";
#endif
        return GT::null_vertex();
      }
    }

    CGAL_assertion(boundary.size() == 2);
    common_incident_edges.erase(edge(boundary[0], tm));
    common_incident_edges.erase(edge(boundary[1], tm));
    if(!is_border(boundary[0], tm) || common_incident_edges.empty())
      queue.push_back(boundary[0]);
    if(!is_border(boundary[1], tm) || common_incident_edges.empty())
      queue.push_back(boundary[1]);
  }
  while(!common_incident_edges.empty());

  // hk1 and hk2 are bounding the region that will be removed.
  // The edge of hk2 will be removed and hk2 will be replaced
  // by the opposite edge of hk1
  halfedge_descriptor hk1 = queue.front();
  halfedge_descriptor hk2 = queue.back();
  if(target(hk1, tm)!=source(hk2, tm))
    std::swap(hk1, hk2);

  CGAL_assertion(target(hk1, tm) == source(hk2, tm));
  CGAL_assertion(source(hk1, tm) == source(h, tm));
  CGAL_assertion(target(hk2, tm) == target(h, tm));

  CGAL_assertion(is_valid_polygon_mesh(tm));
  if(!is_selection_a_topological_disk(marked_faces, tm))
  {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
    std::cout << "The region to be removed is not a topological disk, not handled for now.\n";
#endif
    return GT::null_vertex();
  }

  if(is_border(hk1, tm) && is_border(hk2, tm))
  {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
    std::cout << "The region to be removed is an isolated region, not handled for now.\n";
#endif
    return GT::null_vertex();
  }

  // collect vertices and edges to remove and do remove faces
  std::set<edge_descriptor> edges_to_remove;
  std::set<vertex_descriptor> vertices_to_remove;
  for(face_descriptor fd : marked_faces)
  {
    halfedge_descriptor hd = halfedge(fd, tm);
    for(int i=0; i<3; ++i)
    {
      edges_to_remove.insert(edge(hd, tm));
      vertices_to_remove.insert(target(hd, tm));
      hd = next(hd, tm);
    }
  }

  vertex_descriptor vkept = source(hk1, tm);

  //back-up next, prev halfedge pointers to be restored after removal
  halfedge_descriptor hp = prev(opp_h, tm);
  halfedge_descriptor hn = next(opp_h, tm);
  halfedge_descriptor hk1_opp_next = next(hk2, tm);
  halfedge_descriptor hk1_opp_prev = prev(hk2, tm);
  face_descriptor hk1_opp_face = face(hk2, tm);

  // we will remove the target of hk2, update vertex pointers
  for(halfedge_descriptor hot : halfedges_around_target(hk2, tm))
    set_target(hot, vkept, tm);

  // update halfedge pointers since hk2 will be removed
  set_halfedge(vkept, opposite(hk1, tm), tm);
  set_halfedge(target(hk1, tm), hk1, tm);

  // do not remove hk1 and its vertices
  vertices_to_remove.erase(vkept);
  vertices_to_remove.erase(target(hk1, tm));
  edges_to_remove.erase(edge(hk1, tm));

  bool hk2_equals_hp = (hk2 == hp);
  CGAL_assertion(is_border(hk2, tm) == hk2_equals_hp);

  /*
  - case hk2!=hp

         /\      /
     hk1/  \hk2 /
       /    \  /
  ____/______\/____
  hn   h_opp   hp

  - case hk2==hp

         /\
     hk1/  \hk2 == hp
       /    \
  ____/______\
  hn   h_opp
  */

  // remove vertices
  for(vertex_descriptor vd : vertices_to_remove)
    remove_vertex(vd, tm);

  // remove edges
  for(edge_descriptor ed : edges_to_remove)
  {
    edge_set.erase(ed);
    input_range.erase(ed);
    remove_edge(ed, tm);
  }

  // remove faces
  for(face_descriptor fd : marked_faces)
  {
    face_set.erase(fd);
    remove_face(fd, tm);
  }

  // now update pointers
  set_face(opposite(hk1, tm), hk1_opp_face, tm);
  if(!hk2_equals_hp)
  {
    set_next(hp, hn, tm);
    set_next(opposite(hk1, tm), hk1_opp_next, tm);
    set_next(hk1_opp_prev, opposite(hk1, tm), tm);
    set_halfedge(hk1_opp_face, opposite(hk1, tm), tm);
  }
  else
  {
    set_next(hk1_opp_prev, opposite(hk1, tm), tm);
    set_next(opposite(hk1, tm), hn, tm);
  }

  CGAL_assertion(is_valid_polygon_mesh(tm));
  return vkept;
}

template <typename TriangleMesh>
typename boost::graph_traits<TriangleMesh>::vertex_descriptor
remove_a_border_edge(typename boost::graph_traits<TriangleMesh>::edge_descriptor ed,
                     TriangleMesh& tm)
{
  std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor> input_range;
  std::set<typename boost::graph_traits<TriangleMesh>::edge_descriptor> edge_set;
  std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> face_set;

  return remove_a_border_edge(ed, tm, input_range, edge_set, face_set);
}

template <typename EdgeRange, typename TriangleMesh, typename NamedParameters, typename FaceSet>
bool remove_degenerate_edges(const EdgeRange& edge_range,
                             TriangleMesh& tmesh,
                             FaceSet& face_set,
                             const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tmesh));
  CGAL_assertion(CGAL::is_valid_polygon_mesh(tmesh));

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef TriangleMesh                                                            TM;
  typedef typename boost::graph_traits<TriangleMesh>                              GT;
  typedef typename GT::vertex_descriptor                                          vertex_descriptor;
  typedef typename GT::halfedge_descriptor                                        halfedge_descriptor;
  typedef typename GT::edge_descriptor                                            edge_descriptor;
  typedef typename GT::face_descriptor                                            face_descriptor;

  typedef typename GetVertexPointMap<TM, NamedParameters>::type                   VertexPointMap;
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_property_map(vertex_point, tmesh));

  typedef typename GetGeomTraits<TM, NamedParameters>::type                       Traits;

  std::size_t nb_deg_faces = 0;
  bool all_removed = false;
  bool some_removed = true;
  bool preserve_genus = choose_parameter(get_parameter(np, internal_np::preserve_genus), true);

  // The input edge range needs to be kept up-to-date
  std::set<edge_descriptor> local_edge_range(std::begin(edge_range), std::end(edge_range));

  // collect edges of length 0
  while(some_removed && !all_removed)
  {
    some_removed = false;
    all_removed = true;
    std::set<edge_descriptor> degenerate_edges_to_remove;
    degenerate_edges(local_edge_range, tmesh, std::inserter(degenerate_edges_to_remove,
                                                            degenerate_edges_to_remove.end()));

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
    std::cout << "Found " << degenerate_edges_to_remove.size() << " null edges.\n";
#endif

    // first try to remove all collapsable edges
    typename std::set<edge_descriptor>::iterator it = degenerate_edges_to_remove.begin();
    while(it != degenerate_edges_to_remove.end())
    {
      edge_descriptor e = *it;
      if(CGAL::Euler::does_satisfy_link_condition(e, tmesh))
      {
        const halfedge_descriptor h = halfedge(e, tmesh);
        local_edge_range.erase(*it);
        degenerate_edges_to_remove.erase(it);

        // remove edges that could also be set for removal
        if(face(h, tmesh) != GT::null_face())
        {
          ++nb_deg_faces;
          const edge_descriptor prev_e = edge(prev(h, tmesh), tmesh);
          degenerate_edges_to_remove.erase(prev_e);
          local_edge_range.erase(prev_e);
          face_set.erase(face(h, tmesh));
        }

        if(face(opposite(h, tmesh), tmesh) != GT::null_face())
        {
          ++nb_deg_faces;
          const edge_descriptor prev_opp_e = edge(prev(opposite(h, tmesh), tmesh), tmesh);
          degenerate_edges_to_remove.erase(prev_opp_e);
          local_edge_range.erase(prev_opp_e);
          face_set.erase(face(opposite(h, tmesh), tmesh));
        }

        //now remove the edge
        CGAL::Euler::collapse_edge(e, tmesh);

        // some_removed is not updated on purpose because if nothing
        //  happens below then nothing can be done
        it = degenerate_edges_to_remove.begin();
      }
      else // link condition not satisfied
      {
        ++it;
      }
    }

    CGAL_assertion(is_valid_polygon_mesh(tmesh));
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
    std::cout << "Remaining " << degenerate_edges_to_remove.size() << " null edges to be handled.\n";
#endif

    while(!degenerate_edges_to_remove.empty())
    {
      auto eb = degenerate_edges_to_remove.begin();
      const edge_descriptor ed = *eb;
      degenerate_edges_to_remove.erase(eb);
      local_edge_range.erase(ed);

      const halfedge_descriptor h = halfedge(ed, tmesh);

      if(CGAL::Euler::does_satisfy_link_condition(ed, tmesh))
      {
        // remove edges that could also be set for removal
        if(face(h, tmesh) != GT::null_face())
        {
          ++nb_deg_faces;
          const edge_descriptor prev_e = edge(prev(h, tmesh), tmesh);
          degenerate_edges_to_remove.erase(prev_e);
          local_edge_range.erase(prev_e);
          face_set.erase(face(h, tmesh));
        }

        if(face(opposite(h, tmesh), tmesh)!=GT::null_face())
        {
          ++nb_deg_faces;
          const edge_descriptor prev_opp_e = edge(prev(opposite(h, tmesh), tmesh), tmesh);
          degenerate_edges_to_remove.erase(prev_opp_e);
          local_edge_range.erase(prev_opp_e);
          face_set.erase(face(opposite(h, tmesh), tmesh));
        }

        //now remove the edge
        CGAL::Euler::collapse_edge(ed, tmesh);
        some_removed = true;
      }
      else // link condition not satisfied
      {
        // handle the case when the edge is incident to a triangle hole
        // we first fill the hole and try again
        if(is_border(ed, tmesh))
        {
          halfedge_descriptor hd = halfedge(ed, tmesh);
          if(!is_border(hd, tmesh))
            hd = opposite(hd, tmesh);

          if(is_triangle(hd, tmesh))
          {
            if(!preserve_genus)
            {
              Euler::fill_hole(hd, tmesh);
              degenerate_edges_to_remove.insert(ed); // reinsert the edge for future processing
              local_edge_range.insert(ed);
            }
            else
            {
              all_removed=false;
            }

            continue;
          }

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
          std::cout << "Calling remove_a_border_edge\n";
#endif

          vertex_descriptor vd = remove_a_border_edge(ed, tmesh, local_edge_range,
                                                      degenerate_edges_to_remove, face_set);
          if(vd == GT::null_vertex())
          {
            // @todo: if some border edges are later removed, the edge might be processable later
            // for example if it belongs to  boundary cycle of edges where the number of non-degenerate
            // edges is 2. That's what happen with fused_vertices.off in the testsuite where the edges
            // are not processed the same way with Polyhedron and Surface_mesh. In the case of Polyhedron
            // more degenerate edges could be removed.
            all_removed = false;
          }
          else
            some_removed = true;

          continue;
        }
        else
        {
          halfedge_descriptor hd = halfedge(ed, tmesh);
          // if both vertices are boundary vertices we can't do anything
          bool impossible = false;
          for(halfedge_descriptor h : halfedges_around_target(hd, tmesh))
          {
            if(is_border(h, tmesh))
            {
              impossible = true;
              break;
            }
          }

          if(impossible)
          {
            impossible = false;
            for(halfedge_descriptor h : halfedges_around_source(hd, tmesh))
            {
              if(is_border(h, tmesh))
              {
                impossible = true;
                break;
              }
            }

            if(impossible)
            {
              all_removed = false;
              continue;
            }
          }
        }

        // When the edge does not satisfy the link condition, it means that it cannot be
        // collapsed as is. In the following if there is a topological issue
        // with contracting the edge (component or geometric feature that disappears),
        // nothing is done.
        // We start by marking the faces that are incident to an edge endpoint.
        // If the set of marked faces is a topologically disk, then we simply remove all the simplicies
        // inside the disk and star the hole with the edge vertex kept.
        // If the set of marked faces is not a topological disk, it has some non-manifold vertices
        // on its boundary. We need to mark additional faces to make it a topological disk.
        // We can then apply the star hole procedure.
        // Right now we additionally mark the smallest connected components of non-marked faces
        // (using the number of faces)

        //backup central point
        typename Traits::Point_3 pt = get(vpmap, source(ed, tmesh));

        // mark faces of the link of each endpoints of the edge which collapse is not topologically valid
        std::set<face_descriptor> marked_faces;

        //   first endpoint
        for(halfedge_descriptor hd : CGAL::halfedges_around_target(halfedge(ed, tmesh), tmesh))
          if(!is_border(hd, tmesh))
            marked_faces.insert(face(hd, tmesh));

        //   second endpoint
        for(halfedge_descriptor hd : CGAL::halfedges_around_target(opposite(halfedge(ed, tmesh), tmesh), tmesh))
          if(!is_border(hd, tmesh))
            marked_faces.insert(face(hd, tmesh));

        // extract the halfedges on the boundary of the marked region
        std::vector<halfedge_descriptor> border;
        for(face_descriptor fd : marked_faces)
        {
          for(halfedge_descriptor hd : CGAL::halfedges_around_face(halfedge(fd, tmesh), tmesh))
          {
            halfedge_descriptor hd_opp = opposite(hd, tmesh);
            if(is_border(hd_opp, tmesh) || marked_faces.count(face(hd_opp, tmesh)) == 0)
            {
              border.push_back(hd);
            }
          }
        }

        if(border.empty())
        {
          // a whole connected component (without boundary) got selected and will disappear (not handled for now)
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
          std::cout << "Trying to remove a whole connected component, not handled yet\n";
#endif
          all_removed = false;
          continue;
        }

        // define cc of border halfedges: two halfedges are in the same cc
        // if they are on the border of the cc of non-marked faces.
        typedef CGAL::Union_find<halfedge_descriptor>                     UF_ds;
        UF_ds uf;
        std::map<halfedge_descriptor, typename UF_ds::handle> handles;

        // one cc per border halfedge
        for(halfedge_descriptor hd : border)
          handles.insert(std::make_pair(hd, uf.make_set(hd)));

        // join cc's
        for(halfedge_descriptor hd : border)
        {
          CGAL_assertion(marked_faces.count(face(hd, tmesh)) > 0);
          CGAL_assertion(marked_faces.count(face(opposite(hd, tmesh), tmesh)) == 0);
          halfedge_descriptor candidate = hd;

          do
          {
            candidate = prev(opposite(candidate, tmesh), tmesh);
          }
          while(!marked_faces.count(face(opposite(candidate, tmesh), tmesh)));

          uf.unify_sets(handles[hd], handles[opposite(candidate, tmesh)]);
        }

        std::size_t nb_cc = uf.number_of_sets();
        if(nb_cc != 1)
        {
          // if more than one connected component is found then the patch
          // made of marked faces contains "non-manifold" vertices.
          // The smallest components need to be marked so that the patch
          // made of marked faces is a topological disk

          // we will explore in parallel the connected components and will stop
          // when all but one connected component have been entirely explored.
          // We add one face at a time for each cc in order to not explore a
          // potentially very large cc.
          std::vector<std::vector<halfedge_descriptor> > stacks_per_cc(nb_cc);
          std::vector<std::set<face_descriptor> > faces_per_cc(nb_cc);
          std::vector<bool> exploration_finished(nb_cc, false);

          // init the stacks of halfedges using the cc of the boundary
          std::size_t index = 0;
          std::map<halfedge_descriptor, std::size_t > ccs;

          typedef std::pair<const halfedge_descriptor, typename UF_ds::handle> Pair_type;
          for(const Pair_type& p : handles)
          {
            halfedge_descriptor opp_hedge = opposite(p.first, tmesh);
            if(is_border(opp_hedge, tmesh))
              continue; // nothing to do on the boundary

            typedef typename std::map<halfedge_descriptor, std::size_t >::iterator Map_it;
            std::pair<Map_it, bool> insert_res = ccs.insert(std::make_pair(*uf.find(p.second), index));

            if(insert_res.second)
              ++index;

            stacks_per_cc[insert_res.first->second].push_back(prev(opp_hedge, tmesh));
            stacks_per_cc[insert_res.first->second].push_back(next(opp_hedge, tmesh));
            faces_per_cc[insert_res.first->second].insert(face(opp_hedge, tmesh));
          }

          if(index != nb_cc)
          {
            // most probably, one cc is a cycle of border edges
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
            std::cout << "Trying to remove a component with a cycle of halfedges (nested hole or whole component), not handled yet.\n";
#endif
            all_removed = false;
            continue;
          }

          std::size_t nb_ccs_to_be_explored = nb_cc;
          index = 0;
          //explore the cc's
          do
          {
            // try to extract one more face for a given cc
            do
            {
              CGAL_assertion(!exploration_finished[index]);
              CGAL_assertion(!stacks_per_cc.empty());
              CGAL_assertion(!stacks_per_cc[index].empty());

              halfedge_descriptor hd = stacks_per_cc[index].back();
              stacks_per_cc[index].pop_back();
              hd = opposite(hd, tmesh);

              if(!is_border(hd, tmesh) && !marked_faces.count(face(hd, tmesh)))
              {
                if(faces_per_cc[index].insert(face(hd, tmesh)).second)
                {
                  stacks_per_cc[index].push_back(next(hd, tmesh));
                  stacks_per_cc[index].push_back(prev(hd, tmesh));
                  break;
                }
              }

              if(stacks_per_cc[index].empty())
                break;
            }
            while(true);

            // the exploration of a cc is finished when its stack is empty
            exploration_finished[index]=stacks_per_cc[index].empty();
            if(exploration_finished[index])
              --nb_ccs_to_be_explored;

            if(nb_ccs_to_be_explored==1)
              break;

            while(exploration_finished[(++index)%nb_cc]);

            index = index%nb_cc;
          }
          while(true);

          // @todo use the area criteria? this means maybe continue exploration of larger cc
          // mark faces of completetly explored cc
          for(index=0; index<nb_cc; ++index)
          {
            if(exploration_finished[index])
            {
              for(face_descriptor fd : faces_per_cc[index])
                marked_faces.insert(fd);
            }
          }
        }

        // make sure the selection is a topological disk (otherwise we need another treatment)
        if(!is_selection_a_topological_disk(marked_faces, tmesh))
        {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
          std::cout << "Trying to handle a non-topological disk, do nothing\n";
#endif
          all_removed = false;
          continue;
        }

        some_removed = true;

        // collect simplices to be removed
        std::set<vertex_descriptor> vertices_to_keep;
        std::set<halfedge_descriptor> halfedges_to_keep;
        for(halfedge_descriptor hd : border)
        {
          if(!marked_faces.count(face(opposite(hd, tmesh), tmesh)))
          {
            halfedges_to_keep.insert(hd);
            vertices_to_keep.insert(target(hd, tmesh));
          }
        }

        // backup next,prev relationships to set after patch removal
        std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> > next_prev_halfedge_pairs;
        halfedge_descriptor first_border_hd = *(halfedges_to_keep.begin());
        halfedge_descriptor current_border_hd = first_border_hd;
        do
        {
          halfedge_descriptor prev_border_hd = current_border_hd;
          current_border_hd = next(current_border_hd, tmesh);

          while(marked_faces.count(face(opposite(current_border_hd, tmesh), tmesh)))
            current_border_hd=next(opposite(current_border_hd, tmesh), tmesh);

          next_prev_halfedge_pairs.push_back(std::make_pair(prev_border_hd, current_border_hd));
        }
        while(current_border_hd!=first_border_hd);

        // collect vertices and edges to remove and do remove faces
        std::set<edge_descriptor> edges_to_remove;
        std::set<vertex_descriptor> vertices_to_remove;
        for(face_descriptor fd : marked_faces)
        {
          halfedge_descriptor hd=halfedge(fd, tmesh);
          for(int i=0; i<3; ++i)
          {
            if(!halfedges_to_keep.count(hd))
              edges_to_remove.insert(edge(hd, tmesh));

            if(!vertices_to_keep.count(target(hd, tmesh)))
              vertices_to_remove.insert(target(hd, tmesh));

            hd = next(hd, tmesh);
          }

          remove_face(fd, tmesh);
          face_set.erase(fd);
        }

        // remove vertices
        for(vertex_descriptor vd : vertices_to_remove)
          remove_vertex(vd, tmesh);

        // remove edges
        for(edge_descriptor ed : edges_to_remove)
        {
          degenerate_edges_to_remove.erase(ed);
          local_edge_range.erase(ed);
          remove_edge(ed, tmesh);
        }

        // add a new face, set all border edges pointing to it
        // and update halfedge vertex of patch boundary vertices
        face_descriptor new_face = add_face(tmesh);
        typedef std::pair<halfedge_descriptor, halfedge_descriptor> Pair_type;
        for(const Pair_type& p : next_prev_halfedge_pairs)
        {
          set_face(p.first, new_face, tmesh);
          set_next(p.first, p.second, tmesh);
          set_halfedge(target(p.first, tmesh), p.first, tmesh);
        }

        set_halfedge(new_face, first_border_hd, tmesh);

        // triangulate the new face and update the coordinate of the central vertex
        halfedge_descriptor new_hd=Euler::add_center_vertex(first_border_hd, tmesh);
        put(vpmap, target(new_hd, tmesh), pt);

        for(halfedge_descriptor hd : halfedges_around_target(new_hd, tmesh))
        {
          const edge_descriptor inc_e = edge(hd, tmesh);
          if(is_degenerate_edge(inc_e, tmesh, np))
          {
            degenerate_edges_to_remove.insert(inc_e);
            local_edge_range.insert(inc_e);
          }

          if(face(hd, tmesh) != GT::null_face() && is_degenerate_triangle_face(face(hd, tmesh), tmesh))
            face_set.insert(face(hd, tmesh));
        }

        CGAL_assertion(is_valid_polygon_mesh(tmesh));
      }
    }
  }

  return all_removed;
}

template <typename EdgeRange, typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_degenerate_edges(const EdgeRange& edge_range,
                             TriangleMesh& tmesh,
                             const CGAL_PMP_NP_CLASS& np)
{
  std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> face_set;
  return remove_degenerate_edges(edge_range, tmesh, face_set, np);
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_degenerate_edges(TriangleMesh& tmesh,
                             const CGAL_PMP_NP_CLASS& np)
{
  std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> face_set;
  return remove_degenerate_edges(edges(tmesh), tmesh, face_set, np);
}

template <typename EdgeRange, typename TriangleMesh>
bool remove_degenerate_edges(const EdgeRange& edge_range,
                             TriangleMesh& tmesh)
{
  std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> face_set;
  return remove_degenerate_edges(edge_range, tmesh, face_set, parameters::all_default());
}

template <typename TriangleMesh>
bool remove_degenerate_edges(TriangleMesh& tmesh)
{
  std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> face_set;
  return remove_degenerate_edges(edges(tmesh), tmesh, face_set, parameters::all_default());
}

// \ingroup PMP_repairing_grp
// removes the degenerate faces from a triangulated surface mesh.
// A face is considered degenerate if two of its vertices share the same location,
// or more generally if all its vertices are collinear.
//
// @pre `CGAL::is_triangle_mesh(tmesh)`
//
// @tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
//
// @param tmesh the  triangulated surface mesh to be repaired
// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
//
// \cgalNamedParamsBegin
//   \cgalParamNBegin{vertex_point_map}
//     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
//     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
//                    as key type and `%Point_3` as value type}
//     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
//     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
//                     must be available in `TriangleMesh`.}
//   \cgalParamNEnd
//
//   \cgalParamNBegin{geom_traits}
//     \cgalParamDescription{an instance of a geometric traits class}
//     \cgalParamType{The traits class must provide the nested type `Point_3`,
//                    and the nested functors:
//                    - `Compare_distance_3` to compute the distance between 2 points
//                    - `Collinear_3` to check whether 3 points are collinear
//                    - `Less_xyz_3` to compare lexicographically two points
//                    - `Equal_3` to check whether 2 points are identical.
//                    For each functor Foo, a function `Foo foo_object()` must be provided.}
//     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
//     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
//   \cgalParamNEnd
// \cgalNamedParamsEnd
//
// \todo the function might not be able to remove all degenerate faces.
//       We should probably do something with the return type.
//
// \return `true` if all degenerate faces were successfully removed, and `false` otherwise.
template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_degenerate_faces(const FaceRange& face_range,
                             TriangleMesh& tmesh,
                             const NamedParameters& np)
{
  CGAL_assertion(CGAL::is_triangle_mesh(tmesh));

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef TriangleMesh                                                              TM;
  typedef typename boost::graph_traits<TriangleMesh>                                GT;
  typedef typename GT::vertex_descriptor                                            vertex_descriptor;
  typedef typename GT::halfedge_descriptor                                          halfedge_descriptor;
  typedef typename GT::edge_descriptor                                              edge_descriptor;
  typedef typename GT::face_descriptor                                              face_descriptor;

  typedef typename GetVertexPointMap<TM, NamedParameters>::type                     VertexPointMap;
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_property_map(vertex_point, tmesh));
  typedef typename GetGeomTraits<TM, NamedParameters>::type                         Traits;
  Traits traits = choose_parameter(get_parameter(np, internal_np::geom_traits), Traits());

  typedef typename boost::property_traits<VertexPointMap>::value_type               Point_3;
  typedef typename boost::property_traits<VertexPointMap>::reference                Point_ref;

  std::set<face_descriptor> degenerate_face_set;
  degenerate_faces(face_range, tmesh, std::inserter(degenerate_face_set, degenerate_face_set.begin()), np);

  const std::size_t faces_size = faces(tmesh).size();

  if(degenerate_face_set.empty())
    return true;

  if(degenerate_face_set.size() == faces_size)
  {
    clear(tmesh);
    return true;
  }

  // Sanitize the face range by adding adjacent degenerate faces
  const std::size_t range_size = face_range.size();
  bool is_range_full_mesh = (range_size == faces_size);
  if(!is_range_full_mesh)
  {
    std::list<face_descriptor> faces_to_visit(degenerate_face_set.begin(), degenerate_face_set.end());

    while(!faces_to_visit.empty())
    {
      face_descriptor fd = faces_to_visit.front();
      faces_to_visit.pop_front();

      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
        for(halfedge_descriptor inc_hd : halfedges_around_target(hd, tmesh))
        {
          face_descriptor adj_fd = face(inc_hd, tmesh);
          if(adj_fd == GT::null_face() || adj_fd == fd)
            continue;

          if(is_degenerate_triangle_face(adj_fd, tmesh))
          {
            if(degenerate_face_set.insert(adj_fd).second)
            {
              // successful insertion means we did not know about this face before
              faces_to_visit.push_back(adj_fd);
            }
          }
        }
      }
    }
  }

  // Note that there can't be any null edge incident to the degenerate faces range,
  // otherwise we would have a null face incident to the face range, and that is not possible
  // after the sanitization above
  std::set<edge_descriptor> edge_range;
  for(face_descriptor fd : degenerate_face_set)
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      edge_range.insert(edge(hd, tmesh));

  // First remove edges of length 0
  bool all_removed = remove_degenerate_edges(edge_range, tmesh, degenerate_face_set, np);

  CGAL_assertion_code(for(face_descriptor fd : degenerate_face_set) {)
  CGAL_assertion(is_degenerate_triangle_face(fd, tmesh));
  CGAL_assertion_code(})

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
  {
    std::cout <<"Done with null edges.\n";
    std::ofstream output("/tmp/no_null_edges.off");
    output << std::setprecision(17) << tmesh << "\n";
    output.close();
  }
#endif

  // Then, remove triangles made of 3 collinear points

  // start by filtering out border faces
  // @todo shall we avoid doing that in case a non-manifold vertex on the boundary or if a whole component disappear?
  std::set<face_descriptor> border_deg_faces;
  for(face_descriptor f : degenerate_face_set)
  {
    halfedge_descriptor h = halfedge(f, tmesh);
    for(int i=0; i<3; ++i)
    {
      if(is_border(opposite(h, tmesh), tmesh))
      {
        border_deg_faces.insert(f);
        break;
      }

      h = next(h, tmesh);
    }
  }

  while(!border_deg_faces.empty())
  {
    face_descriptor f_to_rm = *border_deg_faces.begin();
    border_deg_faces.erase(border_deg_faces.begin());

    halfedge_descriptor h = halfedge(f_to_rm, tmesh);
    for(int i=0; i<3; ++i)
    {
      face_descriptor f = face(opposite(h, tmesh), tmesh);
      if(f!=GT::null_face())
      {
        if(is_degenerate_triangle_face(f, tmesh, np))
          border_deg_faces.insert(f);
      }

      h = next(h, tmesh);
    }

    while(!is_border(opposite(h, tmesh), tmesh))
    {
      h = next(h, tmesh);
    }

    degenerate_face_set.erase(f_to_rm);
    Euler::remove_face(h, tmesh);
  }

  // Ignore faces with null edges
  if(!all_removed)
  {
    typename std::set<face_descriptor>::iterator it = degenerate_face_set.begin();
    while(it != degenerate_face_set.end())
    {
      bool has_degenerate_edge = false;
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(*it, tmesh), tmesh))
      {
        const edge_descriptor ed = edge(hd, tmesh);
        if(is_degenerate_edge(ed, tmesh, np))
        {
          has_degenerate_edge = true;
          it = degenerate_face_set.erase(it);
          break;
        }
      }

      if(!has_degenerate_edge)
        ++it;
    }
  }

  // first remove degree 3 vertices that are part of a cap
  // (only the vertex in the middle of the opposite edge)
  // This removal does not change the shape of the mesh.
  while(!degenerate_face_set.empty())
  {
    std::set<vertex_descriptor> vertices_to_remove;
    for(face_descriptor fd : degenerate_face_set)
    {
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, tmesh), tmesh))
      {
        vertex_descriptor vd = target(hd, tmesh);
        if(degree(vd, tmesh) == 3 && is_border(vd, tmesh)==GT::null_halfedge())
        {
          vertices_to_remove.insert(vd);
          break;
        }
      }
    }

    for(vertex_descriptor vd : vertices_to_remove)
    {
      for(halfedge_descriptor hd2 : halfedges_around_target(vd, tmesh))
        degenerate_face_set.erase(face(hd2, tmesh));

      // remove the central vertex and check if the new face is degenerated
      halfedge_descriptor hd = CGAL::Euler::remove_center_vertex(halfedge(vd, tmesh), tmesh);
      if(is_degenerate_triangle_face(face(hd, tmesh), tmesh, np))
      {
        degenerate_face_set.insert(face(hd, tmesh));
      }
    }

    if(vertices_to_remove.empty())
      break;
  }

  while(!degenerate_face_set.empty())
  {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
    std::cout << "Loop on removing deg faces\n";

    // ensure the mesh is not broken
    {
      std::ofstream out("/tmp/out.off");
      out << tmesh;
      out.close();

      std::vector<typename Traits::Point_3> points;
      std::vector<std::vector<std::size_t> > triangles;
      std::ifstream in("/tmp/out.off");
      CGAL::read_OFF(in, points, triangles);
      if(!CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(triangles))
      {
        std::cerr << "Warning: got a polygon soup (may simply be a non-manifold vertex)!\n";
      }
    }
#endif

    face_descriptor fd = *degenerate_face_set.begin();

    // look whether an incident triangle is also degenerate
    bool detect_cc_of_degenerate_triangles = false;
    for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, tmesh), tmesh))
    {
      face_descriptor adjacent_face = face(opposite(hd, tmesh), tmesh);
      if(adjacent_face!=GT::null_face() && degenerate_face_set.count(adjacent_face))
      {
        detect_cc_of_degenerate_triangles = true;
        break;
      }
    }

    if(!detect_cc_of_degenerate_triangles)
    {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      std::cout << "  no degenerate neighbors, using a flip.\n";
#endif
      degenerate_face_set.erase(degenerate_face_set.begin());

      // flip the longest edge of the triangle
      Point_ref p1 = get(vpmap, target(halfedge(fd, tmesh), tmesh));
      Point_ref p2 = get(vpmap, target(next(halfedge(fd, tmesh), tmesh), tmesh));
      Point_ref p3 = get(vpmap, source(halfedge(fd, tmesh), tmesh));

      CGAL_assertion(p1!=p2 && p1!=p3 && p2!=p3);

      typename Traits::Compare_distance_3 compare_distance = traits.compare_distance_3_object();

      halfedge_descriptor edge_to_flip;
      if(compare_distance(p1,p2, p1,p3) != CGAL::SMALLER) // p1p2 > p1p3
      {
        if(compare_distance(p1,p2, p2,p3) != CGAL::SMALLER) // p1p2 > p2p3
          // flip p1p2
          edge_to_flip = next(halfedge(fd, tmesh), tmesh);
        else
          // flip p2p3
          edge_to_flip = prev(halfedge(fd, tmesh), tmesh);
      }
      else
      {
        if(compare_distance(p1,p3, p2,p3) != CGAL::SMALLER) // p1p3>p2p3
          //flip p3p1
          edge_to_flip = halfedge(fd, tmesh);
        else
          //flip p2p3
          edge_to_flip = prev(halfedge(fd, tmesh), tmesh);
      }

      face_descriptor opposite_face=face(opposite(edge_to_flip, tmesh), tmesh);
      if(opposite_face == GT::null_face())
      {
        // simply remove the face
        Euler::remove_face(edge_to_flip, tmesh);
      }
      else
      {
        // condition for the flip to be valid (the edge to be created do not already exists)
        if(!halfedge(target(next(edge_to_flip, tmesh), tmesh),
                       target(next(opposite(edge_to_flip, tmesh), tmesh), tmesh),
                       tmesh).second)
        {
          Euler::flip_edge(edge_to_flip, tmesh);
        }
        else
        {
          all_removed = false;
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
          std::cout << "  WARNING: flip is not possible\n";
          // @todo Let p and q be the vertices opposite to `edge_to_flip`, and let
          //       r be the vertex of `edge_to_flip` that is the furthest away from
          //       the edge `pq`. In that case I think we should remove all the triangles
          //       so that the triangle pqr is in the mesh.
#endif
        }
      }
    }
    else
    {
      // Process a connected component of degenerate faces
      // get all the faces from the connected component
      // and the boundary edges
      std::set<face_descriptor> cc_faces;
      std::vector<face_descriptor> queue;
      std::vector<halfedge_descriptor> boundary_hedges;
      std::vector<halfedge_descriptor> inside_hedges;
      queue.push_back(fd);
      cc_faces.insert(fd);

      while(!queue.empty())
      {
        face_descriptor top=queue.back();
        queue.pop_back();
        for(halfedge_descriptor hd : halfedges_around_face(halfedge(top, tmesh), tmesh))
        {
          face_descriptor adjacent_face = face(opposite(hd, tmesh), tmesh);
          if(adjacent_face==GT::null_face() || degenerate_face_set.count(adjacent_face)==0)
          {
            boundary_hedges.push_back(hd);
          }
          else
          {
            if(cc_faces.insert(adjacent_face).second)
              queue.push_back(adjacent_face);

            if(hd < opposite(hd, tmesh))
              inside_hedges.push_back(hd);
          }
        }
      }

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      std::cout << "  Deal with a cc of " << cc_faces.size() << " degenerate faces.\n";
      /// dump cc_faces
      {
        int id = 0;
        std::map<vertex_descriptor, int> vids;
        for(face_descriptor f : cc_faces)
        {
          if(vids.insert(std::make_pair(target(halfedge(f, tmesh), tmesh), id)).second) ++id;
          if(vids.insert(std::make_pair(target(next(halfedge(f, tmesh), tmesh), tmesh), id)).second) ++id;
          if(vids.insert(std::make_pair(target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh), id)).second) ++id;
        }

        std::ofstream output("/tmp/cc_faces.off");
        output << std::setprecision(17);
        output << "OFF\n" << vids.size() << " " << cc_faces.size() << " 0\n";
        std::vector<typename Traits::Point_3> points(vids.size());
        typedef std::pair<const vertex_descriptor, int> Pair_type;
        for(const Pair_type& p : vids)

          points[p.second] = get(vpmap, p.first);
        for(const Point_3& p : points)
          output << p << "\n";

        for(face_descriptor f : cc_faces)
        {
          output << "3 "
                 << vids[target(halfedge(f, tmesh), tmesh)] << " "
                 << vids[target(next(halfedge(f, tmesh), tmesh), tmesh)] << " "
                 << vids[target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh)] << "\n";
        }

        for(std::size_t pid=2; pid!=points.size(); ++pid)
        {
          CGAL_assertion(collinear(points[0], points[1], points[pid]));
        }
      }
#endif

      // find vertices strictly inside the cc
      std::set<vertex_descriptor> boundary_vertices;
      for(halfedge_descriptor hd : boundary_hedges)
        boundary_vertices.insert(target(hd, tmesh));

      std::set<vertex_descriptor> inside_vertices;
      for(halfedge_descriptor hd : inside_hedges)
      {
        if(!boundary_vertices.count(target(hd, tmesh)))
          inside_vertices.insert(target(hd, tmesh));
        if(!boundary_vertices.count(source(hd, tmesh)))
          inside_vertices.insert(source(hd, tmesh));
      }

      // v-e+f = 1 for a topological disk and e = (3f+#boundary_edges)/2
      if(boundary_vertices.size()+inside_vertices.size() -
          (cc_faces.size()+boundary_hedges.size())/2 != 1)
      {
        //cc_faces does not define a topological disk
        // @todo Find to way to handle that case
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "  WARNING: Cannot remove the component of degenerate faces: not a topological disk.\n";
#endif

        for(face_descriptor f : cc_faces)
          degenerate_face_set.erase(f);

        continue;
      }

      // preliminary step to check if the operation is possible
      // sort the boundary points along the common supporting line
      //    we first need a reference point
      typedef internal::Less_vertex_point<TriangleMesh, VertexPointMap, Traits> Less_vertex;
      std::pair<
          typename std::set<vertex_descriptor>::iterator,
          typename std::set<vertex_descriptor>::iterator > ref_vertices =
          boost::minmax_element(boundary_vertices.begin(),
                                boundary_vertices.end(),
                                Less_vertex(traits, vpmap));

      //    and then we sort the vertices using this reference point
      typedef std::set<Point_3>                                Sorted_point_set;
      Sorted_point_set sorted_points;
      for(vertex_descriptor v : boundary_vertices)
        sorted_points.insert(get(vpmap, v));

      CGAL_assertion(sorted_points.size()==
                     std::set<typename Traits::Point_3>(sorted_points.begin(),
                                                        sorted_points.end()).size());

      CGAL_assertion(get(vpmap, *ref_vertices.first) == *sorted_points.begin());
      CGAL_assertion(get(vpmap, *ref_vertices.second) == *std::prev(sorted_points.end()));

      const Point_3& xtrm1 = *sorted_points.begin();
      const Point_3& xtrm2 = *std::prev(sorted_points.end());

      // recover halfedges on the hole, bounded by the reference vertices
      std::vector<halfedge_descriptor> side_one, side_two;

      // look for the outgoing border halfedge of the first extreme point
      for(halfedge_descriptor hd : boundary_hedges)
      {
        if(get(vpmap, source(hd, tmesh)) == xtrm1)
        {
          side_one.push_back(hd);
          break;
        }
      }

      CGAL_assertion(side_one.size() == 1);

      bool non_monotone_border = false;

      while(get(vpmap, target(side_one.back(), tmesh)) != xtrm2)
      {
        vertex_descriptor prev_vertex = target(side_one.back(), tmesh);
        for(halfedge_descriptor hd : boundary_hedges)
        {
          if(source(hd, tmesh) == prev_vertex)
          {
            if(get(vpmap, target(hd, tmesh)) < get(vpmap, prev_vertex))
              non_monotone_border = true;

            side_one.push_back(hd);
            break;
          }
        }

        if(non_monotone_border)
          break;
      }

      if(non_monotone_border)
      {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "  WARNING: Cannot remove the component of degenerate faces: border not a monotonic cycle.\n";
#endif

        for(face_descriptor f : cc_faces)
          degenerate_face_set.erase(f);

        continue;
      }

      // look for the outgoing border halfedge of second extreme vertex
      for(halfedge_descriptor hd : boundary_hedges)
      {
        if(source(hd, tmesh) == target(side_one.back(), tmesh))
        {
          side_two.push_back(hd);
          break;
        }
      }

      CGAL_assertion(side_two.size() == 1);

      while(target(side_two.back(), tmesh) != source(side_one.front(), tmesh))
      {
        vertex_descriptor prev_vertex = target(side_two.back(), tmesh);
        for(halfedge_descriptor hd : boundary_hedges)
        {
          if(source(hd, tmesh) == prev_vertex)
          {
            if(get(vpmap, target(hd, tmesh)) > get(vpmap, prev_vertex))
              non_monotone_border = true;

            side_two.push_back(hd);
            break;
          }
        }

        if(non_monotone_border)
          break;
      }

      if(non_monotone_border)
      {
#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "  WARNING: Cannot remove the component of degenerate faces: border not a monotonic cycle.\n";
#endif

        for(face_descriptor f : cc_faces)
          degenerate_face_set.erase(f);

        continue;
      }

      CGAL_assertion(side_one.size()+side_two.size()==boundary_hedges.size());

      // reverse the order of the second side so as to follow
      // the same order than side one
      std::reverse(side_two.begin(), side_two.end());
      for(halfedge_descriptor& h : side_two)
        h = opposite(h, tmesh);

      //make sure the points of the vertices along side_one are correctly sorted
      std::vector<Point_3> side_points;
      side_points.reserve(side_one.size()+1);
      side_points.push_back(get(vpmap, source(side_one.front(), tmesh)));

      for(halfedge_descriptor h : side_one)
        side_points.push_back(get(vpmap, target(h, tmesh)));

      CGAL_assertion(get(vpmap,source(side_one.front(), tmesh)) == side_points.front());
      CGAL_assertion(get(vpmap,target(side_one.back(), tmesh)) == side_points.back());

      // @todo the reordering could lead to the apparition of null edges.
      std::sort(side_points.begin(), side_points.end());

      CGAL_assertion(std::unique(side_points.begin(), side_points.end())==side_points.end());
      for(std::size_t i=0; i<side_one.size()-1; ++i)
        put(vpmap, target(side_one[i], tmesh), side_points[i+1]);

      //same thing for side_two
      side_points.clear();
      side_points.reserve(side_two.size()+1);
      side_points.push_back(get(vpmap, source(side_two.front(), tmesh)));

      for(halfedge_descriptor h : side_two)
        side_points.push_back(get(vpmap,target(h, tmesh)));

      CGAL_assertion(get(vpmap, source(side_two.front(), tmesh)) == side_points.front());
      CGAL_assertion(get(vpmap, target(side_two.back(), tmesh)) == side_points.back());

      // @todo the reordering could lead to the apparition of null edges.
      std::sort(side_points.begin(), side_points.end());

      CGAL_assertion(std::unique(side_points.begin(), side_points.end())==side_points.end());
      for(std::size_t i=0; i<side_two.size()-1; ++i)
        put(vpmap, target(side_two[i], tmesh), side_points[i+1]);

      CGAL_assertion(source(side_one.front(), tmesh) == *ref_vertices.first);
      CGAL_assertion(source(side_two.front(), tmesh) == *ref_vertices.first);
      CGAL_assertion(target(side_one.back(), tmesh) == *ref_vertices.second);
      CGAL_assertion(target(side_two.back(), tmesh) == *ref_vertices.second);

      typename Sorted_point_set::iterator it_pt = std::next(sorted_points.begin()),
                                          it_pt_end = std::prev(sorted_points.end());

      bool non_collapsable = false;
      typename std::vector<halfedge_descriptor>::iterator side_one_it = side_one.begin();
      typename std::vector<halfedge_descriptor>::iterator side_two_it = side_two.begin();
      for(;it_pt!=it_pt_end;++it_pt)
      {
        // check if it_pt is the point of the target of one or two halfedges
        bool target_of_side_one = (get(vpmap, target(*side_one_it, tmesh)) == *it_pt);
        bool target_of_side_two = (get(vpmap, target(*side_two_it, tmesh)) == *it_pt);

        if(target_of_side_one && target_of_side_two)
        {
          for(halfedge_descriptor h : halfedges_around_target(*side_one_it, tmesh))
          {
            if(source(h, tmesh) == target(*side_two_it, tmesh))
            {
              non_collapsable = true;
              break;
            }
          }
        }
        else
        {
          CGAL_assertion(target_of_side_one || target_of_side_two);
          vertex_descriptor v1 = target_of_side_one ? target(*side_one_it, tmesh)
                                                    : target(*side_two_it, tmesh);
          vertex_descriptor v2 = target_of_side_two ? target(next(opposite(*side_one_it, tmesh), tmesh), tmesh)
                                                    : target(next(*side_two_it, tmesh), tmesh);
          for(halfedge_descriptor h : halfedges_around_target(v1, tmesh))
          {
            if(source(h, tmesh)==v2)
            {
              non_collapsable=true;
              break;
            }
          }
        }

        if(non_collapsable) break;
        if(target_of_side_one) ++side_one_it;
        if(target_of_side_two) ++side_two_it;
      }

      if(non_collapsable)
      {
        for(face_descriptor f : cc_faces)
          degenerate_face_set.erase(f);

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
        std::cout << "  WARNING: cannot remove a connected components of degenerate faces.\n";
#endif
        continue;
      }

      // now proceed to the fix
      // update the face and halfedge vertex pointers on the boundary
      for(halfedge_descriptor h : boundary_hedges)
      {
        set_face(h, GT::null_face(), tmesh);
        set_halfedge(target(h, tmesh), h, tmesh);
      }

      // update next/prev pointers of boundary_hedges
      for(halfedge_descriptor h : boundary_hedges)
      {
        halfedge_descriptor next_candidate = next(h, tmesh);
        while(face(next_candidate, tmesh)!=GT::null_face())
          next_candidate = next(opposite(next_candidate, tmesh), tmesh);

        set_next(h, next_candidate, tmesh);
      }

      // remove degenerate faces
      for(face_descriptor f : cc_faces)
      {
        degenerate_face_set.erase(f);
        remove_face(f, tmesh);
      }

      // remove interior edges
      for(halfedge_descriptor h : inside_hedges)
        remove_edge(edge(h, tmesh), tmesh);

      // remove interior vertices
      for(vertex_descriptor v : inside_vertices)
        remove_vertex(v, tmesh);

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      std::cout << "  side_one.size() " << side_one.size() << "\n";
      std::cout << "  side_two.size() " << side_two.size() << "\n";
#endif

      CGAL_assertion(source(side_one.front(), tmesh) == *ref_vertices.first);
      CGAL_assertion(source(side_two.front(), tmesh) == *ref_vertices.first);
      CGAL_assertion(target(side_one.back(), tmesh) == *ref_vertices.second);
      CGAL_assertion(target(side_two.back(), tmesh) == *ref_vertices.second);

      // now split each side to contains the same sequence of points
      //    first side
      int hi = 0;

      for(typename Sorted_point_set::iterator it=std::next(sorted_points.begin()),
                                              it_end=sorted_points.end(); it!=it_end; ++it)
      {
        CGAL_assertion(*std::prev(it) == get(vpmap, source(side_one[hi], tmesh)));

        if(*it != get(vpmap, target(side_one[hi], tmesh)))
        {
          // split the edge and update the point
          halfedge_descriptor h1 = next(opposite(side_one[hi], tmesh), tmesh);
          put(vpmap, target(Euler::split_edge(side_one[hi], tmesh), tmesh), *it);

          // split_edge updates the halfedge of the source vertex of h,
          // since we reuse later the halfedge of the first refernce vertex
          // we must set it as we need.
          if(source(h1, tmesh) == *ref_vertices.first)
            set_halfedge(*ref_vertices.first, prev(prev(side_one[hi], tmesh), tmesh), tmesh);

          // retriangulate the opposite face
          if(face(h1, tmesh) != GT::null_face())
            Euler::split_face(h1, opposite(side_one[hi], tmesh), tmesh);
        }
        else
        {
          ++hi;
        }
      }

      //    second side
      hi = 0;
      for(typename Sorted_point_set::iterator it=std::next(sorted_points.begin()),
                                              it_end=sorted_points.end(); it!=it_end; ++it)
      {
        CGAL_assertion(*std::prev(it) == get(vpmap, source(side_two[hi], tmesh)));
        if(*it != get(vpmap, target(side_two[hi], tmesh)))
        {
          // split the edge and update the point
          halfedge_descriptor h2 = Euler::split_edge(side_two[hi], tmesh);
          put(vpmap, target(h2, tmesh), *it);

          // split_edge updates the halfedge of the source vertex of h,
          // since we reuse later the halfedge of the first refernce vertex
          // we must set it as we need.
          if(source(h2, tmesh) == *ref_vertices.first)
            set_halfedge(*ref_vertices.first, opposite(h2, tmesh), tmesh);

          // retriangulate the face
          if(face(h2, tmesh) != GT::null_face())
            Euler::split_face(h2, next(side_two[hi], tmesh), tmesh);
        }
        else
        {
          ++hi;
        }
      }

      CGAL_assertion(target(halfedge(*ref_vertices.first, tmesh), tmesh) == *ref_vertices.first);
      CGAL_assertion(face(halfedge(*ref_vertices.first, tmesh), tmesh) == GT::null_face());

#ifdef CGAL_PMP_REMOVE_DEGENERATE_FACES_DEBUG
      {
        halfedge_descriptor h_side2 = halfedge(*ref_vertices.first, tmesh);
        halfedge_descriptor h_side1 = next(h_side2, tmesh);

        do
        {
          CGAL_assertion(get(vpmap, source(h_side1, tmesh)) == get(vpmap, target(h_side2, tmesh)));
          CGAL_assertion(get(vpmap, target(h_side1, tmesh)) == get(vpmap, source(h_side2, tmesh)));

          if(target(next(opposite(h_side1, tmesh), tmesh), tmesh) ==
               target(next(opposite(h_side2, tmesh), tmesh), tmesh))
          {
            CGAL_assertion(!"Forbidden simplification");
          }

          h_side2 = prev(h_side2, tmesh);
          h_side1 = next(h_side1, tmesh);
        }
        while(target(h_side1, tmesh) != *ref_vertices.second);
      }
#endif

      // remove side1 and replace its opposite hedges by those of side2
      halfedge_descriptor h_side2 = halfedge(*ref_vertices.first, tmesh);
      halfedge_descriptor h_side1 = next(h_side2, tmesh);
      for(;;)
      {
        CGAL_assertion(get(vpmap, source(h_side1, tmesh)) == get(vpmap, target(h_side2, tmesh)));
        CGAL_assertion(get(vpmap, target(h_side1, tmesh)) == get(vpmap, source(h_side2, tmesh)));

        // backup target vertex
        vertex_descriptor vertex_to_remove = target(h_side1, tmesh);
        if(vertex_to_remove!=*ref_vertices.second)
        {
          vertex_descriptor replacement_vertex = source(h_side2, tmesh);
          // replace the incident vertex
          for(halfedge_descriptor hd : halfedges_around_target(h_side1, tmesh))
            set_target(hd, replacement_vertex, tmesh);
        }

        // prev side2 hedge for next loop
        halfedge_descriptor h_side2_for_next_turn = prev(h_side2, tmesh);

        // replace the opposite of h_side1 by h_side2
        halfedge_descriptor opposite_h_side1 = opposite(h_side1, tmesh);
        face_descriptor the_face = face(opposite_h_side1, tmesh);
        set_face(h_side2, the_face, tmesh);

        if(the_face!=GT::null_face())
          set_halfedge(the_face, h_side2, tmesh);

        set_next(h_side2, next(opposite_h_side1, tmesh), tmesh);
        set_next(prev(opposite_h_side1, tmesh), h_side2, tmesh);

        // take the next hedges
        edge_descriptor edge_to_remove = edge(h_side1, tmesh);
        h_side1 = next(h_side1, tmesh);

        // now remove the extra edge
        remove_edge(edge_to_remove, tmesh);

        // ... and the extra vertex if it's not the second reference
        if(vertex_to_remove == *ref_vertices.second)
        {
          // update the halfedge pointer of the last vertex (others were already from side 2)
          CGAL_assertion(target(opposite(h_side2, tmesh), tmesh) == vertex_to_remove);
          set_halfedge(vertex_to_remove, opposite(h_side2, tmesh), tmesh);
          break;
        }
        else
        {
          remove_vertex(vertex_to_remove , tmesh);
        }

        h_side2 = h_side2_for_next_turn;
      }
    }
  }

  return all_removed;
}

template <typename FaceRange, typename TriangleMesh>
bool remove_degenerate_faces(const FaceRange& face_range,
                             TriangleMesh& tmesh)
{
  return remove_degenerate_faces(face_range, tmesh, CGAL::parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_degenerate_faces(TriangleMesh& tmesh,
                             const CGAL_PMP_NP_CLASS& np)
{
  return remove_degenerate_faces(faces(tmesh), tmesh, np);
}

template<typename TriangleMesh>
bool remove_degenerate_faces(TriangleMesh& tmesh)
{
  return remove_degenerate_faces(tmesh, CGAL::parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_DEGENERACIES_H
