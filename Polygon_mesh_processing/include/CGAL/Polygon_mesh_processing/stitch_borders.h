// Copyright (c) 2014 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot


#ifndef CGAL_STITCH_POLYGON_MESH_H
#define CGAL_STITCH_POLYGON_MESH_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/array.h>
#include <CGAL/Union_find.h>
#include <CGAL/utility.h>

#include <map>
#include <vector>
#include <utility>
#include <boost/range.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL{

namespace Polygon_mesh_processing{

namespace internal{

template <typename PM, typename VertexPointMap>
struct Less_for_halfedge
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor
    halfedge_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::reference Point;

  Less_for_halfedge(const PM& pmesh_,
                    const VertexPointMap& vpmap_)
    : pmesh(pmesh_),
      vpmap(vpmap_)
  {}

  bool operator()(halfedge_descriptor h1,
                  halfedge_descriptor h2) const
  {
    Point s1 = get(vpmap,target(opposite(h1, pmesh), pmesh));
    Point t1 = get(vpmap,target(h1, pmesh));
    Point s2 = get(vpmap,target(opposite(h2, pmesh), pmesh));
    Point t2 = get(vpmap,target(h2, pmesh));
    return
    ( s1 < t1?  std::make_pair(s1,t1) : std::make_pair(t1, s1) )
    <
    ( s2 < t2?  std::make_pair(s2,t2) : std::make_pair(t2, s2) );
  }

  const PM& pmesh;
  const VertexPointMap& vpmap;
};

template <typename PM, typename OutputIterator, typename LessHedge, typename VertexPointMap>
OutputIterator
collect_duplicated_stitchable_boundary_edges
(PM& pmesh, OutputIterator out, LessHedge less_hedge, const VertexPointMap& vpmap)
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef std::map<halfedge_descriptor, std::pair<int, std::size_t>, LessHedge> Border_halfedge_map;
  Border_halfedge_map border_halfedge_map(less_hedge);

  std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > halfedge_pairs;
  std::vector< bool > manifold_halfedge_pairs;

  BOOST_FOREACH(halfedge_descriptor he, halfedges(pmesh))
  {
    if ( !CGAL::is_border(he, pmesh) )
      continue;

    typename Border_halfedge_map::iterator set_it;
    bool insertion_ok;
    CGAL::cpp11::tie(set_it, insertion_ok)
      = border_halfedge_map.insert(std::make_pair(he,std::make_pair(1,0)));

    if ( !insertion_ok ){ // we found already a halfedge with the points
      ++set_it->second.first; // increase the multiplicity
      if(set_it->second.first == 2)
      {
        if ( get(vpmap, source(he,pmesh))==get(vpmap, target(set_it->first,pmesh)) &&
             get(vpmap, target(he,pmesh))==get(vpmap, source(set_it->first,pmesh)) )
        {
          set_it->second.second = halfedge_pairs.size(); // set the id of the pair in the vector
          halfedge_pairs.push_back( std::make_pair(set_it->first, he) );
          manifold_halfedge_pairs.push_back(true);
        }
      }
      else
        if ( set_it->second.first > 2 )
          manifold_halfedge_pairs[ set_it->second.second ] = false;
    }
  }

  // put in `out` only manifold edges from the set of edges to stitch.
  // We choose not to allow only a pair out of the whole set to be stitched
  // as we can produce inconsistent stitching along a sequence of non-manifold edges
  std::size_t nb_pairs=halfedge_pairs.size();
  for (std::size_t i=0; i<nb_pairs; ++i)
  {
    if( manifold_halfedge_pairs[i] )
      *out++=halfedge_pairs[i];
  }

  return out;
}


template <class PM>
void update_target_vertex(typename boost::graph_traits<PM>::halfedge_descriptor h,
                          typename boost::graph_traits<PM>::vertex_descriptor v_kept,
                          PM& pmesh)
{
  typename boost::graph_traits<PM>::halfedge_descriptor start = h;
  do{
    set_target(h, v_kept, pmesh);
    h = opposite(next(h, pmesh), pmesh);
  } while( h != start );
}

template <class vertex_descriptor, class Handle_map>
typename Union_find<vertex_descriptor>::handle
uf_get_handle(vertex_descriptor v,
              Union_find<vertex_descriptor>& uf_vertices,
              Handle_map& handles)
{
  std::pair<typename Handle_map::iterator, bool> insert_res =
  handles.insert( std::make_pair(v, typename Union_find<vertex_descriptor>::handle()) );
  if( insert_res.second )
  {
    insert_res.first->second=uf_vertices.make_set(v);
  }
  return insert_res.first->second;
}

template <class vertex_descriptor, class Handle_map>
void uf_join_vertices(vertex_descriptor v1, vertex_descriptor v2,
                      Union_find<vertex_descriptor>& uf_vertices,
                      Handle_map& handles)
{
  typename Union_find<vertex_descriptor>::handle
    h1 = uf_get_handle(v1, uf_vertices, handles),
    h2 = uf_get_handle(v2, uf_vertices, handles);
  uf_vertices.unify_sets(h1, h2);
}



// main functions (vertices to keep selected and halfedge pairs filtered)
template <class PM, typename HalfedgePairsRange,
          typename Uf_vertices, typename Uf_handles>
void run_stitch_borders(PM& pmesh,
                        const HalfedgePairsRange& to_stitch,
                        Uf_vertices& uf_vertices,
                        Uf_handles& uf_handles)
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename std::pair<halfedge_descriptor, halfedge_descriptor> halfedges_pair;

  std::vector<vertex_descriptor> vertices_to_delete;
  BOOST_FOREACH(const halfedges_pair hk, to_stitch)
  {
    halfedge_descriptor h1 = hk.first;
    halfedge_descriptor h2 = hk.second;

    vertex_descriptor h1_tgt = target(h1, pmesh);
    vertex_descriptor h2_src = source(h2, pmesh);

    //update vertex pointers: target of h1 vs source of h2
    vertex_descriptor v_to_keep =
      *uf_vertices.find(uf_get_handle(h1_tgt, uf_vertices, uf_handles));

    if (v_to_keep!=h1_tgt){
      vertices_to_delete.push_back(h1_tgt);
      update_target_vertex(h1, v_to_keep, pmesh);
    }
    if (v_to_keep != h2_src && h1_tgt!=h2_src)
    {
      vertices_to_delete.push_back( h2_src );
      update_target_vertex(opposite(h2, pmesh), v_to_keep, pmesh);
    }
    set_halfedge(v_to_keep, h1, pmesh);

    vertex_descriptor h1_src = source(h1, pmesh);
    vertex_descriptor h2_tgt = target(h2, pmesh);

    //update vertex pointers: target of h1 vs source of h2
    v_to_keep = *uf_vertices.find(uf_get_handle(h2_tgt, uf_vertices, uf_handles));
    if (v_to_keep!=h2_tgt)
    {
      vertices_to_delete.push_back( h2_tgt );
      update_target_vertex(h2, v_to_keep, pmesh);
    }
    if (v_to_keep!=h1_src && h1_src!=h2_tgt)
    {
      vertices_to_delete.push_back( h1_src );
      update_target_vertex(opposite(h1, pmesh), v_to_keep, pmesh);
    }
    set_halfedge(v_to_keep, opposite(h1,pmesh), pmesh);
  }

  /// Update next/prev of neighbor halfedges (that are not set for stiching)
  /// _______   _______
  ///        | |
  ///        | |
  /// In order to avoid having to maintain a set with halfedges to stitch
  /// we do on purpose next-prev linking that might not be useful but that
  /// is harmless and still less expensive than doing queries in a set
  BOOST_FOREACH(const halfedges_pair hk, to_stitch)
  {
    halfedge_descriptor h1 = hk.first;
    halfedge_descriptor h2 = hk.second;

    //link h2->prev() to h1->next()
    halfedge_descriptor pr = prev(h2, pmesh);
    halfedge_descriptor nx = next(h1, pmesh);
    set_next(pr, nx, pmesh);

    //link h1->prev() to h2->next()
    pr = prev(h1, pmesh);
    nx = next(h2, pmesh);
    set_next(pr, nx, pmesh);
  }

  /// update HDS connectivity, removing the second halfedge
  /// of each the pair and its opposite
  BOOST_FOREACH(const halfedges_pair hk, to_stitch)
  {
    halfedge_descriptor h1 = hk.first;
    halfedge_descriptor h2 = hk.second;

  ///Set face-halfedge relationship
    //h2 and its opposite will be removed
    set_face(h1, face(opposite(h2, pmesh), pmesh), pmesh);
    set_halfedge(face(h1, pmesh), h1, pmesh);
    //update next/prev pointers
    halfedge_descriptor tmp = prev(opposite(h2, pmesh), pmesh);
    set_next(tmp, h1, pmesh);
    tmp = next(opposite(h2, pmesh), pmesh);
    set_next(h1, tmp, pmesh);

  /// remove the extra halfedges
    remove_edge(edge(h2, pmesh), pmesh);
  }

  //remove the extra vertices
  BOOST_FOREACH(vertex_descriptor vd, vertices_to_delete)
  {
    remove_vertex(vd, pmesh);
  }
}

template <class PM, typename HalfedgePairsRange>
void stitch_borders_impl(PM& pmesh,
                         const HalfedgePairsRange& to_stitch)
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename std::pair<halfedge_descriptor, halfedge_descriptor> halfedges_pair;

  // The first step of the algorithm is to filter halfedges to be stitched so that
  // after stitching no edges will be present more than once.

  // Merge the vertices
  typedef CGAL::Union_find<vertex_descriptor> Uf_vertices;
  Uf_vertices uf_vertices;
  typedef boost::unordered_map<vertex_descriptor, typename Uf_vertices::handle> Uf_handles;
  Uf_handles uf_handles;

  BOOST_FOREACH(const halfedges_pair hk, to_stitch)
  {
    halfedge_descriptor h1 = hk.first;
    halfedge_descriptor h2 = hk.second;

    CGAL_assertion(CGAL::is_border(h1, pmesh));
    CGAL_assertion(CGAL::is_border(h2, pmesh));
    CGAL_assertion(!CGAL::is_border(opposite(h1, pmesh), pmesh));
    CGAL_assertion(!CGAL::is_border(opposite(h2, pmesh), pmesh));

    vertex_descriptor tgt1 = target(h1, pmesh), src1 = source(h1, pmesh);
    vertex_descriptor src2 = source(h2, pmesh), tgt2 = target(h2, pmesh);
    uf_join_vertices(tgt1, src2, uf_vertices, uf_handles);
    uf_join_vertices(src1, tgt2, uf_vertices, uf_handles);
  }

  // detect vertices that cannot be stitched because it would produce a non-manifold edge
  // We look for vertex to be stitched and collect all incident edges with another endpoint
  // to be stitched (that is not an edge scheduled for stitching). That way we can detect
  // if more that one edge will share the same two "master" endpoints.
  typedef boost::unordered_map< std::pair<vertex_descriptor, vertex_descriptor>,
                                std::vector<halfedge_descriptor> > Halfedges_after_stitching;
  Halfedges_after_stitching halfedges_after_stitching;

  typedef std::pair<const vertex_descriptor, typename Uf_vertices::handle> Pair_type;
  BOOST_FOREACH(const Pair_type p, uf_handles)
  {
    vertex_descriptor vd=p.first;
    typename Uf_vertices::handle tgt_handle = uf_vertices.find(uf_handles[vd]);
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(vd, pmesh))
    {
      vertex_descriptor other_vd = source(hd, pmesh);

      typename Uf_handles::iterator it_res = uf_handles.find(other_vd);

      if (it_res!=uf_handles.end()) // if the other vertex is also involved in a merge
      {
        if (other_vd < vd) continue; // avoid reporting twice the same edge
        typename Uf_vertices::handle src_handle=uf_vertices.find(it_res->second);
        halfedges_after_stitching[make_sorted_pair(*tgt_handle, *src_handle)].push_back(hd);
      }
      else
        halfedges_after_stitching[make_sorted_pair(*tgt_handle, other_vd)].push_back(hd);
    }
  }

  // look for edges that will be present more than once after the stitching
  // (no edges scheduled for stitching involved)
  boost::unordered_set<vertex_descriptor> unstitchable_vertices;
  for (typename Halfedges_after_stitching::iterator it=halfedges_after_stitching.begin(),
                                                    it_end=halfedges_after_stitching.end();
                                                    it!=it_end; ++it)
  {
    switch (it->second.size())
    {
      case 1:
       break; // nothing to do
      case 2:
      {
        if (is_border_edge(it->second.front(), pmesh) &&
            is_border_edge(it->second.back(), pmesh))
          break; // these are edges that are most possibly scheduled for stitching or will create a two halfedge loop
        CGAL_FALLTHROUGH;
      }
      default:
      {
        // this is a bit extreme as maybe some could be stitched
        // (but safer because the master could be one of them)
        BOOST_FOREACH(halfedge_descriptor hd, it->second)
        {
          unstitchable_vertices.insert(source(hd, pmesh));
          unstitchable_vertices.insert(target(hd, pmesh));
        }
      }
    }
  }

  // filter halfedges to stitch
  if (!unstitchable_vertices.empty())
  {
    std::vector<halfedges_pair> to_stitch_filtered;
    to_stitch_filtered.reserve( to_stitch.size());
    BOOST_FOREACH(const halfedges_pair hk, to_stitch)
    {
      // We test both halfedges because the previous test
      // might involve only one of the two halfedges
      if ( unstitchable_vertices.count( source(hk.first, pmesh) )== 0 &&
           unstitchable_vertices.count( target(hk.first, pmesh) )== 0 &&
           unstitchable_vertices.count( source(hk.second, pmesh) )== 0 &&
           unstitchable_vertices.count( target(hk.second, pmesh) )== 0 )
      {
        to_stitch_filtered.push_back(hk);
      }
    }

    // redo union find as some "master" vertex might be unstitchable
    uf_vertices.clear();
    uf_handles.clear();
    BOOST_FOREACH(const halfedges_pair hk, to_stitch_filtered)
    {
      halfedge_descriptor h1 = hk.first;
      halfedge_descriptor h2 = hk.second;

      vertex_descriptor tgt1 = target(h1, pmesh), src1 = source(h1, pmesh);
      vertex_descriptor src2 = source(h2, pmesh), tgt2 = target(h2, pmesh);
      uf_join_vertices(tgt1, src2, uf_vertices, uf_handles);
      uf_join_vertices(src1, tgt2, uf_vertices, uf_handles);
    }

    run_stitch_borders(pmesh, to_stitch_filtered, uf_vertices, uf_handles);
  }
  else
    run_stitch_borders(pmesh, to_stitch, uf_vertices, uf_handles);
}

template <class PM>
void stitch_boundary_cycle_2(PM& pmesh)
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;

  std::vector<halfedge_descriptor> cycles;
  BOOST_FOREACH(halfedge_descriptor hd, halfedges(pmesh))
  {
    if ( is_border(hd, pmesh) )
    {
      if ( hd < next(hd, pmesh) && next(next(hd, pmesh), pmesh) == hd )
      {
        cycles.push_back(hd);
      }
    }
  }

  BOOST_FOREACH(halfedge_descriptor hd, cycles)
  {
    halfedge_descriptor nhd = next(hd, pmesh);

    //nhd and its opposite will be removed
    //update face pointer
    set_face(hd, face(opposite(nhd, pmesh), pmesh), pmesh);
    set_halfedge(face(hd, pmesh), hd, pmesh);
    //update next/prev pointers
    halfedge_descriptor tmp = prev(opposite(nhd, pmesh), pmesh);
    set_next(tmp, hd, pmesh);
    tmp = next(opposite(nhd, pmesh), pmesh);
    set_next(hd, tmp, pmesh);
    //update vertex pointers
    set_halfedge(source(hd, pmesh), opposite(hd, pmesh), pmesh);
    set_halfedge(target(hd, pmesh), hd, pmesh);

    // remove the extra halfedges
    remove_edge(edge(nhd, pmesh), pmesh);
  }
}

} //end of namespace internal



/*!
* \ingroup PMP_repairing_grp
* Stitches together border halfedges in a polygon mesh.
* The halfedges to be stitched are provided in `hedge_pairs_to_stitch`.
* For each pair `p` in this vector, `p.second` and its opposite will be removed
* from `pmesh`.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
* @tparam HalfedgePairsRange a range of
*         `std::pair<boost::graph_traits<PolygonMesh>::%halfedge_descriptor,
*         boost::graph_traits<PolygonMesh>::%halfedge_descriptor>`,
*         model of `Range`.
*         Its iterator type is `InputIterator`.
*
* @param pmesh the polygon mesh to be modified by stitching
* @param hedge_pairs_to_stitch a range of `std::pair` of halfedges to be stitched together
*
*/
template <typename PolygonMesh,
          typename HalfedgePairsRange>
void stitch_borders(PolygonMesh& pmesh,
                    const HalfedgePairsRange& hedge_pairs_to_stitch)
{
  using boost::choose_param;
  using boost::get_param;

  internal::stitch_borders_impl(pmesh, hedge_pairs_to_stitch);
}

/// \ingroup PMP_repairing_grp
/// Same as the other overload but the pairs of halfedges to be stitched
/// are automatically found amongst all border halfedges.
/// Two border halfedges `h1` and `h2` are set to be stitched
/// if the points associated to the source and target vertices of `h1` are
/// the same as those of the target and source vertices of `h2` respectively.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
/// @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
///
/// @param pmesh the polygon mesh to be modified by stitching
/// @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
/// If this parameter is omitted, an internal property map for
/// `CGAL::vertex_point_t` should be available in `PolygonMesh`\cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename PolygonMesh, class CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void stitch_borders(PolygonMesh& pmesh, const CGAL_PMP_NP_CLASS& np)
{
  using boost::choose_param;
  using boost::get_param;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
    halfedge_descriptor;
  std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > hedge_pairs_to_stitch;

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_PMP_NP_CLASS>::const_type VPMap;
  VPMap vpm = choose_param(get_param(np, internal_np::vertex_point),
                           get_const_property_map(vertex_point, pmesh));

  internal::collect_duplicated_stitchable_boundary_edges(pmesh,
    std::back_inserter(hedge_pairs_to_stitch),
    internal::Less_for_halfedge<PolygonMesh, VPMap>(pmesh, vpm), vpm);

  stitch_borders(pmesh, hedge_pairs_to_stitch);
  internal::stitch_boundary_cycle_2(pmesh);
}


///\cond SKIP_IN_MANUAL
template <typename PolygonMesh>
void stitch_borders(PolygonMesh& pmesh)
{
  stitch_borders(pmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

//for backward compatibility
template <typename PolygonMesh,
          typename HalfedgePairsRange,
          typename NamedParameters>
void stitch_borders(PolygonMesh& pmesh,
                    const HalfedgePairsRange& hedge_pairs_to_stitch,
                    NamedParameters)
{
  stitch_borders(pmesh, hedge_pairs_to_stitch);
}
///\endcond

} //end of namespace Polygon_mesh_processing

} //end of namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_STITCH_POLYGON_MESH_H
