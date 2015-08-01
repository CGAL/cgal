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
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_STITCH_POLYGON_MESH_H
#define CGAL_STITCH_POLYGON_MESH_H

#include <CGAL/Modifier_base.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <set>
#include <vector>
#include <utility>
#include <boost/range.hpp>
#include <boost/foreach.hpp>

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
detect_duplicated_boundary_edges
(PM& pmesh, OutputIterator out, LessHedge less_hedge, const VertexPointMap& vpmap)
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef std::set<halfedge_descriptor, LessHedge> Border_halfedge_set;

  Border_halfedge_set border_halfedge_set(less_hedge);
  BOOST_FOREACH(halfedge_descriptor he, halfedges(pmesh))
  {
    if ( !CGAL::is_border(he, pmesh) )
      continue;
    typename Border_halfedge_set::iterator set_it;
    bool insertion_ok;
    CGAL::cpp11::tie(set_it, insertion_ok)
      = border_halfedge_set.insert(he);

    if ( !insertion_ok )
      if ( get(vpmap, source(he,pmesh))==get(vpmap, target(*set_it,pmesh)) &&
           get(vpmap, target(he,pmesh))==get(vpmap, source(*set_it,pmesh)) )
        *out++ = std::make_pair(*set_it, he);
  }
  return out;
}

template <class PM, typename VertexPointMap, typename HalfedgePairsRange>
struct Naive_border_stitching_modifier
  : CGAL::Modifier_base<PM>
{
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename std::pair<halfedge_descriptor, halfedge_descriptor> halfedges_pair;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  typedef HalfedgePairsRange To_stitch;

  Naive_border_stitching_modifier(const To_stitch& to_stitch_,
                                  VertexPointMap vpmap_)
    : to_stitch(to_stitch_)
    , vpmap(vpmap_)
  {}

  void update_target_vertex(halfedge_descriptor h,
                            vertex_descriptor v_kept,
                            PM& pmesh)
  {
    halfedge_descriptor start = h;
    do{
      set_target(h, v_kept, pmesh);
      h = opposite(next(h, pmesh), pmesh);
    } while( h != start );
  }

  void operator()(PM& pmesh)
  {
    /// Merge the vertices
    std::vector<vertex_descriptor> vertices_to_delete;
    // since there might be several vertices with identical point
    // we use the following map to choose one vertex per point
    std::map<Point, vertex_descriptor> vertices_kept;

    BOOST_FOREACH(const halfedges_pair hk, to_stitch)
    {
      halfedge_descriptor h1 = hk.first;
      halfedge_descriptor h2 = hk.second;

      CGAL_assertion(CGAL::is_border(h1, pmesh));
      CGAL_assertion(CGAL::is_border(h2, pmesh));
      CGAL_assertion(!CGAL::is_border(opposite(h1, pmesh), pmesh));
      CGAL_assertion(!CGAL::is_border(opposite(h2, pmesh), pmesh));

      vertex_descriptor h1_tgt = target(h1, pmesh);
      vertex_descriptor h2_src = source(h2, pmesh);

      //update vertex pointers: target of h1 vs source of h2
      vertex_descriptor v_to_keep = h1_tgt;
      std::pair<typename std::map<Point, vertex_descriptor>::iterator, bool >
        insert_res =
        vertices_kept.insert( std::make_pair(get(vpmap, v_to_keep), v_to_keep) );

      if (!insert_res.second && v_to_keep != insert_res.first->second)
      {
        v_to_keep = insert_res.first->second;
        //we remove h1->vertex()
        vertices_to_delete.push_back( h1_tgt );
        update_target_vertex(h1, v_to_keep, pmesh);
      }
      if (v_to_keep != h2_src)
      {
        //we remove h2->opposite()->vertex()
        vertices_to_delete.push_back( h2_src );
        update_target_vertex(opposite(h2, pmesh), v_to_keep, pmesh);
      }
      set_halfedge(v_to_keep, h1, pmesh);

      vertex_descriptor h1_src = source(h1, pmesh);
      vertex_descriptor h2_tgt = target(h2, pmesh);

      //update vertex pointers: target of h1 vs source of h2
      v_to_keep = h2_tgt;
      insert_res =
          vertices_kept.insert( std::make_pair(get(vpmap, v_to_keep), v_to_keep) );
      if (!insert_res.second && v_to_keep != insert_res.first->second)
      {
        v_to_keep = insert_res.first->second;
        //we remove h2->vertex()
        vertices_to_delete.push_back( h2_tgt );
        update_target_vertex(h2, v_to_keep, pmesh);
      }
      if (v_to_keep!=h1_src)
      {
        //we remove h1->opposite()->vertex()
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

private:
  const To_stitch& to_stitch;
  VertexPointMap vpmap;
};

} //end of namespace internal



/*!
* \ingroup stitching_grp
* Stitches together border halfedges in a polygon mesh.
* The halfedges to be stitched are provided in `hedge_pairs_to_stitch`.
* For each pair `p` in this vector, p.second and its opposite will be removed
* from `pmesh`.
*
* The vertices that get removed from `pmesh` are selected as follows:
* The pair of halfedges in `hedge_pairs_to_stitch` are processed linearly.
* Let `p` be such a pair.
* If the target of `p.first` has not been marked for deletion,
* then the source of `p.second` is.
* If the target of `p.second` has not been marked for deletion,
* then the source of `p.first` is.
* A vertex is marked for deletion if its corresponding point has already be seen
* with a vertex previously handled.
*
* \pre For each halfedge in a pair of `hedge_pairs_to_stitch`, the corresponding
*      edge is neither degenerated nor incident to a degenerate edge.
*
* @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
*        that has a property map for `boost::vertex_point_t`
* @tparam HalfedgePairsRange a range of
*         `std::pair<boost::graph_traits<PolygonMesh>::%halfedge_descriptor,
*         boost::graph_traits<PolygonMesh>::%halfedge_descriptor>`,
*         model of `Range`.
*         Its iterator type is `InputIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh the polygon mesh to be modified by stitching
* @param hedge_pairs_to_stitch a range of `std::pair` of halfedges to be stitched together
* @param np optional \ref namedparameters described below
*
* \cgalNamedParamsBegin
*    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
* \cgalNamedParamsEnd
*/
template <typename PolygonMesh,
          typename HalfedgePairsRange,
          typename NamedParameters>
void stitch_borders(PolygonMesh& pmesh,
                    const HalfedgePairsRange& hedge_pairs_to_stitch,
                    const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VPMap;
  VPMap vpm = choose_const_pmap(get_param(np, boost::vertex_point),
                                pmesh,
                                boost::vertex_point);

  internal::Naive_border_stitching_modifier<PolygonMesh, VPMap, HalfedgePairsRange>
    modifier(hedge_pairs_to_stitch, vpm);

  modifier(pmesh);
}

///\cond SKIP_IN_MANUAL
template <typename PolygonMesh, typename HalfedgePairsRange>
void stitch_borders(PolygonMesh& pmesh,
                    const HalfedgePairsRange& hedge_pairs_to_stitch)
{
  stitch_borders(pmesh, hedge_pairs_to_stitch,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}
///\endcond

/// \ingroup stitching_grp
/// Same as the other overload but the pairs of halfedges to be stitched
/// are automatically found amongst all border halfedges.
/// Two border halfedges `h1` and `h2` are set to be stitched
/// if the points associated to the source and target vertices of `h1` are
/// the same as those of the target and source vertices of `h2` respectively.
///
/// \pre `pmesh` does not contains any degenerate border edge.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
///        that has a property map for `boost::vertex_point_t`
/// @tparam NamedParameters a sequence of \ref namedparameters
///
/// @param pmesh the polygon mesh to be modified by stitching
/// @param np optional sequence of \ref namedparameters among the ones listed below
///
/// \cgalNamedParamsBegin
///    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
/// \cgalNamedParamsEnd
///
template <typename PolygonMesh, class CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void stitch_borders(PolygonMesh& pmesh, const CGAL_PMP_NP_CLASS& np)
{
  using boost::choose_param;
  using boost::choose_const_pmap;
  using boost::get_param;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
    halfedge_descriptor;
  std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > hedge_pairs_to_stitch;

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_PMP_NP_CLASS>::const_type VPMap;
  VPMap vpm = choose_const_pmap(get_param(np, boost::vertex_point),
                                pmesh,
                                boost::vertex_point);

  internal::detect_duplicated_boundary_edges(pmesh,
    std::back_inserter(hedge_pairs_to_stitch),
    internal::Less_for_halfedge<PolygonMesh, VPMap>(pmesh, vpm), vpm);

  stitch_borders(pmesh, hedge_pairs_to_stitch, np);
}


///\cond SKIP_IN_MANUAL
template <typename PolygonMesh>
void stitch_borders(PolygonMesh& pmesh)
{
  stitch_borders(pmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}
///\endcond

} //end of namespace Polygon_mesh_processing

} //end of namespace CGAL


#endif //CGAL_STITCH_POLYGON_MESH_H
