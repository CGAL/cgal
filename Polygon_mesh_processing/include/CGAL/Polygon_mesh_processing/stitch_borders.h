// Copyright (c) 2014 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_STITCH_BORDERS_H
#define CGAL_POLYGON_MESH_PROCESSING_STITCH_BORDERS_H

#include <CGAL/license/Polygon_mesh_processing/combinatorial_repair.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Union_find.h>
#include <CGAL/utility.h>
#include <CGAL/use.h>

#include <boost/range.hpp>
#include <boost/functional/hash.hpp>

#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <type_traits>

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
# ifndef CGAL_PMP_STITCHING_DEBUG
#  define CGAL_PMP_STITCHING_DEBUG
# endif
#endif

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

////// Helper structs

// Used to compare halfedges based on their geometry
template <typename PolygonMesh, typename VertexPointMap, typename GeomTraits>
struct Less_for_halfedge
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::reference     Point;

  Less_for_halfedge(const PolygonMesh& pmesh_, const VertexPointMap vpm_, const GeomTraits& gt_)
    : pmesh(pmesh_), vpm(vpm_), gt(gt_)
  {}

  bool operator()(const halfedge_descriptor h1, const halfedge_descriptor h2) const
  {
    typename GeomTraits::Equal_3 equal = gt.equal_3_object();
    typename GeomTraits::Less_xyz_3 less = gt.less_xyz_3_object();

    vertex_descriptor vm1 = source(h1, pmesh);
    vertex_descriptor vM1 = target(h1, pmesh);
    vertex_descriptor vm2 = source(h2, pmesh);
    vertex_descriptor vM2 = target(h2, pmesh);

    if(less(get(vpm, vM1), get(vpm, vm1)))
      std::swap(vM1, vm1);
    if(less(get(vpm, vM2), get(vpm, vm2)))
      std::swap(vM2, vm2);

    Point pm1 = get(vpm, vm1);
    Point pM1 = get(vpm, vM1);
    Point pm2 = get(vpm, vm2);
    Point pM2 = get(vpm, vM2);

    if(equal(pm1, pm2))
      return less(pM1, pM2);

    return less(pm1, pm2);
  }

  const PolygonMesh& pmesh;
  const VertexPointMap vpm;
  const GeomTraits& gt;
};

// The following structs determine which of the two halfedges is kept when a pair is merged
template <typename Mesh>
struct Default_halfedges_keeper
{
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor       halfedge_descriptor;

  halfedge_descriptor operator()(const halfedge_descriptor h1, const halfedge_descriptor h2) const
  {
    return (h1 < h2) ? h1 : h2; // Arbitrary preference
  }
};

// Provided for convenience. If passed to stitch_borders(), this ensures that if any of the two
// edges is marked, then the preserver edge is also marked.
template <typename EMM, typename PolygonMesh>
struct Halfedges_keeper_with_marked_edge_priority
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor       halfedge_descriptor;

  Halfedges_keeper_with_marked_edge_priority(const EMM& emm, const PolygonMesh& pmesh)
    : emm(emm), pmesh(pmesh)
  { }

  halfedge_descriptor operator()(const halfedge_descriptor h1, const halfedge_descriptor h2) const
  {
    // If h2 is marked, then whether h1 is marked or not, we can just return h2
    if(get(emm, edge(h2, pmesh)))
      return h2;

    // If only h1 is marked, return h1;
    // If both or none are marked, it does not matter which we are keeping, so also return h1
    return h1;
  }

private:
  const EMM& emm;
  const PolygonMesh& pmesh;
};

// The following structs are visitors used to maintain a set of cycle representatives (i.e. halfedges)
// when stitching only a subset of the boundary cycles of the mesh
//
// This is for the (default) version, where all cycles are being considered hence
// there is nothing to maintain
template <typename PolygonMesh>
struct Dummy_cycle_rep_maintainer
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_iterator             halfedge_iterator;

  Dummy_cycle_rep_maintainer(const PolygonMesh& pmesh) : m_pmesh(pmesh) { }

  Iterator_range<halfedge_iterator> halfedges_to_consider() const { return halfedges(m_pmesh); }

  std::vector<halfedge_descriptor> cycle_representatives() const
  {
    std::vector<halfedge_descriptor> boundary_cycle_representatives;
    extract_boundary_cycles(m_pmesh, std::back_inserter(boundary_cycle_representatives));

    return boundary_cycle_representatives;
  }

  // Dummies just to fit the API
  void add_representative(const halfedge_descriptor) const { }
  void remove_representative(const halfedge_descriptor) const { }
  void clear_representatives() const { }

  template <typename CycleHalfedgeRange, typename FilteredHalfedgePairsRange, typename VPM>
  void update_representatives(const CycleHalfedgeRange&,
                              const FilteredHalfedgePairsRange&,
                              const VPM) const { }

private:
  const PolygonMesh& m_pmesh;
};

// This is the version used when a specific (sub)range of cycles are being used
template <typename PolygonMesh>
struct Boundary_cycle_rep_maintainer
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor             vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;
  typedef typename std::pair<halfedge_descriptor, halfedge_descriptor>             halfedges_pair;

  typedef CGAL::dynamic_halfedge_property_t<bool>                                 Candidate_tag;
  typedef typename boost::property_map<PolygonMesh, Candidate_tag>::type          Candidate_map;

  Boundary_cycle_rep_maintainer(PolygonMesh& pmesh)
    : m_pmesh(pmesh)
  {
    m_candidate_halfedges = get(Candidate_tag(), pmesh);
  }

public:
  std::vector<halfedge_descriptor> halfedges_to_consider() const
  {
    std::vector<halfedge_descriptor> boundaries;
    for(const halfedge_descriptor bh : m_cycle_reps)
      for(halfedge_descriptor h : CGAL::halfedges_around_face(bh, m_pmesh))
        boundaries.push_back(h);

    // There should be only one representative per cycle
    CGAL_assertion(std::set<halfedge_descriptor>(boundaries.begin(), boundaries.end()).size() == boundaries.size());

    return boundaries;
  }

  const std::set<halfedge_descriptor>& cycle_representatives() const { return m_cycle_reps; }

  void add_representative(const halfedge_descriptor h) { m_cycle_reps.insert(h); }
  void remove_representative(const halfedge_descriptor h) { m_cycle_reps.erase(h); }
  void clear_representatives( ) { m_cycle_reps.clear(); }

  // Pick a single cycle representative for each cycle that appears in 'cycle_halfedges'
  // The representative must not appear in 'filtered_stitchable_halfedges' since this contains
  // halfedges that will not be border (or even exist) after stitching
  //
  // Stitching matching pairs on the same boundary can split a single hole into multiple holes,
  // each with their representative
  template <typename CycleHalfedgeRange,
            typename FilteredHalfedgePairsRange,
            typename VPM>
  void update_representatives(const CycleHalfedgeRange& cycle_halfedges,
                              const FilteredHalfedgePairsRange& filtered_stitchable_halfedges,
                              const VPM vpm)
  {
    typedef typename boost::property_traits<VPM>::reference                 Point_ref;

#ifdef CGAL_PMP_STITCHING_DEBUG
    std::cout << "update_representatives(" << cycle_halfedges.size() << ", "
                                           << filtered_stitchable_halfedges.size() << ")" << std::endl;
#endif

    CGAL_assertion(!cycle_halfedges.empty());

    for(const halfedge_descriptor h : cycle_halfedges)
      put(m_candidate_halfedges, h, true);

    for(const halfedges_pair& hp : filtered_stitchable_halfedges)
    {
      put(m_candidate_halfedges, hp.first, false);
      put(m_candidate_halfedges, hp.second, false);
    }

    for(const halfedge_descriptor h : cycle_halfedges)
    {
      if(!is_border(h, m_pmesh) || !get(m_candidate_halfedges, h))
        continue;

      // This halfedge is now the representative of a boundary cycle
      // --> walk the cycle to unmark halfedges as potential representatives
      // --> when encountering a halfedge that is already not a potential representative,
      //     that means this halfedge is going to be stitched and so the part of the old cycle
      //     starting from that halfedge till we are back to the that vertex should be ignored.
      add_representative(h);

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
      int border_length = 0;
#endif

      halfedge_descriptor walker_h = h;
      for(;;)
      {
        put(m_candidate_halfedges, walker_h, false);

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
        std::cout << "walker at : " << edge(walker_h, m_pmesh) << std::endl;
        ++border_length;
#endif

        // Everything below is just to find the next halfedge of the cycle post-stitching (things
        // aren't stitched yet because we need to keep valid halfedges)
        vertex_descriptor curr_v = target(walker_h, m_pmesh);
        Point_ref curr_p = get(vpm, curr_v);

        walker_h = next(walker_h, m_pmesh);
        CGAL_assertion(is_border(walker_h, m_pmesh));

        if(walker_h == h)
          break;

        bool ignore_till_back_at_curr_p = !get(m_candidate_halfedges, walker_h);
        while(ignore_till_back_at_curr_p) // can have multiple loops at curr_v
        {
#ifdef CGAL_PMP_STITCHING_DEBUG_PP
          std::cout << "Ignoring a cycle starting at " << curr_v << " pos: " << curr_p << std::endl;
#endif

          // walk the cycle to be ignored
          do
          {
#ifdef CGAL_PMP_STITCHING_DEBUG_PP
            std::cout << "Ignoring " << edge(walker_h, m_pmesh)
                      << "(" << source(walker_h, m_pmesh) << " " << target(walker_h, m_pmesh) << ")" << std::endl;
#endif
            CGAL_assertion(walker_h != h);
            walker_h = next(walker_h, m_pmesh);
          }
          while(get(vpm, source(walker_h, m_pmesh)) != curr_p);

          ignore_till_back_at_curr_p = (walker_h != h) && !get(m_candidate_halfedges, walker_h);
        }

        if(walker_h == h)
          break;

        CGAL_assertion(get(m_candidate_halfedges, walker_h));
      }

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
      std::cout << "new cycle rep: " << edge(h, m_pmesh)
                << "\n\t" << source(h, m_pmesh) << "\t(" << get(vpm, source(h, m_pmesh)) << ")"
                << "\n\t" << target(h, m_pmesh) << "\t(" << get(vpm, target(h, m_pmesh)) << ")\n"
                << " length: " << border_length << std::endl;
      CGAL_assertion(border_length >= 2); // length 2 can happen because not everything is stitchable
#endif
    }
  }

private:
  std::set<halfedge_descriptor> m_cycle_reps;
  Candidate_map m_candidate_halfedges; // candidacy to be a representative
  PolygonMesh& m_pmesh;
};

////// Functions

//add a pair of border halfedges to be stitched.
//Specifies if they are manifold or not in the map.
template<typename Halfedge,
         typename Border_halfedge_map,
         typename Halfedge_pair,
         typename Manifold_halfedge_pair,
         typename Mesh,
         typename VPM,
         typename GT>
void fill_pairs(const Halfedge& he,
                Border_halfedge_map& border_halfedge_map,
                Halfedge_pair& halfedge_pairs,
                Manifold_halfedge_pair& manifold_halfedge_pairs,
                const Mesh& pmesh,
                VPM vpm,
                const GT& gt)
{
  typename GT::Equal_3 equal = gt.equal_3_object();

  typename Border_halfedge_map::iterator set_it;
  bool insertion_ok;
  std::tie(set_it, insertion_ok) = border_halfedge_map.emplace(he, std::make_pair(1,0));

  if(!insertion_ok) // there is already a halfedge with the points
  {
    ++set_it->second.first; // increase the multiplicity
    if(set_it->second.first == 2)
    {
      const Halfedge other_he = set_it->first;
      set_it->second.second = halfedge_pairs.size(); // set the id of the pair in the vector
      halfedge_pairs.emplace_back(other_he, he);
      if(equal(get(vpm, source(he,pmesh)), get(vpm, target(other_he, pmesh))) &&
         equal(get(vpm, target(he,pmesh)), get(vpm, source(other_he, pmesh))))
      {
        // Even if the halfedges are compatible, refuse to stitch if that would break the graph
        if(face(opposite(he, pmesh), pmesh) == face(opposite(other_he, pmesh), pmesh))
          manifold_halfedge_pairs.push_back(false);
        else
          manifold_halfedge_pairs.push_back(true);
      }
      else
      {
        manifold_halfedge_pairs.push_back(false);
      }
    }
    else if(set_it->second.first > 2)
    {
      manifold_halfedge_pairs[ set_it->second.second ] = false;
    }
  }
}

template <typename HalfedgeRange,
          typename PolygonMesh,
          typename HalfedgeKeeper,
          typename OutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator collect_duplicated_stitchable_boundary_edges(const HalfedgeRange& halfedge_range,
                                                            PolygonMesh& pmesh,
                                                            const HalfedgeKeeper& hd_kpr,
                                                            const bool per_cc,
                                                            OutputIterator out,
                                                            const CGAL_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor          halfedge_descriptor;

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_NP_CLASS>::const_type  VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef CGAL::dynamic_face_property_t<int>                                      Face_property_tag;
  typedef typename boost::property_map<PolygonMesh, Face_property_tag>::type      Face_cc_map;

  Face_cc_map cc;
  std::size_t num_cc = 0;
  std::vector<std::vector<halfedge_descriptor> > border_edges_per_cc;

  typedef Less_for_halfedge<PolygonMesh, VPM, GT>                                 Less_hedge;
  typedef std::map<halfedge_descriptor, std::pair<int, std::size_t>, Less_hedge>  Border_halfedge_map;

  Less_hedge less_hedge(pmesh, vpm, gt);
  Border_halfedge_map border_halfedge_map(less_hedge);

  std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> > halfedge_pairs;
  std::vector<bool> manifold_halfedge_pairs;

#ifdef CGAL_PMP_STITCHING_DEBUG
  std::cout << "Collecting stitchable pair from a hrange of size: " << halfedge_range.size()
            << " (total: " << halfedges(pmesh).size() << ")" << std::endl;
#endif

  if(per_cc)
  {
    // 'PMP::connected_component' is not local, but it is cheap.
    // One could write a local version of PMP::connected_components with seed faces and which
    // would not use a global dynamic pmap (but rather an unordered set) to mark visited faces,
    // but it would have to be very small CCs for it to yield a better runtime.
    // Consequently, leaving the global version till it is the limiting factor.
    cc = get(Face_property_tag(), pmesh);
    num_cc = connected_components(pmesh, cc, np);
    border_edges_per_cc.resize(num_cc);
  }

  for(halfedge_descriptor he : halfedge_range)
  {
    if(!CGAL::is_border(he, pmesh))
      continue;

    if(per_cc)
      border_edges_per_cc[get(cc, face(opposite(he, pmesh), pmesh))].push_back(he);
    else
      fill_pairs(he, border_halfedge_map, halfedge_pairs, manifold_halfedge_pairs, pmesh, vpm, gt);
  }

  if(per_cc)
  {
    for(std::size_t i=0; i<num_cc; ++i)
    {
      CGAL_assertion(halfedge_pairs.empty());
      CGAL_assertion(manifold_halfedge_pairs.empty());

      Border_halfedge_map border_halfedge_map_in_cc(less_hedge);
      for(std::size_t j=0; j<border_edges_per_cc[i].size(); ++j)
      {
        halfedge_descriptor he = border_edges_per_cc[i][j];
        fill_pairs(he, border_halfedge_map_in_cc, halfedge_pairs,
                   manifold_halfedge_pairs, pmesh, vpm, gt);
      }

      // put in `out` only manifold edges from the set of edges to stitch.
      // We choose not to allow only a pair out of the whole set to be stitched
      // as we can produce inconsistent stitching along a sequence of non-manifold edges
      const std::size_t nb_pairs = halfedge_pairs.size();
      for(std::size_t k=0; k<nb_pairs; ++k)
      {
        if(manifold_halfedge_pairs[k])
        {
          // the first halfedge of the pair is kept
          if(hd_kpr(halfedge_pairs[k].first, halfedge_pairs[k].second) == halfedge_pairs[k].second)
            std::swap(halfedge_pairs[k].first, halfedge_pairs[k].second);

          *out++ = halfedge_pairs[k];

#ifdef CGAL_PMP_STITCHING_DEBUG
          const halfedge_descriptor h = halfedge_pairs[k].first;
          const halfedge_descriptor hn = halfedge_pairs[k].second;
          std::cout << "Stitch "
                    << edge(h, pmesh) << " (" << get(vpm, source(h, pmesh)) << ") - (" << get(vpm, target(h, pmesh)) << ") and "
                    << edge(hn, pmesh) << " (" << get(vpm, source(hn, pmesh)) << ") - (" << get(vpm, target(hn, pmesh)) << ")" << std::endl;
#endif
        }
      }

      halfedge_pairs.clear();
      manifold_halfedge_pairs.clear();
    }
  }
  else
  {
    // put in `out` only manifold edges from the set of edges to stitch.
    // We choose not to allow only a pair out of the whole set to be stitched
    // as we can produce inconsistent stitching along a sequence of non-manifold edges
    const std::size_t nb_pairs=halfedge_pairs.size();
    for(std::size_t i=0; i<nb_pairs; ++i)
    {
      if(manifold_halfedge_pairs[i])
      {
        // the first halfedge of the pair is kept
        if(hd_kpr(halfedge_pairs[i].first, halfedge_pairs[i].second) == halfedge_pairs[i].second)
          std::swap(halfedge_pairs[i].first, halfedge_pairs[i].second);

        *out++ = halfedge_pairs[i];

#ifdef CGAL_PMP_STITCHING_DEBUG
        const halfedge_descriptor h = halfedge_pairs[i].first;
        const halfedge_descriptor hn = halfedge_pairs[i].second;
        std::cout << "Stitch "
                  << edge(h, pmesh) << " (" << get(vpm, source(h, pmesh)) << ") - (" << get(vpm, target(h, pmesh)) << ") and "
                  << edge(hn, pmesh) << " (" << get(vpm, source(hn, pmesh)) << ") - (" << get(vpm, target(hn, pmesh)) << ")" << std::endl;
#endif
      }
    }
  }

  return out;
}

template <class PolygonMesh>
void update_target_vertex(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                          typename boost::graph_traits<PolygonMesh>::vertex_descriptor v_kept,
                          PolygonMesh& pmesh)
{
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor start = h;
  do
  {
    set_target(h, v_kept, pmesh);
    h = opposite(next(h, pmesh), pmesh);
  }
  while(h != start);
}

template <class vertex_descriptor, class Handle_map>
typename Union_find<vertex_descriptor>::handle
uf_get_handle(vertex_descriptor v,
              Union_find<vertex_descriptor>& uf_vertices,
              Handle_map& handles)
{
  std::pair<typename Handle_map::iterator, bool> insert_res =
    handles.emplace(v, typename Union_find<vertex_descriptor>::handle());
  if(insert_res.second)
    insert_res.first->second = uf_vertices.make_set(v);

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
template <typename PolygonMesh, typename HalfedgePairsRange, typename VertexPointMap,
          typename Uf_vertices, typename Uf_handles>
void run_stitch_borders(PolygonMesh& pmesh,
                        const HalfedgePairsRange& to_stitch,
                        const VertexPointMap& vpm,
                        Uf_vertices& uf_vertices,
                        Uf_handles& uf_handles)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor    halfedge_descriptor;
  typedef typename std::pair<halfedge_descriptor, halfedge_descriptor>      halfedges_pair;

  CGAL_USE(vpm);

  std::vector<vertex_descriptor> vertices_to_delete;
  for(const halfedges_pair& hk : to_stitch)
  {
    halfedge_descriptor h1 = hk.first;
    halfedge_descriptor h2 = hk.second;

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
      std::cout << "Actually stitching:\n";
      std::cout << edge(h1, pmesh) << "\n\t" << source(h1, pmesh) << "\t(" << get(vpm, source(h1, pmesh)) << ")"
                                   << "\n\t" << target(h1, pmesh) << "\t(" << get(vpm, target(h1, pmesh)) << ")\n"
                << edge(h2, pmesh) << "\n\t" << source(h2, pmesh) << "\t(" << get(vpm, source(h2, pmesh)) << ")"
                                   << "\n\t" << target(h2, pmesh) << "\t(" << get(vpm, target(h2, pmesh)) << ")" << std::endl;
#endif

    vertex_descriptor h1_tgt = target(h1, pmesh);
    vertex_descriptor h2_src = source(h2, pmesh);

    // update vertex pointers: target of h1 vs source of h2
    vertex_descriptor v_to_keep =
      *uf_vertices.find(uf_get_handle(h1_tgt, uf_vertices, uf_handles));

    if(v_to_keep!=h1_tgt)
    {
      vertices_to_delete.push_back(h1_tgt);
      update_target_vertex(h1, v_to_keep, pmesh);
    }

    if(v_to_keep != h2_src && h1_tgt!=h2_src)
    {
      vertices_to_delete.push_back( h2_src );
      update_target_vertex(opposite(h2, pmesh), v_to_keep, pmesh);
    }
    set_halfedge(v_to_keep, h1, pmesh);

    vertex_descriptor h1_src = source(h1, pmesh);
    vertex_descriptor h2_tgt = target(h2, pmesh);

    //update vertex pointers: target of h1 vs source of h2
    v_to_keep = *uf_vertices.find(uf_get_handle(h2_tgt, uf_vertices, uf_handles));
    if(v_to_keep!=h2_tgt)
    {
      vertices_to_delete.push_back( h2_tgt );
      update_target_vertex(h2, v_to_keep, pmesh);
    }

    if(v_to_keep!=h1_src && h1_src!=h2_tgt)
    {
      vertices_to_delete.push_back( h1_src );
      update_target_vertex(opposite(h1, pmesh), v_to_keep, pmesh);
    }

    set_halfedge(v_to_keep, opposite(h1, pmesh), pmesh);
  }

  /// Update next/prev of neighbor halfedges (that are not set for stiching)
  /// _______   _______
  ///        | |
  ///        | |
  /// In order to avoid having to maintain a set with halfedges to stitch
  /// we do on purpose next-prev linking that might not be useful but that
  /// is harmless and still less expensive than doing queries in a set
  for(const halfedges_pair& hk : to_stitch)
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
  for(const halfedges_pair& hk : to_stitch)
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
  for(vertex_descriptor vd : vertices_to_delete)
  {
#ifdef CGAL_PMP_STITCHING_DEBUG_PP
    std::cout << "Delete vertex: " << vd << std::endl;
#endif
    remove_vertex(vd, pmesh);
  }
}

template <typename PolygonMesh, typename HalfedgePair,
          typename Uf_vertices, typename Uf_handles>
const std::vector<HalfedgePair>&
filter_stitchable_pairs(PolygonMesh& pmesh,
                        const std::vector<HalfedgePair>& to_stitch,
                        std::vector<HalfedgePair>& to_stitch_filtered,
                        Uf_vertices& uf_vertices,
                        Uf_handles& uf_handles)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;

  // Merge the vertices
  for(const HalfedgePair& hk : to_stitch)
  {
    const halfedge_descriptor h1 = hk.first;
    const halfedge_descriptor h2 = hk.second;

    CGAL_assertion(is_border(h1, pmesh));
    CGAL_assertion(is_border(h2, pmesh));
    CGAL_assertion(!is_border(opposite(h1, pmesh), pmesh));
    CGAL_assertion(!is_border(opposite(h2, pmesh), pmesh));

    vertex_descriptor tgt1 = target(h1, pmesh), src1 = source(h1, pmesh);
    vertex_descriptor src2 = source(h2, pmesh), tgt2 = target(h2, pmesh);
    uf_join_vertices(tgt1, src2, uf_vertices, uf_handles);
    uf_join_vertices(src1, tgt2, uf_vertices, uf_handles);
  }

  // detect vertices that cannot be stitched because it would produce a non-manifold edge
  // We look for vertex to be stitched and collect all incident edges with another endpoint
  // to be stitched (that is not an edge scheduled for stitching). That way we can detect
  // if more that one edge will share the same two "master" endpoints.
  typedef std::pair<vertex_descriptor, vertex_descriptor> Vertex_pair;
  typedef std::unordered_map<Vertex_pair,
                             std::vector<halfedge_descriptor>,
                             boost::hash<Vertex_pair>>           Halfedges_after_stitching;
  Halfedges_after_stitching halfedges_after_stitching;

  typedef std::pair<const vertex_descriptor, typename Uf_vertices::handle> Pair_type;
  for(const Pair_type& p : uf_handles)
  {
    vertex_descriptor vd = p.first;
    typename Uf_vertices::handle tgt_handle = uf_vertices.find(uf_handles[vd]);
    for(halfedge_descriptor hd : halfedges_around_target(vd, pmesh))
    {
      const vertex_descriptor other_vd = source(hd, pmesh);

      typename Uf_handles::iterator it_res = uf_handles.find(other_vd);

      if(it_res!=uf_handles.end()) // if the other vertex is also involved in a merge
      {
        if(other_vd < vd)
          continue; // avoid reporting twice the same edge

        typename Uf_vertices::handle src_handle=uf_vertices.find(it_res->second);
        halfedges_after_stitching[make_sorted_pair(*tgt_handle, *src_handle)].push_back(hd);
      }
      else
        halfedges_after_stitching[make_sorted_pair(*tgt_handle, other_vd)].push_back(hd);
    }
  }

  // look for edges that will be present more than once after the stitching
  // (no edges scheduled for stitching involved)
  std::unordered_set<vertex_descriptor> unstitchable_vertices;
  for(typename Halfedges_after_stitching::iterator it=halfedges_after_stitching.begin(),
                                                   it_end=halfedges_after_stitching.end();
                                                   it!=it_end; ++it)
  {
    switch(it->second.size())
    {
      case 1:
       break; // nothing to do
      case 2:
      {
        if(is_border_edge(it->second.front(), pmesh) && is_border_edge(it->second.back(), pmesh))
          break; // these are edges that are most possibly scheduled for stitching or will create a two halfedge loop
        CGAL_FALLTHROUGH;
      }
      default:
      {
        // this is a bit extreme as maybe some could be stitched
        // (but safer because the master could be one of them)
        for(halfedge_descriptor hd : it->second)
        {
          unstitchable_vertices.insert(source(hd, pmesh));
          unstitchable_vertices.insert(target(hd, pmesh));
        }
      }
    }
  }

  // filter halfedges to stitch
  if(!unstitchable_vertices.empty())
  {
    to_stitch_filtered.reserve(to_stitch.size());
    for(const HalfedgePair& hk : to_stitch)
    {
      // We test both halfedges because the previous test
      // might involve only one of the two halfedges
      if(unstitchable_vertices.count( source(hk.first, pmesh) )== 0 &&
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
    for(const HalfedgePair& hk : to_stitch_filtered)
    {
      halfedge_descriptor h1 = hk.first;
      halfedge_descriptor h2 = hk.second;

      vertex_descriptor tgt1 = target(h1, pmesh), src1 = source(h1, pmesh);
      vertex_descriptor src2 = source(h2, pmesh), tgt2 = target(h2, pmesh);
      uf_join_vertices(tgt1, src2, uf_vertices, uf_handles);
      uf_join_vertices(src1, tgt2, uf_vertices, uf_handles);
    }
    return to_stitch_filtered;
  }
  return to_stitch;
}

template <typename HalfedgePair, typename CandidateHalfedgeRange, typename PolygonMesh,
          typename CycleRepMaintainer, typename VertexPointMap>
std::size_t stitch_halfedge_range(const std::vector<HalfedgePair>& to_stitch,
                                  const CandidateHalfedgeRange& representative_candidates,
                                  PolygonMesh& pmesh,
                                  const VertexPointMap& vpm,
                                  CycleRepMaintainer& cycle_reps_maintainer)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  std::cout << to_stitch.size() << " tentative halfedge pair(s) to stitch:\n";
  for(const auto& hp : to_stitch)
  {
    const halfedge_descriptor h = hp.first;
    const halfedge_descriptor hn = hp.second;

    std::cout << edge(h, pmesh) << "\n\t" << source(h, pmesh) << "\t(" << get(vpm, source(h, pmesh)) << ")"
                                << "\n\t" << target(h, pmesh) << "\t(" << get(vpm, target(h, pmesh)) << ")\n"
              << edge(hn, pmesh) << "\n\t" << source(hn, pmesh) << "\t(" << get(vpm, source(hn, pmesh)) << ")"
                                 << "\n\t" << target(hn, pmesh) << "\t(" << get(vpm, target(hn, pmesh)) << ")" << std::endl;
  }
#endif

  // The first step of the algorithm is to filter halfedges to be stitched so that
  // after stitching no edges will be present more than once.
  typedef CGAL::Union_find<vertex_descriptor> Uf_vertices;
  Uf_vertices uf_vertices;
  typedef std::unordered_map<vertex_descriptor, typename Uf_vertices::handle> Uf_handles;
  Uf_handles uf_handles;

  std::vector<HalfedgePair> to_stitch_local;
  const std::vector<HalfedgePair>& to_stitch_filtered =
    filter_stitchable_pairs(pmesh, to_stitch, to_stitch_local, uf_vertices, uf_handles);

  cycle_reps_maintainer.update_representatives(representative_candidates, to_stitch_filtered, vpm);

  // Actually stitching
  run_stitch_borders(pmesh, to_stitch_filtered, vpm, uf_vertices, uf_handles);

  return to_stitch_filtered.size();
}

template <typename HalfedgePair, typename PolygonMesh, typename VertexPointMap>
std::size_t stitch_halfedge_range(const std::vector<HalfedgePair>& to_stitch,
                                  PolygonMesh& pmesh,
                                  const VertexPointMap& vpm)
{
  Dummy_cycle_rep_maintainer<PolygonMesh> cycle_reps_maintainer(pmesh);
  return stitch_halfedge_range(to_stitch, halfedges(pmesh), pmesh, vpm, cycle_reps_maintainer);
}

// overload to avoid a useless copy
template <typename HalfedgePair, typename PolygonMesh, typename VertexPointMap>
std::size_t stitch_halfedge_range_dispatcher(const std::vector<HalfedgePair>& to_stitch,
                                             PolygonMesh& pmesh,
                                             const VertexPointMap& vpm)
{
  return stitch_halfedge_range(to_stitch, pmesh, vpm);
}

// overload making a copy
template <typename HalfedgePairRange, typename PolygonMesh, typename VertexPointMap>
std::size_t stitch_halfedge_range_dispatcher(const HalfedgePairRange& to_stitch_const,
                                             PolygonMesh& pmesh,
                                             const VertexPointMap& vpm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor HD;
  std::vector< std::pair<HD,HD> > to_stitch(to_stitch_const.begin(), to_stitch_const.end());
  return stitch_halfedge_range(to_stitch, pmesh, vpm);
}

// collect_duplicated_stitchable_boundary_edges() cannot handle configurations with non-manifoldness.
// However, even if non-manifoldness exists within a loop, it is safe choice to stitch consecutive
// stitchable halfedges
template <typename HalfedgeRange,
          typename HalfedgeKeeper,
          typename PolygonMesh,
          typename VPM,
          typename GT>
std::size_t zip_boundary_cycle(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor& bh,
                               const HalfedgeRange& cycle_halfedges,
                               const HalfedgeKeeper& hd_kpr,
                               PolygonMesh& pmesh,
                               const VPM vpm,
                               const GT& gt)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;

  typename GT::Equal_3 equal = gt.equal_3_object();

  std::size_t stitched_boundary_cycles_n = 0;

  // Zipping cannot change the topology of the hole so the maintenance is trivial
  internal::Dummy_cycle_rep_maintainer<PolygonMesh> dummy_maintainer(pmesh);

  // A boundary cycle might need to be stitched starting from different extremities
  //
  //                        v11 ------ v10
  //                         |          |
  //   v0 --- v1(v13) === v2(v12)     v5(v9) === v6(v8) --- v7
  //                         |          |
  //                        v3 ------- v4
  //
  // As long as we find vertices on the boundary with both incident halfedges being compatible,
  // we zip it up as much as possible.

  // not everything is always stitchable
  std::set<halfedge_descriptor> unstitchable_halfedges;

  const halfedge_descriptor null_h = boost::graph_traits<PolygonMesh>::null_halfedge();
  for(;;) // until there is nothing to stitch anymore
  {
    if(bh == null_h) // the complete boundary cycle is stitched
      break;

#ifdef CGAL_PMP_STITCHING_DEBUG
    std::cout << "Walking border from halfedge: " << edge(bh, pmesh) << std::endl;
#endif

    CGAL_assertion(is_border(bh, pmesh));

    halfedge_descriptor hn = next(bh, pmesh), start_h = null_h;
    do
    {
      halfedge_descriptor hnn = next(hn, pmesh);
      CGAL_assertion(equal(get(vpm, target(hn, pmesh)), get(vpm, source(hnn, pmesh))));

      if(equal(get(vpm, source(hn, pmesh)), get(vpm, target(hnn, pmesh))) &&
         !is_degenerate_edge(edge(hn, pmesh), pmesh,
                             parameters::vertex_point_map(vpm).geom_traits(gt)))
      {
        if(unstitchable_halfedges.count(hn) == 0)
        {
          start_h = hn;
          break;
        }
      }

      hn = hnn;
    }
    while(hn != bh);

    if(start_h == null_h) // nothing to be stitched on this boundary cycle
      break;

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
    std::cout << "Starting stitching from halfedge: "
              << get(vpm, source(edge(start_h, pmesh), pmesh)) << " "
              << get(vpm, target(edge(start_h, pmesh), pmesh)) << std::endl;
#endif

    CGAL_assertion(is_border(start_h, pmesh));

    // Associate as many consecutive halfedge pairs as possible ("zipping")
    std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> > hedges_to_stitch;

    halfedge_descriptor curr_h = start_h;
    halfedge_descriptor curr_hn = next(curr_h, pmesh);
    for(;;) // while we can expand the zipping range
    {
      // Don't create an invalid polygon mesh, even if the geometry allows it
      if(face(opposite(curr_h, pmesh), pmesh) == face(opposite(curr_hn, pmesh), pmesh))
      {
        unstitchable_halfedges.insert(curr_h);
        bh = curr_hn;
        break;
      }

      CGAL_assertion(is_border(curr_h, pmesh));
      CGAL_assertion(is_border(curr_hn, pmesh));

      if(hd_kpr(curr_h, curr_hn) == curr_h)
        hedges_to_stitch.emplace_back(curr_h, curr_hn);
      else
        hedges_to_stitch.emplace_back(curr_hn, curr_h);

#ifdef CGAL_PMP_STITCHING_DEBUG_PP
      std::cout << "expand zip with:\n"
                << edge(curr_h, pmesh) << "\n\t" << source(curr_h, pmesh) << "\t(" << get(vpm, source(curr_h, pmesh)) << ")"
                                    << "\n\t" << target(curr_h, pmesh) << "\t(" << get(vpm, target(curr_h, pmesh)) << ")\n"
                << edge(curr_hn, pmesh) << "\n\t" << source(curr_hn, pmesh) << "\t(" << get(vpm, source(curr_hn, pmesh)) << ")"
                                     << "\n\t" << target(curr_hn, pmesh) << "\t(" << get(vpm, target(curr_hn, pmesh)) << ")" << std::endl;
#endif

      // check if we have reached the end of the boundary cycle
      if(prev(curr_h, pmesh) == curr_hn || prev(curr_h, pmesh) == next(curr_hn, pmesh))
      {
        bh = null_h;
        break;
      }

      curr_h = prev(curr_h, pmesh);
      curr_hn = next(curr_hn, pmesh);

      // check if the next two halfedges are not geometrically compatible
      if(!equal(get(vpm, source(curr_h, pmesh)), get(vpm, target(curr_hn, pmesh))) ||
         is_degenerate_edge(edge(curr_hn, pmesh), pmesh,
                            parameters::vertex_point_map(vpm).geom_traits(gt)))
      {
        bh = curr_hn;
        break;
      }
    }

    // bh must be a boundary halfedge on the border that will not be impacted by any stitching
    CGAL_assertion_code(if(bh != null_h) {)
    CGAL_assertion_code(  for(const auto& hp : hedges_to_stitch) {)
    CGAL_assertion(         bh != hp.first && bh != hp.second);
    CGAL_assertion_code(}})

    if(!hedges_to_stitch.empty())
    {
#ifdef CGAL_PMP_STITCHING_DEBUG_PP
      std::cout << hedges_to_stitch.size() " halfedge pairs to stitch on border containing:\n"
                << edge(h, pmesh) << "\n\t" << source(h, pmesh) << "\t(" << get(vpm, source(h, pmesh)) << ")"
                                  << "\n\t" << target(h, pmesh) << "\t(" << get(vpm, target(h, pmesh)) << ")" << std::endl;
#endif

      std::size_t local_stitches = internal::stitch_halfedge_range(hedges_to_stitch, cycle_halfedges,
                                                                   pmesh, vpm, dummy_maintainer);
      stitched_boundary_cycles_n += local_stitches;

      if(local_stitches == 0) // refused to stitch this halfedge pair range due to manifold issue
      {
#ifdef CGAL_PMP_STITCHING_DEBUG_PP
        std::cout << "Failed to stitch this range!" << std::endl;
#endif

        for(const auto& hp : hedges_to_stitch)
        {
          unstitchable_halfedges.insert(hp.first);
          unstitchable_halfedges.insert(hp.second);
        }
      }
    }
  }

  return stitched_boundary_cycles_n;
}

/// High-level functions

template <typename PolygonMesh, typename CycleRepMaintainer, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_boundary_cycle(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                  PolygonMesh& pmesh,
                                  CycleRepMaintainer& cycle_reps_maintainer,
                                  const CGAL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;
  typedef typename std::pair<halfedge_descriptor, halfedge_descriptor>             halfedges_pair;

  CGAL_precondition(is_valid_halfedge_descriptor(h, pmesh));
  CGAL_precondition(is_border(h, pmesh));
  CGAL_precondition(is_valid(pmesh));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_NP_CLASS>::const_type   VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  typedef typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename internal_np::Lookup_named_param_def<internal_np::halfedges_keeper_t,
                                                       CGAL_NP_CLASS,
                                                       Default_halfedges_keeper<PolygonMesh> >::type  Halfedge_keeper;
  const Halfedge_keeper hd_kpr = choose_parameter(get_parameter(np, internal_np::halfedges_keeper),
                                                  Default_halfedges_keeper<PolygonMesh>());

  halfedge_descriptor bh = h, bh_mem = bh;

  std::vector<halfedge_descriptor> cycle_halfedges;
  for(halfedge_descriptor h : halfedges_around_face(bh, pmesh))
    cycle_halfedges.push_back(h);

  std::size_t res = internal::zip_boundary_cycle(bh, cycle_halfedges, hd_kpr, pmesh, vpm, gt);
  if(bh == boost::graph_traits<PolygonMesh>::null_halfedge()) // stitched everything
  {
    cycle_reps_maintainer.remove_representative(bh);
    return res;
  }

  // Re-compute the range if something was stitched
  if(res != 0)
  {
    cycle_reps_maintainer.remove_representative(bh_mem);
    cycle_reps_maintainer.add_representative(bh);

    cycle_halfedges.clear();
    for(halfedge_descriptor h : halfedges_around_face(bh, pmesh))
      cycle_halfedges.push_back(h);
  }

  std::vector<halfedges_pair> to_stitch;
  internal::collect_duplicated_stitchable_boundary_edges(cycle_halfedges, pmesh,
                                                         hd_kpr, false /*per cc*/,
                                                         std::back_inserter(to_stitch), np);

  res += stitch_halfedge_range(to_stitch, cycle_halfedges, pmesh, vpm, cycle_reps_maintainer);

  return res;
}

} //end of namespace internal

/// \ingroup PMP_combinatorial_repair_grp
///
/// \brief stitches together, whenever possible, two halfedges belonging to the boundary cycle
/// described by the halfedge `h`.
///
/// Two border halfedges `h1` and `h2` can be stitched
/// if the points associated to the source and target vertices of `h1` are
/// the same as those of the target and source vertices of `h2`, respectively.
///
/// \tparam PolygonMesh a model of `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param h a border halfedge of the polygon mesh `pmesh`
/// \param pmesh the polygon mesh to be stitched
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pm`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `PolygonMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functors:
///                    - `Less_xyz_3` to compare lexicographically two points
///                    - `Equal_3` to check whether two points are identical.
///                    For each functor `Foo`, a function `Foo foo_object()` must be provided.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns the number of pairs of halfedges that were stitched.
///
/// \sa `stitch_boundary_cycles()`
/// \sa `stitch_borders()`
///
template <typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_boundary_cycle(const typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                                  PolygonMesh& pmesh,
                                  const CGAL_NP_CLASS& np = parameters::default_values())
{
  internal::Dummy_cycle_rep_maintainer<PolygonMesh> dummy_maintainer(pmesh);
  return internal::stitch_boundary_cycle(h, pmesh, dummy_maintainer, np);
}

namespace internal {

template <typename BorderHalfedgeRange, typename PolygonMesh,
          typename CycleRepMaintainer, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_boundary_cycles(const BorderHalfedgeRange& boundary_cycle_representatives,
                                   PolygonMesh& pmesh,
                                   CycleRepMaintainer& cycle_reps_maintainer,
                                   const CGAL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;

  std::size_t stitched_boundary_cycles_n = 0;
  for(const halfedge_descriptor h : boundary_cycle_representatives)
    stitched_boundary_cycles_n += stitch_boundary_cycle(h, pmesh, cycle_reps_maintainer, np);

  return stitched_boundary_cycles_n;
}

} // namespace internal

/// \ingroup PMP_combinatorial_repair_grp
///
/// \brief stitches together, whenever possible, two halfedges belonging to the same boundary cycle.
///
/// Two border halfedges `h1` and `h2` can be stitched
/// if the points associated to the source and target vertices of `h1` are
/// the same as those of the target and source vertices of `h2`, respectively.
///
/// \tparam BorderHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
/// \tparam PolygonMesh a model of `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param boundary_cycle_representatives a range of border halfedges, each describing a boundary cycle of the mesh `pmesh`
/// \param pmesh the polygon mesh to be modified by stitching
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pm`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `PolygonMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functors:
///                    - `Less_xyz_3` to compare lexicographically two points
///                    - `Equal_3` to check whether two points are identical.
///                    For each functor `Foo`, a function `Foo foo_object()` must be provided.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \returns the number of pairs of halfedges that were stitched.
///
/// \sa `stitch_boundary_cycle()`
/// \sa `stitch_borders()`
///
template <typename BorderHalfedgeRange, typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_boundary_cycles(const BorderHalfedgeRange& boundary_cycle_representatives,
                                   PolygonMesh& pmesh,
                                   const CGAL_NP_CLASS& np = parameters::default_values())
{
  // If this API is called, we are not from stitch_borders() (otherwise there would be a maintainer)
  // so there is only one pass and we don't carea bout maintaining the cycle subset
  internal::Dummy_cycle_rep_maintainer<PolygonMesh> dummy_maintainer(pmesh);
  return stitch_boundary_cycles(boundary_cycle_representatives, pmesh, dummy_maintainer, np);
}

///\cond SKIP_IN_MANUAL

// convenience overloads
template <typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_boundary_cycles(PolygonMesh& pmesh,
                                   const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;

  std::vector<halfedge_descriptor> boundary_cycle_representatives;
  extract_boundary_cycles(pmesh, std::back_inserter(boundary_cycle_representatives));

  return stitch_boundary_cycles(boundary_cycle_representatives, pmesh, np);
}

/// \endcond

// The VPM is only used here for debugging info purposes as in this overload, the halfedges
// to stitch are already provided and all further checks are combinatorial and not geometrical.
/*!
* \ingroup PMP_combinatorial_repair_grp
*
* \brief stitches together border halfedges in a polygon mesh.
*
* The halfedges to be stitched are provided in `hedge_pairs_to_stitch`.
* For each pair `p` in this vector, `p.second` and its opposite will be removed
* from `pmesh`.
*
* \tparam PolygonMesh a model of `MutableFaceGraph`
* \tparam HalfedgePairsRange a range of
*         `std::pair<boost::graph_traits<PolygonMesh>::%halfedge_descriptor,
*         boost::graph_traits<PolygonMesh>::%halfedge_descriptor>`,
*         model of `Range`.
*         Its iterator type is `InputIterator`.
*
* \param pmesh the polygon mesh to be modified by stitching
* \param hedge_pairs_to_stitch a range of `std::pair` of halfedges to be stitched together
* \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pm`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \return the number of pairs of halfedges that were stitched.
*
* \sa `stitch_boundary_cycle()`
* \sa `stitch_boundary_cycles()`
*/
template <typename PolygonMesh,
          typename HalfedgePairsRange,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_borders(PolygonMesh& pmesh,
                           const HalfedgePairsRange& hedge_pairs_to_stitch,
                           const CGAL_NP_CLASS& np = parameters::default_values(),
                           std::enable_if_t<
                             boost::has_range_iterator<HalfedgePairsRange>::value
                           >* = 0)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_NP_CLASS>::const_type  VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  return internal::stitch_halfedge_range_dispatcher(hedge_pairs_to_stitch, pmesh, vpm);
}

namespace internal {

template <typename BorderHalfedgeRange, typename PolygonMesh,
          typename CycleRepMaintainer,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_borders(const BorderHalfedgeRange& boundary_cycle_representatives,
                           PolygonMesh& pmesh,
                           CycleRepMaintainer& cycle_maintainer,
                           const CGAL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  if(boundary_cycle_representatives.size() == 0)
    return 0;

  typedef typename GetVertexPointMap<PolygonMesh, CGAL_NP_CLASS>::const_type  VPM;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pmesh));

  typedef typename internal_np::Lookup_named_param_def<internal_np::halfedges_keeper_t,
                                                       CGAL_NP_CLASS,
                                                       Default_halfedges_keeper<PolygonMesh> >::type  Halfedge_keeper;
  const Halfedge_keeper hd_kpr = choose_parameter(get_parameter(np, internal_np::halfedges_keeper),
                                                  Default_halfedges_keeper<PolygonMesh>());

  bool per_cc = choose_parameter(get_parameter(np, internal_np::apply_per_connected_component), false);

#ifdef CGAL_PMP_STITCHING_DEBUG
  std::cout << "------- Stitch cycles (#1)... (" << boundary_cycle_representatives.size() << " cycle(s))" << std::endl;
#endif

  std::size_t res = stitch_boundary_cycles(boundary_cycle_representatives, pmesh, cycle_maintainer, np);

#ifdef CGAL_PMP_STITCHING_DEBUG
  std::cout << "------- Stitched " << res << " halfedge pairs in boundary cycles" << std::endl;
  std::cout << "------- Stitch all..." << std::endl;
#endif

  const auto& to_consider = cycle_maintainer.halfedges_to_consider();
  cycle_maintainer.clear_representatives();

  std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> > to_stitch;
  internal::collect_duplicated_stitchable_boundary_edges(to_consider, pmesh, hd_kpr, per_cc,
                                                         std::back_inserter(to_stitch), np);
  res += stitch_halfedge_range(to_stitch, to_consider, pmesh, vpm, cycle_maintainer);

#ifdef CGAL_PMP_STITCHING_DEBUG
  std::cout << "------- Stitched " << res << " halfedge pairs after cycles & general" << std::endl;
  std::cout << "------- Stitch cycles (#2)... (" << new_representatives.size() << " cycle(s))" << std::endl;
#endif

  const auto& new_representatives = cycle_maintainer.cycle_representatives();

  // Don't care about keeping track of the sub-cycles as this is the last pass
  internal::Dummy_cycle_rep_maintainer<PolygonMesh> dummy_cycle_maintainer(pmesh);
  res += stitch_boundary_cycles(new_representatives, pmesh, dummy_cycle_maintainer, np);

#ifdef CGAL_PMP_STITCHING_DEBUG
  std::cout << "------- Stitched " << res << " (total)" << std::endl;
#endif

  return res;
}

} // namespace internal

/// \ingroup PMP_combinatorial_repair_grp
///
/// \brief Same as the other overload, but the pairs of halfedges to be stitched
/// are automatically found amongst all border halfedges.
///
/// Two border halfedges `h1` and `h2` are set to be stitched
/// if the points associated to the source and target vertices of `h1` are
/// the same as those of the target and source vertices of `h2`, respectively.
///
/// \tparam BorderHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
/// \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param pmesh the polygon mesh to be modified by the stitching procedure
/// \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `PolygonMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{apply_per_connected_component}
///     \cgalParamDescription{specifies if the borders should only be stitched only within their own connected component.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{face_index_map}
///     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
///                    as key type and `std::size_t` as value type}
///     \cgalParamDefault{an automatically indexed internal map}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return the number of pairs of halfedges that were stitched.
///
/// \sa `stitch_boundary_cycle()`
/// \sa `stitch_boundary_cycles()`
///
template <typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_borders(PolygonMesh& pmesh,
                           const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor           halfedge_descriptor;

  std::vector<halfedge_descriptor> boundary_cycle_representatives;
  extract_boundary_cycles(pmesh, std::back_inserter(boundary_cycle_representatives));

  // We are working on all boundary cycles, so there is no need to keep track of any subset
  internal::Dummy_cycle_rep_maintainer<PolygonMesh> dummy_maintainer(pmesh);
  return stitch_borders(boundary_cycle_representatives, pmesh, dummy_maintainer, np);
}

/// \ingroup PMP_combinatorial_repair_grp
///
/// \brief Same as the other overload, but the pairs of halfedges to be stitched
/// are automatically found amongst halfedges in cycles described by `boundary_cycle_representatives`.
///
/// Two border halfedges `h1` and `h2` are set to be stitched
/// if the points associated to the source and target vertices of `h1` are
/// the same as those of the target and source vertices of `h2`, respectively.
///
/// \tparam BorderHalfedgeRange a model of `Range` with value type `boost::graph_traits<PolygonMesh>::%halfedge_descriptor`
/// \tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param boundary_cycle_representatives a range of border halfedges, each describing a boundary cycle whose halfedges
///                                       will be considered for stitching
/// \param pmesh the polygon mesh to be modified by the stitching procedure
/// \param np optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{apply_per_connected_component}
///     \cgalParamDescription{specifies if the borders should only be stitched only within their own connected component.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{face_index_map}
///     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
///                    as key type and `std::size_t` as value type}
///     \cgalParamDefault{an automatically indexed internal map}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `PolygonMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functors:
///                    - `Less_xyz_3` to compare lexicographically two points
///                    - `Equal_3` to check whether two points are identical.
///                    For each functor `Foo`, a function `Foo foo_object()` must be provided.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return the number of pairs of halfedges that were stitched.
///
/// \sa `stitch_boundary_cycle()`
/// \sa `stitch_boundary_cycles()`
///
template <typename BorderHalfedgeRange, typename PolygonMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t stitch_borders(const BorderHalfedgeRange& boundary_cycle_representatives,
                           PolygonMesh& pmesh,
                           const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                           , std::enable_if_t<
                               boost::has_range_iterator<BorderHalfedgeRange>::value
                           >* = 0
#endif
                           )
{
  // Need to keep track of the cycles since we are working on a subset of all the boundary cycles
  internal::Boundary_cycle_rep_maintainer<PolygonMesh> cycle_reps_maintainer(pmesh);
  return stitch_borders(boundary_cycle_representatives, pmesh, cycle_reps_maintainer, np);
}


} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_STITCH_BORDERS_H
