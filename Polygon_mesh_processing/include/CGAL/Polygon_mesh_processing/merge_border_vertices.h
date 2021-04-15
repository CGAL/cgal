// Copyright (c) 2018 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H
#define CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>
#include <utility>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template <typename PM, typename VertexPointMap>
struct Less_on_point_of_target
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor
    halfedge_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::reference Point;

  Less_on_point_of_target(const PM& pm,
                          const VertexPointMap& vpm)
    : pm(pm),
      vpm(vpm)
  {}

  bool operator()(const std::pair<halfedge_descriptor, std::size_t>& h1,
                  const std::pair<halfedge_descriptor, std::size_t>& h2) const
  {
    if ( get(vpm, target(h1.first, pm)) < get(vpm, target(h2.first, pm)) )
      return true;
    if ( get(vpm, target(h1.first, pm)) > get(vpm, target(h2.first, pm)) )
      return false;
    return h1.second < h2.second;
  }

  const PM& pm;
  const VertexPointMap& vpm;
};

// Given a container of vectors of halfedges whose target are geometrically indentical,
// check that the intervals described by these pairs are either disjoint or nested.
// This is done to ensure valid combinatorics when we merge the vertices.
// If incompatible (overlapping) intervals are found, the pair representating the longest
// interval (arbitrary choice) is removed from the candidate list.
template <typename VPM, typename PolygonMesh>
void sanitize_candidates(const std::vector<std::pair<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor, std::size_t> >& cycle_hedges,
                         std::vector<std::vector<std::size_t> >& candidate_hedges_with_id,
                         const VPM vpm,
                         const PolygonMesh& pm)
{
  if(candidate_hedges_with_id.empty())
    return;

  std::size_t nm_vertices_n = candidate_hedges_with_id.size();
  for(std::size_t fr_id=0, fr_end=nm_vertices_n-1; fr_id<fr_end; ++fr_id)
  {
    std::vector<std::size_t>& first_candidates = candidate_hedges_with_id[fr_id];
    CGAL_assertion(first_candidates.size() >= 2);

    for(std::size_t i=0, ie=first_candidates.size()-1; i<ie; ++i)
    {
      const std::size_t first_left = cycle_hedges[first_candidates[i]].second;
      const std::size_t first_right = cycle_hedges[first_candidates[i+1]].second;
      CGAL_assertion(first_left < first_right);

      for(std::size_t sr_id=i+1, sr_end=nm_vertices_n; sr_id<sr_end; ++sr_id)
      {
        std::vector<std::size_t>& second_candidates = candidate_hedges_with_id[sr_id];
        CGAL_assertion(second_candidates.size() >= 2);

        for(std::size_t j=0, je=second_candidates.size()-1; j<je; ++j)
        {
          const std::size_t second_left = cycle_hedges[second_candidates[j]].second;
          const std::size_t second_right = cycle_hedges[second_candidates[j+1]].second;
          CGAL_assertion(second_left < second_right);

          // The pair of intervals should be either disjoint or nested
          // so reject:
          // sl -- fl -- sr -- fr and fl -- sl -- fr -- sr
          if((second_left < first_left && first_left < second_right && second_right < first_right) ||
             (first_left < second_left && second_left < first_right && first_right < second_right))
          {
            // Remove the candidate with largest range
            const std::size_t first_candidates_range =
                cycle_hedges[first_candidates.back()].second - cycle_hedges[first_candidates.front()].second;
            const std::size_t second_candidates_range =
                cycle_hedges[second_candidates.back()].second - cycle_hedges[second_candidates.front()].second;

            CGAL_assertion(first_candidates_range <= cycle_hedges.size());
            CGAL_assertion(second_candidates_range <= cycle_hedges.size());

#ifdef CGAL_PMP_MERGE_BORDER_VERTICES_DEBUG
            std::cout << "Incompatible ranges:\n";
            std::cout << "first range: " << first_left << " to " << first_right << std::endl;
            std::cout << "second range: " << second_left << " to " << second_right << std::endl;
            std::cout << "Full ranges:" << std::endl;
            std::cout << cycle_hedges[first_candidates.front()].second << " to " << cycle_hedges[first_candidates.back()].second;
            std::cout << " (" << first_candidates.size() << " halfedges)";
            std::cout << " at " << get(vpm, target(cycle_hedges[first_candidates.front()].first, pm)) << std::endl;
            std::cout << cycle_hedges[second_candidates.front()].second << " to " << cycle_hedges[second_candidates.back()].second;
            std::cout << " (" << second_candidates.size() << " halfedges)";
            std::cout << " at " << get(vpm, target(cycle_hedges[second_candidates.front()].first, pm)) << std::endl;
#endif

            std::vector<std::vector<std::size_t> >::iterator to_remove_iter = candidate_hedges_with_id.begin();
            if(first_candidates_range > second_candidates_range)
              std::advance(to_remove_iter, fr_id);
            else
              std::advance(to_remove_iter, sr_id);

            candidate_hedges_with_id.erase(to_remove_iter);

            // restart the whole thing
            return sanitize_candidates(cycle_hedges, candidate_hedges_with_id, vpm, pm);
          }
        } // entries of the second range
      } // second range
    } // entries of the first range
  } // first range
}

// warning: cycle_hedges will be altered (sorted)
template <class PolygonMesh, class Vpm, class halfedge_descriptor>
void detect_identical_mergeable_vertices(
        std::vector< std::pair<halfedge_descriptor, std::size_t> >& cycle_hedges,
        std::vector< std::vector<halfedge_descriptor> >& hedges_with_identical_point_target,
        const PolygonMesh& pm,
        Vpm vpm)
{
  // sort vertices using their point to ease the detection
  // of vertices with identical points
  Less_on_point_of_target<PolygonMesh, Vpm> less(pm, vpm);
  std::sort( cycle_hedges.begin(), cycle_hedges.end(), less);

  std::size_t nbv=cycle_hedges.size();
  std::size_t i=1;

  // IDs of cycle_hedges
  std::vector<std::vector<std::size_t> > candidate_hedges_with_id;

  while(i != nbv)
  {
    if(get(vpm, target(cycle_hedges[i].first, pm)) ==
       get(vpm, target(cycle_hedges[i-1].first, pm)) )
    {
      candidate_hedges_with_id.resize(candidate_hedges_with_id.size() + 1);
      candidate_hedges_with_id.back().push_back(i-1);
      candidate_hedges_with_id.back().push_back(i);
      while(++i != nbv)
      {
        if(get(vpm, target(cycle_hedges[i].first, pm)) ==
           get(vpm, target(cycle_hedges[i-1].first, pm)))
        {
          candidate_hedges_with_id.back().push_back(i);
        }
        else
        {
          ++i;
          break;
        }
      }
    }
    else
    {
      ++i;
    }
  }

  // Check that intervals are disjoint or strictly nested
  sanitize_candidates(cycle_hedges, candidate_hedges_with_id, vpm, pm);

  for(const std::vector<std::size_t>& candidates : candidate_hedges_with_id)
  {
    hedges_with_identical_point_target.resize(hedges_with_identical_point_target.size() + 1);
    for(const std::size_t hid : candidates)
    {
      hedges_with_identical_point_target.back().push_back(cycle_hedges[hid].first);
    }
  }
}

// \ingroup PMP_repairing_grp
// merges target vertices of a list of halfedges.
// Halfedges must be sorted in the list.
//
// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`.
// @tparam HalfedgeRange a range of halfedge descriptors of `PolygonMesh`, model of `Range`.
//
// @param sorted_hedges a sorted list of halfedges.
// @param pm the polygon mesh which contains the list of halfedges.
//
template <typename PolygonMesh, class HalfedgeRange>
void merge_vertices_in_range(const HalfedgeRange& sorted_hedges,
                             PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  halfedge_descriptor in_h_kept = *boost::begin(sorted_hedges);
  halfedge_descriptor out_h_kept = next(in_h_kept, pm);
  vertex_descriptor v_kept = target(in_h_kept, pm);

  std::vector<vertex_descriptor> vertices_to_rm;

  for(halfedge_descriptor in_h_rm : sorted_hedges)
  {
    vertex_descriptor vd = target(in_h_rm, pm);
    if (vd==v_kept) continue; // skip identical vertices (in particular this skips the first halfedge)
    if (edge(vd, v_kept, pm).second) continue; // skip null edges
    bool shall_continue=false;
    for(halfedge_descriptor h : halfedges_around_target(v_kept, pm))
    {
      if (edge(vd, source(h, pm), pm).second)
      {
        shall_continue=true;
        break;
      }
    }
    if (shall_continue) continue; // skip vertices already incident to the same vertex
    // update the vertex of the halfedges incident to the vertex to remove
    internal::update_target_vertex(in_h_rm, v_kept, pm);
    // update next/prev pointers around the 2 vertices to be merged
    halfedge_descriptor out_h_rm = next(in_h_rm, pm);
    set_next(in_h_kept, out_h_rm, pm);
    set_next(in_h_rm, out_h_kept, pm);
    vertices_to_rm.push_back(vd);
    out_h_kept=out_h_rm;
  }

  for(vertex_descriptor vd : vertices_to_rm)
    remove_vertex(vd, pm);
}

} // end of internal

/// \ingroup PMP_repairing_grp
/// merges identical vertices around a cycle of boundary edges.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`.
/// @tparam NamedParameter a sequence of \ref bgl_namedparameters "Named Parameters".
///
/// @param h a halfedge that belongs to a boundary cycle.
/// @param pm the polygon mesh which contains the boundary cycle.
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pm`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
template <class PolygonMesh, class NamedParameter>
void merge_duplicated_vertices_in_boundary_cycle(
        typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
        PolygonMesh& pm,
        const NamedParameter& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameter>::const_type Vpm;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(vertex_point, pm));

  // collect all the halfedges of the cycle
  std::vector< std::pair<halfedge_descriptor, std::size_t> > cycle_hedges;
  halfedge_descriptor start=h;
  std::size_t index=0;
  do{
    cycle_hedges.push_back( std::make_pair(h, index) );
    h=next(h, pm);
    ++index;
  }while(start!=h);

  std::vector< std::vector<halfedge_descriptor> > hedges_with_identical_point_target;
  internal::detect_identical_mergeable_vertices(cycle_hedges, hedges_with_identical_point_target, pm, vpm);

  for(const std::vector<halfedge_descriptor>& hedges :
                hedges_with_identical_point_target)
  {
    start=hedges.front();
    // hedges are sorted along the cycle
    internal::merge_vertices_in_range(hedges, pm);
  }
}

/// \ingroup PMP_repairing_grp
/// extracts boundary cycles and merges the duplicated vertices of each cycle.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`.
/// @tparam NamedParameter a sequence of \ref bgl_namedparameters "Named Parameters".
///
/// @param pm the polygon mesh which contains the cycles.
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pm`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \sa `merge_duplicated_vertices_in_boundary_cycle()`
template <class PolygonMesh, class NamedParameter>
void merge_duplicated_vertices_in_boundary_cycles(      PolygonMesh& pm,
                                                  const NamedParameter& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  std::vector<halfedge_descriptor> cycles;
  extract_boundary_cycles(pm, std::back_inserter(cycles));

  for(halfedge_descriptor h : cycles)
    merge_duplicated_vertices_in_boundary_cycle(h, pm, np);
}

template <class PolygonMesh>
void merge_duplicated_vertices_in_boundary_cycles(PolygonMesh& pm)
{
  merge_duplicated_vertices_in_boundary_cycles(pm, parameters::all_default());
}

template <class PolygonMesh>
void merge_duplicated_vertices_in_boundary_cycle(
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
  PolygonMesh& pm)
{
  merge_duplicated_vertices_in_boundary_cycle(h, pm, parameters::all_default());
}

} } // end of CGAL::Polygon_mesh_processing

#endif //CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H
