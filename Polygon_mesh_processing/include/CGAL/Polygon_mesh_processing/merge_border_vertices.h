// Copyright (c) 2018 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H
#define CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <boost/unordered_set.hpp>
#include <boost/bind.hpp>

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

  std::set< std::pair<std::size_t, std::size_t> > intervals;

  while(i!=nbv)
  {
    if ( get(vpm, target(cycle_hedges[i].first, pm)) ==
         get(vpm, target(cycle_hedges[i-1].first, pm)) )
    {
      hedges_with_identical_point_target.push_back( std::vector<halfedge_descriptor>() );
      hedges_with_identical_point_target.back().push_back(cycle_hedges[i-1].first);
      hedges_with_identical_point_target.back().push_back(cycle_hedges[i].first);
      intervals.insert( std::make_pair(cycle_hedges[i-1].second, cycle_hedges[i].second) );
      std::size_t previous = cycle_hedges[i].second;
      while(++i!=nbv)
      {
        if ( get(vpm, target(cycle_hedges[i].first, pm)) ==
             get(vpm, target(cycle_hedges[i-1].first, pm)) )
        {
          hedges_with_identical_point_target.back().push_back(cycle_hedges[i].first);
          intervals.insert( std::make_pair(previous, cycle_hedges[i].second) );
          previous = cycle_hedges[i].second;
        }
        else
        {
          ++i;
          break;
        }
      }
    }
    else
      ++i;
  }

  // check that intervals are disjoint or strictly nested
  // if there is only one issue we drop the whole cycle.
  /// \todo shall we try to be more conservative?
  if (hedges_with_identical_point_target.empty()) return;
  std::set< std::pair<std::size_t, std::size_t> >::iterator it1 = intervals.begin(),
                                                            end2 = intervals.end(),
                                                            end1 = cpp11::prev(end2),
                                                            it2;
  for (; it1!=end1; ++it1)
    for(it2=cpp11::next(it1); it2!= end2; ++it2 )
    {
      CGAL_assertion(it1->first<it2->first);
      CGAL_assertion(it1->first < it1->second && it2->first < it2->second);
      if (it1->second > it2->first && it2->second > it1->second)
      {
        std::cerr << "Merging is skipt to avoid bad cycle connections\n";
        hedges_with_identical_point_target.clear();
        return;
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

  BOOST_FOREACH(halfedge_descriptor in_h_rm, sorted_hedges)
  {
    vertex_descriptor vd = target(in_h_rm, pm);
    if (vd==v_kept) continue; // skip identical vertices (in particular this skips the first halfedge)
    if (edge(vd, v_kept, pm).second) continue; // skip null edges
    bool shall_continue=false;
    BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(v_kept, pm))
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

  BOOST_FOREACH(vertex_descriptor vd, vertices_to_rm)
    remove_vertex(vd, pm);
}

} // end of internal

/// \ingroup PMP_repairing_grp
/// merges identical vertices around a cycle of boundary edges.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`.
/// @tparam NamedParameter a sequence of \ref pmp_namedparameters "Named Parameters".
///
/// @param h a halfedge that belongs to a boundary cycle.
/// @param pm the polygon mesh which contains the boundary cycle.
/// @param np optional parameter of \ref pmp_namedparameters "Named Parameters" listed below.
///
/// \cgalNamedParamsBegin
/// \cgalParamBegin{vertex_point_map}
///   the property map with the points associated to the vertices of `pm`.
///    If this parameter is omitted, an internal property map for
///  `CGAL::vertex_point_t` should be available in `PolygonMesh`
/// \cgalParamEnd
/// \cgalNamedParamsEnd
template <class PolygonMesh, class NamedParameter>
void merge_duplicated_vertices_in_boundary_cycle(
        typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
        PolygonMesh& pm,
        const NamedParameter& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameter>::const_type Vpm;

  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
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

  BOOST_FOREACH(const std::vector<halfedge_descriptor>& hedges,
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
/// @tparam NamedParameter a sequence of \ref pmp_namedparameters "Named Parameters".
///
/// @param pm the polygon mesh which contains the cycles.
/// @param np optional parameter of \ref pmp_namedparameters "Named Parameters" listed below.
///
/// \cgalNamedParamsBegin
/// \cgalParamBegin{vertex_point_map}
///   the property map with the points associated to the vertices of `pm`.
///    If this parameter is omitted, an internal property map for
///  `CGAL::vertex_point_t` should be available in `PolygonMesh`
/// \cgalParamEnd
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

  BOOST_FOREACH(halfedge_descriptor h, cycles)
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
