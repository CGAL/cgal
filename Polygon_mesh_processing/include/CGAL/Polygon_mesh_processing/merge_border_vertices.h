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

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <boost/unordered_set.hpp>
#include <boost/bind.hpp>

namespace CGAL{

namespace Polygon_mesh_processing{

namespace internal {

#if 0
// warning: vertices will be altered (sorted)
template <class Vpm, class vertex_descriptor>
void detect_identical_vertices(std::vector<vertex_descriptor>& vertices,
                               std::vector< std::vector<vertex_descriptor> >& identical_vertices,
                               Vpm vpm)
{
  typedef typename boost::property_traits<Vpm>::value_type Point_3;

  // sort vertices using their point to ease the detection
  // of vertices with identical points
  CGAL::Property_map_to_unary_function<Vpm> Get_point(vpm);
  std::sort( vertices.begin(), vertices.end(),
             boost::bind(std::less<Point_3>(), boost::bind(Get_point,_1),
                                               boost::bind(Get_point, _2)) );
  std::size_t nbv=vertices.size();
  std::size_t i=1;

  while(i!=nbv)
  {
    if (get(vpm, vertices[i]) == get(vpm, vertices[i-1]))
    {
      identical_vertices.push_back( std::vector<vertex_descriptor>() );
      identical_vertices.back().push_back(vertices[i-1]);
      identical_vertices.back().push_back(vertices[i]);
      while(++i!=nbv)
      {
        if (get(vpm, vertices[i]) == get(vpm, vertices[i-1]))
          identical_vertices.back().push_back(vertices[i]);
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
}
#endif

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

  bool operator()(halfedge_descriptor h1,
                  halfedge_descriptor h2) const
  {
    return get(vpm, target(h1, pm)) < get(vpm, target(h2, pm));
  }

  const PM& pm;
  const VertexPointMap& vpm;
};


// warning: cycle_hedges will be altered (sorted)
template <class PolygonMesh, class Vpm, class halfedge_descriptor>
void detect_identical_vertices(std::vector<halfedge_descriptor>& cycle_hedges,
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

  while(i!=nbv)
  {
    if ( get(vpm, target(cycle_hedges[i], pm)) ==
         get(vpm, target(cycle_hedges[i-1], pm)) )
    {
      hedges_with_identical_point_target.push_back( std::vector<halfedge_descriptor>() );
      hedges_with_identical_point_target.back().push_back(cycle_hedges[i-1]);
      hedges_with_identical_point_target.back().push_back(cycle_hedges[i]);
      while(++i!=nbv)
      {
        if ( get(vpm, target(cycle_hedges[i], pm)) ==
             get(vpm, target(cycle_hedges[i-1], pm)) )
          hedges_with_identical_point_target.back().push_back(cycle_hedges[i]);
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
}

} // end of internal

/// \todo document me
/// It should probably go into BGL package
/// It should make sense to also return the length of each cycle
template <typename PolygonMesh, typename OutputIterator>
OutputIterator
extract_boundary_cycles(PolygonMesh& pm,
                        OutputIterator out)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  boost::unordered_set<halfedge_descriptor> hedge_handled;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pm))
  {
    if(is_border(h, pm) && hedge_handled.insert(h).second)
    {
      *out++=h;
      BOOST_FOREACH(halfedge_descriptor h2, halfedges_around_face(h, pm))
        hedge_handled.insert(h2);
    }
  }
  return out;
}

/// \ingroup PMP_repairing_grp
/// \todo document me
/// we merge the all the target of the halfedges in `hedges`
/// hedges must be sorted along the cycle
template <typename PolygonMesh, class HalfedgeRange>
void merge_boundary_vertices_in_cycle(const HalfedgeRange& sorted_hedges,
                                            PolygonMesh& pm)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  halfedge_descriptor in_h_kept = *boost::begin(sorted_hedges);
  halfedge_descriptor out_h_kept = next(in_h_kept, pm);
  vertex_descriptor v_kept=target(in_h_kept, pm);

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

/// \ingroup PMP_repairing_grp
/// \todo document me
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
  std::vector<halfedge_descriptor> cycle_hedges;
  halfedge_descriptor start=h;
  do{
    cycle_hedges.push_back(h);
    h=next(h, pm);
  }while(start!=h);

  std::vector< std::vector<halfedge_descriptor> > hedges_with_identical_point_target;
  internal::detect_identical_vertices(cycle_hedges, hedges_with_identical_point_target, pm, vpm);

  BOOST_FOREACH(const std::vector<halfedge_descriptor>& hedges,
                hedges_with_identical_point_target)
  {
    start=hedges.front();
    // collect all halfedges in the cycle
    std::vector<halfedge_descriptor> sorted_hedges;
    h=start;
    do{
      sorted_hedges.push_back(h);
      do
      {
        h=next(h, pm);
      }
      while( get(vpm, target(h, pm)) != get(vpm, target(start, pm)) );
    }
    while(h!=start);

    if (sorted_hedges.size() != hedges.size())
    {
      std::cerr << "WARNING: cycle broken at " << get(vpm, target(start, pm)) << ". Skipped\n";
      std::cout << sorted_hedges.size() << " vs " << hedges.size() << "\n";
      CGAL_assertion(sorted_hedges.size() == hedges.size());
      continue;
    }
    merge_boundary_vertices_in_cycle(sorted_hedges, pm);
  }
}

/// \ingroup PMP_repairing_grp
/// \todo document me
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

#if 0
/// \ingroup PMP_repairing_grp
/// \todo document me
template <class PolygonMesh, class NamedParameter>
void merge_duplicated_boundary_vertices(      PolygonMesh& pm,
                                        const NamedParameter& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameter>::const_type Vpm;

  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(vertex_point, pm));

  std::vector<vertex_descriptor> border_vertices;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pm))
  {
    if(is_border(h, pm))
      border_vertices.push_back(target(h, pm));
  }

  std::vector< std::vector<vertex_descriptor> > identical_vertices;
  internal::detect_identical_vertices(border_vertices, identical_vertices, vpm);

  BOOST_FOREACH(const std::vector<vertex_descriptor>& vrtcs, identical_vertices)
    merge_boundary_vertices(vrtcs, pm);
}
#endif

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

#if 0
template <class PolygonMesh>
void merge_duplicated_boundary_vertices(PolygonMesh& pm)
{
  merge_duplicated_boundary_vertices(pm, parameters::all_default());
}
#endif

} } // end of CGAL::Polygon_mesh_processing

#endif //CGAL_POLYGON_MESH_PROCESSING_MERGE_BORDER_VERTICES_H
