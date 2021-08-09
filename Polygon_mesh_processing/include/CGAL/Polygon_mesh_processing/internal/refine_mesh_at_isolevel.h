// Copyright (c) 2021 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_CLIP_H
#define CGAL_POLYGON_MESH_PROCESSING_CLIP_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace experimental {


template <class PolygonMesh, class ValueMap, class NamedParameters>
void refine_mesh_at_isolevel(PolygonMesh& pm,
                             ValueMap value_map,
                             typename boost::property_traits<ValueMap>::value_type isovalue,
                             const NamedParameters& np)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  typedef Static_boolean_property_map<edge_descriptor, false>  Default_ECM;
  typedef typename internal_np::Lookup_named_param_def<internal_np::edge_is_constrained_t,
                                                       NamedParameters,
                                                       Default_ECM>::type       ECM;
  typedef typename GetVertexPointMap < PolygonMesh, NamedParameters>::type VPM;

  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, pm));

  ECM ecm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained), Default_ECM());

  std::unordered_map<face_descriptor, std::vector<halfedge_descriptor> > faces_to_split;
  std::vector<edge_descriptor> to_split;
  for (edge_descriptor e : edges(pm))
  {
    vertex_descriptor src = source(e, pm), tgt = target(e, pm);
    if (get(value_map, src)==isovalue)
    {
      for (halfedge_descriptor h : halfedges_around_source(halfedge(e, pm), pm))
      {
        face_descriptor f = face(h, pm);
        if (f!=boost::graph_traits<PolygonMesh>::null_face())
          faces_to_split[f].push_back(opposite(h, pm));
      }
      continue;
    }
    if (get(value_map, tgt)==isovalue)
    {
      for (halfedge_descriptor h : halfedges_around_target(halfedge(e, pm), pm))
      {
        face_descriptor f = face(h, pm);
        if (f!=boost::graph_traits<PolygonMesh>::null_face())
          faces_to_split[f].push_back(h);
      }
      continue;
    }
    if ( (get(value_map, tgt) < isovalue) != (get(value_map, src) < isovalue) )
    {
      to_split.push_back(e);
    }
  }

  for (edge_descriptor e : to_split)
  {
    vertex_descriptor src = source(e, pm), tgt = target(e, pm);
    double ds = get(value_map, src);
    double dt = get(value_map, tgt);
    double alpha = (isovalue - dt) / (ds - dt);
    halfedge_descriptor hnew = CGAL::Euler::split_edge(halfedge(e, pm), pm);
    put(vpm, target(hnew, pm), barycenter(get(vpm,src), alpha, get(vpm, tgt), 1-alpha));
    put(value_map,  target(hnew, pm) , isovalue);
    face_descriptor f = face(hnew, pm);
    if (f!=boost::graph_traits<PolygonMesh>::null_face())
      faces_to_split[f].push_back(hnew);
    hnew=pm.prev(opposite(hnew, pm));
    f = face(hnew, pm);
    if (f!=boost::graph_traits<PolygonMesh>::null_face())
      faces_to_split[f].push_back(hnew);
  }

  for (const auto& p : faces_to_split)
  {
    if(p.second.size()!=2) continue;

    std::pair<edge_descriptor ,bool> res = edge(target(p.second[0],pm),
                                                target(p.second[1],pm), pm);
    if (res.second)
    {
      // no split as the edge already exists (the two vertices are on the isolevel)
      put(ecm, res.first, true);
      continue;
    }

    halfedge_descriptor hnew = CGAL::Euler::split_face(p.second[0], p.second[1], pm);
    put(ecm, edge(hnew, pm), true);
  }
}

template <class PolygonMesh, class ValueMap>
void refine_mesh_at_isolevel(PolygonMesh& pm,
                             ValueMap value_map,
                             typename boost::property_traits<ValueMap>::value_type iso_level)
{
  refine_mesh_at_isolevel(pm, value_map, iso_level, parameters::all_default());
}

} } } // end of CGAL::Polygon_mesh_processing::experimental


#endif
