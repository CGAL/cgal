// Copyright (c) 2021-2023 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_REFINE_MESH_AT_ISOLEVEL_H
#define CGAL_POLYGON_MESH_PROCESSING_REFINE_MESH_AT_ISOLEVEL_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <unordered_set>
#include <unordered_map>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
 * \ingroup PkgPolygonMeshProcessingRef
 *
 * refines `pm` by adding new vertices on edges having their incident vertices associated with
 * values respectively larger and smaller than `isovalue` in `value_map`.
 * The placement of new vertices on edges will be done by linear interpolation
 * using the aforementioned values.
 * New vertices will be associated `isovalue` in `value_map` when created.
 * Additionally, new edges will be added by connecting new vertices created sharing
 * a common incident face. Note that in case more than two new vertices are added
 * on a face boundary, no edges will be created in that face.
 *
 * @tparam PolygonMesh a model of the concepts `EdgeListGraph` and `FaceListGraph`
 * @tparam ValueMap a model of the concept `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                  as key type and with its value type being the type of the coordinates of points associated with vertices
 *                  in the vertex map provided to the `vertex_point_map()` named parameter.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters" for `pm`
 *
 * @param pm the polygon mesh to be refined.
 * @param value_map the property map containing a value at each vertex for a given function defined over the mesh.
 * @param isovalue the value used to refine
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{edge_is_constrained_map}
 *     \cgalParamDescription{an output property map associating `true` to all edges connecting vertices on the isolevel,
 *                           and `false` for all other edges.}
 *     \cgalParamType{a class model of `WritablePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
 *                    as key type and `bool` as value type}
 *     \cgalParamDefault{No marks on edges will be put}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `PolygonMesh`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class PolygonMesh, class ValueMap, class NamedParameters = parameters::Default_named_parameters>
void refine_mesh_at_isolevel(PolygonMesh& pm,
                             ValueMap value_map,
                             typename boost::property_traits<ValueMap>::value_type isovalue,
                             const NamedParameters& np = parameters::default_values())
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::property_traits<ValueMap>::value_type FT;

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
  std::unordered_set<vertex_descriptor> vertices_on_isoline;
  for (edge_descriptor e : edges(pm))
  {
    vertex_descriptor src = source(e, pm), tgt = target(e, pm);

    if (get(value_map, src)==isovalue)
    {
      vertices_on_isoline.insert(src);
      if (get(value_map, tgt)==isovalue)
      {
        put(ecm, e, true); // special case for faces entirely on an isovalue
        continue;
      }
      continue;
    }
    if (get(value_map, tgt)==isovalue)
    {
      vertices_on_isoline.insert(tgt);
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
    FT ds = get(value_map, src);
    FT dt = get(value_map, tgt);
    FT alpha = (isovalue - dt) / (ds - dt);
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

  for (vertex_descriptor vh : vertices_on_isoline)
  {
    for (halfedge_descriptor h : halfedges_around_target(vh, pm))
    {
      face_descriptor f = face(h, pm);
      if (f!=boost::graph_traits<PolygonMesh>::null_face())
        faces_to_split[f].push_back(h);
    }
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

} } // end of CGAL::Polygon_mesh_processing


#endif
