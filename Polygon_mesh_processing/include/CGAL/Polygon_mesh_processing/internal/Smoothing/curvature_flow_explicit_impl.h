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
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_EXPLICIT_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_EXPLICIT_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_helpers.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/property_map.h>
#include <CGAL/iterator.h>

#include <boost/graph/graph_traits.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename TriangleMesh,
         typename VertexPointMap,
         typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<TriangleMesh, VertexPointMap> >
class Cotangent_weight
  : public CotangentValue
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor     halfedge_descriptor;

public:
  Cotangent_weight(TriangleMesh& pmesh_, VertexPointMap vpmap_) : CotangentValue(pmesh_, vpmap_) {}

  TriangleMesh& pmesh() { return CotangentValue::pmesh(); }

  double operator()(halfedge_descriptor he)
  {
    vertex_descriptor vs = source(he, pmesh());
    vertex_descriptor vt = target(he, pmesh());
    vertex_descriptor v1 = target(next(he, pmesh()), pmesh());
    vertex_descriptor v2 = source(prev(opposite(he, pmesh()), pmesh()), pmesh());

    return (CotangentValue::operator()(vt, v1, vs) + CotangentValue::operator()(vs, v2, vt));
  }
};

template<typename TriangleMesh,
         typename VertexPointMap, typename VertexConstraintMap,
         typename GeomTraits>
class Curvature_flow
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor     halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor         edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor         face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type         Point;
  typedef typename boost::property_traits<VertexPointMap>::reference          Point_ref;

  typedef typename GeomTraits::Vector_3                                       Vector;
  typedef typename GeomTraits::Triangle_3                                     Triangle;
  typedef std::vector<Triangle>                                               Triangle_list;

  typedef Cotangent_weight<TriangleMesh, VertexPointMap>                      Weight_calculator;

public:
  Curvature_flow(TriangleMesh& pmesh,
                 VertexPointMap vpmap,
                 VertexConstraintMap vcmap)
    :
      mesh_(pmesh),
      vpmap_(vpmap), vcmap_(vcmap),
      weight_calculator_(pmesh, vpmap)
  { }

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    CGAL_precondition(CGAL::is_triangle_mesh(mesh_));
    CGAL_precondition_code(std::set<typename boost::graph_traits<TriangleMesh>::face_descriptor> degen_faces;)
        CGAL_precondition_code(CGAL::Polygon_mesh_processing::degenerate_faces(
                                 mesh_, std::inserter(degen_faces, degen_faces.begin()),
                                 CGAL::parameters::vertex_point_map(vpmap_).geom_traits(traits_));)
        CGAL_precondition(degen_faces.empty());

    set_vertex_range(face_range);
  }

  // @todo should be possible to pass a time step because this is not a very stable process
  void curvature_smoothing()
  {
    // In general, we will smooth whole meshes, so this is better than a map
    typedef CGAL::dynamic_vertex_property_t<Point>                      Vertex_property_tag;
    typedef typename boost::property_map<TriangleMesh,
                                         Vertex_property_tag>::type     Position_map;

    Position_map new_positions = get(Vertex_property_tag(), mesh_);

    for(vertex_descriptor v : vrange_)
    {
      if(is_border(v, mesh_) || is_constrained(v))
        continue;

      // calculate movement
      Vector move = CGAL::NULL_VECTOR;
      double sum_cot_weights = 0;
      for(halfedge_descriptor hi : halfedges_around_source(v, mesh_))
      {
        // weight
        const double weight = weight_calculator_(hi);
        sum_cot_weights += weight;

        // displacement vector
        const Point_ref xi = get(vpmap_, source(hi, mesh_));
        const Point_ref xj = get(vpmap_, target(hi, mesh_));

        Vector vec = traits_.construct_vector_3_object()(xi, xj);
        vec = traits_.construct_scaled_vector_3_object()(vec, weight);

        move = traits_.construct_sum_of_vectors_3_object()(move, vec);
      }

      // divide with total weight
      if(sum_cot_weights != 0.)
        move = traits_.construct_scaled_vector_3_object()(move, 1./sum_cot_weights);

      Point new_pos = traits_.construct_translated_point_3_object()(get(vpmap_, v), move);
      put(new_positions, v, new_pos);
    }

    // update locations
    for(vertex_descriptor v : vrange_)
      put(vpmap_, v, get(new_positions, v));
  }

private:
  bool is_constrained(const vertex_descriptor v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void set_vertex_range(const FaceRange& face_range)
  {
    vrange_.reserve(3 * face_range.size());
    for(face_descriptor f : face_range)
    {
      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
        vrange_.push_back(v);
    }

    // get rid of duplicate vertices
    std::sort(vrange_.begin(), vrange_.end());
    vrange_.erase(std::unique(vrange_.begin(), vrange_.end()), vrange_.end());
  }

private:
  // data members
  // ------------
  TriangleMesh& mesh_;

  VertexPointMap vpmap_;
  VertexConstraintMap vcmap_;

  std::vector<vertex_descriptor> vrange_;

  Weight_calculator weight_calculator_;
  GeomTraits traits_;
};

} // internal
} // Polygon_mesh_processing
} // CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_EXPLICIT_IMPL_H
