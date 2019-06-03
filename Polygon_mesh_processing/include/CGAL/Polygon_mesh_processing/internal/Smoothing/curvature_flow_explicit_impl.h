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

  typedef std::pair<halfedge_descriptor, halfedge_descriptor>                 he_pair;

public:
  Cotangent_weight(TriangleMesh& pmesh_, VertexPointMap vpmap_) : CotangentValue(pmesh_, vpmap_) {}

  TriangleMesh& pmesh() { return CotangentValue::pmesh(); }

  double operator()(halfedge_descriptor he, he_pair incd_edges)
  {
    vertex_descriptor vs = source(he, pmesh());
    vertex_descriptor vt = target(he, pmesh());
    vertex_descriptor v1 = target(incd_edges.first, pmesh());
    vertex_descriptor v2 = source(incd_edges.second, pmesh());

    return ( CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt) );
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

  typedef typename GeomTraits::Point_3                                        Point;
  typedef typename GeomTraits::Vector_3                                       Vector;
  typedef typename GeomTraits::Triangle_3                                     Triangle;
  typedef std::vector<Triangle>                                               Triangle_list;

  typedef std::pair<halfedge_descriptor, halfedge_descriptor>                 he_pair;
  typedef std::map<halfedge_descriptor, he_pair>                              Edges_around_map;

  typedef Cotangent_weight<TriangleMesh, VertexPointMap>                      Weight_calculator;

public:
  Curvature_flow(TriangleMesh& pmesh, VertexPointMap& vpmap, VertexConstraintMap& vcmap)
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

    input_triangles_.clear();
    input_triangles_.reserve(face_range.size());

    for(face_descriptor f : face_range)
    {
      Triangle t;
      construct_triangle(f, mesh_, t);
      input_triangles_.push_back(t);
    }
  }

  void curvature_smoothing()
  {
    boost::unordered_map<vertex_descriptor, Point> barycenters;
    for(vertex_descriptor v : vrange_)
    {
      if(!is_border(v, mesh_) && !is_constrained(v))
      {
        // find incident halfedges
        Edges_around_map he_map;
        typename Edges_around_map::iterator it;
        for(halfedge_descriptor hi : halfedges_around_source(v, mesh_))
          he_map[hi] = std::make_pair(next(hi, mesh_), prev(opposite(hi, mesh_), mesh_));

        // calculate movement
        Vector curvature_normal = CGAL::NULL_VECTOR;
        double sum_cot_weights = 0;
        for(it = he_map.begin(); it!= he_map.end(); ++it)
        {
          halfedge_descriptor hi = it->first;
          he_pair incd_edges = it->second;

          // weight
          double weight = weight_calculator_(hi, incd_edges);;
          sum_cot_weights += weight;

          // displacement vector
          Point xi = get(vpmap_, source(hi, mesh_));
          Point xj = get(vpmap_, target(hi, mesh_));
          Vector vec(xj, xi); // towards the vertex that is being moved

          // add weight
          vec *= weight;

          // sum vecs
          curvature_normal += vec;
        }

        // divide with total weight
        if(sum_cot_weights != 0.)
          curvature_normal /= sum_cot_weights;

        Point weighted_barycenter = get(vpmap_, v) - curvature_normal;
        barycenters[v] = weighted_barycenter;

      } // not on border
    } // all vertices

    // update location
    typedef typename boost::unordered_map<vertex_descriptor, Point>::value_type VP;
    for(const VP& vp : barycenters)
      put(vpmap_, vp.first, vp.second);
  }

private:
  bool is_constrained(const vertex_descriptor& v)
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

  VertexPointMap& vpmap_;
  VertexConstraintMap vcmap_;

  Triangle_list input_triangles_;
  std::vector<vertex_descriptor> vrange_;

  Weight_calculator weight_calculator_;
  GeomTraits traits_;
};

} // internal
} // Polygon_mesh_processing
} // CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_EXPLICIT_IMPL_H
