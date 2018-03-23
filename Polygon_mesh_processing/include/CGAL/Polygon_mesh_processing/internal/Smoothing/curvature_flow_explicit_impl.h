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
//
//
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_CURVATURE_FLOW_EXPLICIT_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_CURVATURE_FLOW_EXPLICIT_IMPL_H

#include <utility>
#include <math.h>

#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/property_map.h>
#include <CGAL/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh, typename VertexPointMap,
      typename CotangentValue = CGAL::internal::Cotangent_value_Meyer<PolygonMesh, VertexPointMap> >
class Cotangent_weight : CotangentValue
{
public:
  Cotangent_weight(PolygonMesh& pmesh_, VertexPointMap vpmap_)
      : CotangentValue(pmesh_, vpmap_)
  {}

  PolygonMesh& pmesh()
  {
      return CotangentValue::pmesh();
  }

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;
  typedef std::pair<halfedge_descriptor, halfedge_descriptor>              he_pair;

  double operator()(halfedge_descriptor he, he_pair incd_edges)
  {
      vertex_descriptor vs = source(he, pmesh());
      vertex_descriptor vt = target(he, pmesh());
      vertex_descriptor v1 = target(incd_edges.first, pmesh());
      vertex_descriptor v2 = source(incd_edges.second, pmesh());

      return ( CotangentValue::operator()(vs, v1, vt) + CotangentValue::operator()(vs, v2, vt) );
  }
};

template<typename PolygonMesh, typename VertexPointMap, typename VertexConstraintMap, typename GeomTraits>
class Curvature_flow
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;

  typedef typename GeomTraits::Point_3  Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Triangle_3 Triangle;
  typedef std::vector<Triangle> Triangle_list;

  typedef std::pair<halfedge_descriptor, halfedge_descriptor> he_pair;
  typedef std::map<halfedge_descriptor, he_pair> Edges_around_map;

  typedef Cotangent_weight<PolygonMesh, VertexPointMap> Weight_calculator;

public:
  Curvature_flow(PolygonMesh& pmesh, VertexPointMap& vpmap, VertexConstraintMap& vcmap) :
    mesh_(pmesh), vpmap_(vpmap), vcmap_(vcmap),
    weight_calculator_(pmesh, vpmap) {}

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    check_vertex_range(face_range);
    BOOST_FOREACH(face_descriptor f, face_range)
      input_triangles_.push_back(triangle(f));
  }

  std::size_t remove_degenerate_faces()
  {
    std::size_t nb_removed_faces = 0;
    nb_removed_faces = CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh_);
    return nb_removed_faces;
  }

  void curvature_smoothing()
  {
    std::map<vertex_descriptor, Point> barycenters;
    BOOST_FOREACH(vertex_descriptor v, vrange)
    {
      if(!is_border(v, mesh_) && !is_constrained(v))
      {
        // find incident halfedges
        Edges_around_map he_map;
        typename Edges_around_map::iterator it;
        BOOST_FOREACH(halfedge_descriptor hi, halfedges_around_source(v, mesh_))
          he_map[hi] = he_pair( next(hi, mesh_), prev(opposite(hi, mesh_), mesh_) );

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
        if(sum_cot_weights != 0)
            curvature_normal /= sum_cot_weights;

        Point weighted_barycenter = get(vpmap_, v) - curvature_normal;
        barycenters[v] = weighted_barycenter;

      } // not on border
    } // all vertices

    // update location
    typedef typename std::map<vertex_descriptor, Point>::value_type VP;
    BOOST_FOREACH(const VP& vp, barycenters)
      put(vpmap_, vp.first, vp.second);
  }

private:
  // helper functions
  // ----------------
  Triangle triangle(face_descriptor f) const
  {
    halfedge_descriptor h = halfedge(f, mesh_);
    vertex_descriptor v1 = target(h, mesh_);
    vertex_descriptor v2 = target(next(h, mesh_), mesh_);
    vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
    return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
  }

  double sqlength(const vertex_descriptor& v1, const vertex_descriptor& v2) const
  {
    return to_double(CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2)));
  }

  double sqlength(const halfedge_descriptor& h) const
  {
      vertex_descriptor v1 = target(h, mesh_);
      vertex_descriptor v2 = source(h, mesh_);
      return sqlength(v1, v2);
  }

  double sqlength(const edge_descriptor& e) const
  {
      return sqlength(halfedge(e, mesh_));
  }

  // degeneracy removal
  // ------------------
  void check_degeneracy(halfedge_descriptor h1)
  {
    halfedge_descriptor h2 = next(h1, mesh_);
    halfedge_descriptor h3 = next(h2, mesh_);

    double a1 = get_angle(h1, h2);
    double a2 = get_angle(h2, h3);
    double a3 = get_angle(h3, h1);

    double angle_min_threshold = 0.05; // rad
    double angle_max_threshold = CGAL_PI - 0.05;

    if(a1 < angle_min_threshold || a2 < angle_min_threshold || a3 < angle_min_threshold)
    {
      Euler::remove_face(h1, mesh_);
    }

    if(a1 > angle_max_threshold || a2 > angle_max_threshold || a3 > angle_max_threshold)
    {
      Euler::remove_face(h1, mesh_);
    }
  }

  double get_angle(halfedge_descriptor ha, halfedge_descriptor hb)
  {
    Vector a(get(vpmap_, source(ha, mesh_)), get(vpmap_, target(ha, mesh_)));
    Vector b(get(vpmap_, source(hb, mesh_)), get(vpmap_, target(hb, mesh_)));
    return get_angle(a, b);
  }

  double get_angle(const Vector& e1, const Vector& e2)
  {
    //double rad_to_deg = 180. / CGAL_PI;
    double cos_angle = (e1 * e2)
        / std::sqrt(e1.squared_length() * e2.squared_length());
    return std::acos(cos_angle); //* rad_to_deg;
  }

  bool is_constrained(const vertex_descriptor& v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void check_vertex_range(const FaceRange& face_range)
  {
    BOOST_FOREACH(face_descriptor f, face_range)
    {
      BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, mesh_), mesh_))
        vrange.insert(v);
    }
  }

private:
  // data members
  // ------------
  PolygonMesh& mesh_;
  VertexPointMap& vpmap_;
  VertexConstraintMap vcmap_;
  Triangle_list input_triangles_;
  GeomTraits traits_;
  std::set<vertex_descriptor> vrange;
  Weight_calculator weight_calculator_;
};


} // internal
} // Polygon_mesh_processing
} // CGAL


#endif // CGAL_POLYGON_MESH_PROCESSING_CURVATURE_FLOW_EXPLICIT_IMPL_H
