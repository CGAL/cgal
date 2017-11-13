// Copyright (c) 2016 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_CALLBACK_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_CALLBACK_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/property_map.h>
#include <CGAL/enum.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <set>

namespace CGAL {
namespace Corefinement {

template<class TriangleMesh, class EdgeToFaces>
class Collect_face_bbox_per_edge_bbox {
protected:
  const TriangleMesh& tm_faces;
  const TriangleMesh& tm_edges;
  EdgeToFaces& edge_to_faces;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor> Box;

public:
  Collect_face_bbox_per_edge_bbox(
    const TriangleMesh& tm_faces,
    const TriangleMesh& tm_edges,
    EdgeToFaces& edge_to_faces)
  : tm_faces(tm_faces)
  , tm_edges(tm_edges)
  , edge_to_faces(edge_to_faces)
  {}

  void operator()( const Box& face_box, const Box& edge_box) const
  {
    halfedge_descriptor fh = face_box.info();
    halfedge_descriptor eh = edge_box.info();

    edge_to_faces[eh].insert(face(fh, tm_faces));
  }

  void operator()( const Box* face_box_ptr, const Box* edge_box_ptr) const
  {
    operator()(*face_box_ptr, *edge_box_ptr);
  }
};

template<class TriangleMesh,
         class VertexPointMap,
         class EdgeToFaces,
         class CoplanarFaceSet>
class Collect_face_bbox_per_edge_bbox_with_coplanar_handling {
protected:
  const TriangleMesh& tm_faces;
  const TriangleMesh& tm_edges;
  const VertexPointMap& vpmap_tmf;
  const VertexPointMap& vpmap_tme;
  EdgeToFaces& edge_to_faces;
  CoplanarFaceSet& coplanar_faces;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor> Box;
  typedef typename boost::property_traits<VertexPointMap>::reference Point;

public:
  Collect_face_bbox_per_edge_bbox_with_coplanar_handling(
    const TriangleMesh& tm_faces,
    const TriangleMesh& tm_edges,
    const VertexPointMap& vpmap_tmf,
    const VertexPointMap& vpmap_tme,
    EdgeToFaces& edge_to_faces,
    CoplanarFaceSet& coplanar_faces)
  : tm_faces(tm_faces)
  , tm_edges(tm_edges)
  , vpmap_tmf(vpmap_tmf)
  , vpmap_tme(vpmap_tme)
  , edge_to_faces(edge_to_faces)
  , coplanar_faces(coplanar_faces)
  {}

  void operator()( const Box& face_box, const Box& edge_box) const {
    halfedge_descriptor fh = face_box.info();
    halfedge_descriptor eh = edge_box.info();
    if(is_border(eh,tm_edges)) eh = opposite(eh, tm_edges);

    //check if the segment intersects the plane of the facet or if it is included in the plane
    Point a = get(vpmap_tmf, source(fh, tm_faces));
    Point b = get(vpmap_tmf, target(fh, tm_faces));
    Point c = get(vpmap_tmf, target(next(fh, tm_faces), tm_faces));
    /// SHOULD_USE_TRAITS_TAG
    const Orientation abcp = orientation(a,b,c, get(vpmap_tme, target(eh, tm_edges)));
    const Orientation abcq = orientation(a,b,c, get(vpmap_tme, source(eh, tm_edges)));
    if (abcp==abcq){
      if (abcp!=COPLANAR){
        return; //no intersection
      }

      if (orientation(a,b,c,get(vpmap_tme, target( next(eh, tm_edges), tm_edges)))==COPLANAR)
      {
        coplanar_faces.insert(
            &tm_edges < &tm_faces // TODO can we avoid by reporting them in only of the two calls to the filter function?
            ? std::make_pair(face(eh, tm_edges), face(fh, tm_faces))
            : std::make_pair(face(fh, tm_faces), face(eh, tm_edges))
          );
      }
      halfedge_descriptor eh_opp=opposite(eh, tm_edges);
      if (!is_border(eh_opp, tm_edges) &&
          orientation(a,b,c,get(vpmap_tme, target(next(eh_opp, tm_edges),tm_edges)))==COPLANAR)
      {
        coplanar_faces.insert(
            &tm_edges < &tm_faces // TODO can we avoid by reporting them in only of the two calls to the filter function?
            ? std::make_pair(face(opposite(eh, tm_edges), tm_edges), face(fh, tm_faces))
            : std::make_pair(face(fh, tm_faces), face(opposite(eh, tm_edges), tm_edges))
          );
      }
      //in case only the edge is coplanar, the intersection points will be detected using an incident facet
      return;
    }
    // non-coplanar case
    edge_to_faces[edge(eh,tm_edges)].insert(face(fh, tm_faces));
  }

  void operator()(const Box* face_box_ptr, const Box* edge_box_ptr) const
  {
    operator()(*face_box_ptr, *edge_box_ptr);
  }
};

template <class TriangleMesh, class Base>
class Callback_with_self_intersection_report
  : public Base
{
  typedef typename Base::face_descriptor face_descriptor;
  typedef typename Base::Box Box;
  boost::shared_ptr< std::set<face_descriptor> > faces_with_bbox_involved_in_intersections;
public:
  Callback_with_self_intersection_report(const Base& base)
  : Base(base), faces_with_bbox_involved_in_intersections(new std::set<face_descriptor>())
  {}

  void operator()( const Box* fb, const Box* eb) {
    faces_with_bbox_involved_in_intersections->insert( face(fb->info(), this->tm_faces) );
    Base::operator()(fb, eb);
  }
  bool self_intersections_found()
  {
    return Polygon_mesh_processing::does_self_intersect(
      *faces_with_bbox_involved_in_intersections,
      this->tm_faces,
      Polygon_mesh_processing::parameters::vertex_point_map(this->vpmap_tmf));
  }
};

} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_CALLBACK_H
