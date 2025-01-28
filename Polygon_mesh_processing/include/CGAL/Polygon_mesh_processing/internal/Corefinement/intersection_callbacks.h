// Copyright (c) 2016 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_CALLBACK_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_CALLBACK_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Box_intersection_d/Box_with_info_d.h>
#include <CGAL/property_map.h>
#include <CGAL/enum.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <boost/graph/graph_traits.hpp>
#include <memory>
#include <set>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace Corefinement {

struct Self_intersection_exception :
  public std::runtime_error
{
  Self_intersection_exception()
    : std::runtime_error("Self-intersection detected in input mesh")
  {}
};

template<class TriangleMesh, class EdgeToFaces>
class Collect_face_bbox_per_edge_bbox {
protected:
  const TriangleMesh& tm_faces;
  const TriangleMesh& tm_edges;
  EdgeToFaces& edge_to_faces;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor, Box_policy> Box;

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

    edge_to_faces[edge(eh,tm_edges)].insert(face(fh, tm_faces));
  }

  void operator()( const Box* face_box_ptr, const Box* edge_box_ptr) const
  {
    operator()(*face_box_ptr, *edge_box_ptr);
  }
};

template<class TriangleMesh,
         class VertexPointMapF, class VertexPointMapE,
         class EdgeToFaces,
         class CoplanarFaceSet,
         class Visitor>
class Collect_face_bbox_per_edge_bbox_with_coplanar_handling {
protected:
  const TriangleMesh& tm_faces;
  const TriangleMesh& tm_edges;
  const VertexPointMapF& vpmap_tmf;
  const VertexPointMapE& vpmap_tme;
  EdgeToFaces& edge_to_faces;
  CoplanarFaceSet& coplanar_faces;
  const Visitor& visitor;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_traits<VertexPointMapF>::reference Point;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor, Box_policy> Box;

public:
  Collect_face_bbox_per_edge_bbox_with_coplanar_handling(
    const TriangleMesh& tm_faces,
    const TriangleMesh& tm_edges,
    const VertexPointMapF& vpmap_tmf,
    const VertexPointMapE& vpmap_tme,
    EdgeToFaces& edge_to_faces,
    CoplanarFaceSet& coplanar_faces,
    const Visitor& visitor)
  : tm_faces(tm_faces)
  , tm_edges(tm_edges)
  , vpmap_tmf(vpmap_tmf)
  , vpmap_tme(vpmap_tme)
  , edge_to_faces(edge_to_faces)
  , coplanar_faces(coplanar_faces)
  , visitor(visitor)
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

  bool is_face_degenerated(halfedge_descriptor fh) const
  {
    Point a = get(vpmap_tmf, source(fh, tm_faces));
    Point b = get(vpmap_tmf, target(fh, tm_faces));
    Point c = get(vpmap_tmf, target(next(fh, tm_faces), tm_faces));

    return collinear(a, b, c);
  }

  bool are_edge_faces_degenerated(halfedge_descriptor eh) const
  {
    Point a = get(vpmap_tme, source(eh, tm_edges));
    Point b = get(vpmap_tme, target(eh, tm_edges));

    if(!is_border(eh,tm_edges))
    {
      Point c = get(vpmap_tme, target(next(eh, tm_edges), tm_edges));
      if (collinear(a, b, c)) return true;
    }

    eh = opposite(eh, tm_edges);
    if(!is_border(eh,tm_edges))
    {
      Point c = get(vpmap_tme, target(next(eh, tm_edges), tm_edges));
      if (collinear(a, b, c)) return true;
    }

    return false;
  }


  void operator()(const Box* face_box_ptr, const Box* edge_box_ptr) const
  {
    operator()(*face_box_ptr, *edge_box_ptr);
  }

  void progress(double d)
  {
    visitor.progress_filtering_intersections(d);
  }

};

template<class TriangleMesh,
         class VertexPointMap,
         class EdgeToFaces,
         class CoplanarFaceSet>
class Collect_face_bbox_per_edge_bbox_with_coplanar_handling_one_mesh {
protected:
  const TriangleMesh& tm;
  const VertexPointMap& vpmap;
  EdgeToFaces& edge_to_faces;
  CoplanarFaceSet& coplanar_faces;

  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, halfedge_descriptor, Box_policy> Box;

  typedef typename boost::property_traits<VertexPointMap>::reference Point;

  bool is_edge_target_incident_to_face(halfedge_descriptor hd,
                                       halfedge_descriptor fhd) const
  {
    vertex_descriptor vd = target(hd,tm);
    for (int i=0;i<3;++i)
      if (target(fhd, tm) == vd)
        return true;
      else
        fhd = next(fhd, tm);

    return false;
  }

  bool does_incident_edge_intersect_face_interior(halfedge_descriptor hd,
                                                  halfedge_descriptor fhd) const
  {
    vertex_descriptor tgt = target(hd,tm);
    vertex_descriptor src = source(hd,tm);
    for (int i=0;i<3;++i)
      if (target(fhd, tm) == tgt)
        break;
      else
        if (target(fhd, tm) == src)
        {
          hd=opposite(hd, tm);
          break;
        }
        else
          fhd = next(fhd, tm);

    CGAL_assertion(target(fhd,tm)==target(hd, tm));

    Point a = get(vpmap, target(fhd, tm));
    Point b = get(vpmap, target(next(fhd, tm), tm));
    Point c = get(vpmap, source(fhd, tm));
    Point q = get(vpmap, source(hd, tm));

    if ( coplanar_orientation(a, b, c, q) == POSITIVE &&
         coplanar_orientation(c, a, b, q) == POSITIVE )
    {
      return true;
    }

    return false;
  }

public:
  Collect_face_bbox_per_edge_bbox_with_coplanar_handling_one_mesh(
    const TriangleMesh& tm,
    const VertexPointMap& vpmap,
    EdgeToFaces& edge_to_faces,
    CoplanarFaceSet& coplanar_faces)
  : tm(tm)
  , vpmap(vpmap)
  , edge_to_faces(edge_to_faces)
  , coplanar_faces(coplanar_faces)
  {}

  void operator()( const Box& face_box, const Box& edge_box) const {
    halfedge_descriptor fh = face_box.info();
    halfedge_descriptor eh = edge_box.info();

    if ( face(eh, tm) == face(fh, tm) || face(opposite(eh,tm), tm) == face(fh, tm) )
      return; //edge incident to the triangle

    if(is_border(eh,tm)) eh = opposite(eh, tm);
    halfedge_descriptor eh_opp=opposite(eh, tm);

    //check if the segment intersects the plane of the facet or if it is included in the plane
    Point a = get(vpmap, source(fh, tm));
    Point b = get(vpmap, target(fh, tm));
    Point c = get(vpmap, target(next(fh, tm), tm));
    /// SHOULD_USE_TRAITS_TAG
    const Orientation abcp = orientation(a,b,c, get(vpmap, target(eh, tm)));
    const Orientation abcq = orientation(a,b,c, get(vpmap, source(eh, tm)));
    if (abcp==abcq){
      if (abcp!=COPLANAR){
        return; //no intersection
      }

      if ( is_edge_target_incident_to_face(eh, fh) ||
           is_edge_target_incident_to_face(eh_opp, fh) )
      {
        if (does_incident_edge_intersect_face_interior(eh,fh))
        {
          if (orientation(a,b,c,get(vpmap, target( next(eh, tm), tm)))==COPLANAR)
          {
            coplanar_faces.insert(make_sorted_pair(face(eh, tm), face(fh, tm)));
          }

          if (!is_border(eh_opp, tm) &&
              orientation(a,b,c,get(vpmap, target(next(eh_opp, tm),tm)))==COPLANAR)
          {
            coplanar_faces.insert(make_sorted_pair(face(eh_opp, tm), face(fh, tm)));
          }
        }
        else
          return; // If there is an intersection, it will be reported
                  // when handling another edge of the face
      }
      else
      {
        // In the following we differ the insertion of intersecting coplanar
        // faces if the third vertex is a vertex of the face so that it is
        // handled in the previous case (when an incident edge will be handled)
        if (!is_edge_target_incident_to_face(next(eh, tm), fh) &&
            orientation(a,b,c,get(vpmap, target( next(eh, tm), tm)))==COPLANAR)
        {
          coplanar_faces.insert(make_sorted_pair(face(eh, tm), face(fh, tm)));
        }

        if (!is_border(eh_opp, tm) &&
            !is_edge_target_incident_to_face(next(eh_opp, tm), fh) &&
            orientation(a,b,c,get(vpmap, target(next(eh_opp, tm),tm)))==COPLANAR)
        {
          coplanar_faces.insert(make_sorted_pair(face(eh_opp, tm), face(fh, tm)));
        }
      }

      // In case only the edge is coplanar, the intersection points will
      // be detected using an incident facet
      return;
    }

    if ( abcp==COPLANAR &&
         is_edge_target_incident_to_face(eh, fh) )
    {
      return; // no intersection (incident edge)
    }

    if ( abcq==COPLANAR &&
         is_edge_target_incident_to_face(opposite(eh, tm), fh) )
    {
      return; // no intersection (incident edge)
    }

    // non-coplanar case
    edge_to_faces[edge(eh,tm)].insert(face(fh, tm));
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
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::Box Box;
  std::set<face_descriptor>* tmf_collected_faces_ptr;
  std::set<face_descriptor>* tme_collected_faces_ptr;
public:
  Callback_with_self_intersection_report(const Base& base,
                                         std::set<face_descriptor>& tmf_collected_faces,
                                         std::set<face_descriptor>& tme_collected_faces)
  : Base(base),
    tmf_collected_faces_ptr(&tmf_collected_faces),
    tme_collected_faces_ptr(&tme_collected_faces)
  {}

  void operator()( const Box* fb, const Box* eb) {
    halfedge_descriptor h = eb->info();
    if (!is_border(h, this->tm_edges))
      tme_collected_faces_ptr->insert( face(h, this->tm_edges) );
    h = opposite(h, this->tm_edges);
    if (!is_border(h, this->tm_edges))
      tme_collected_faces_ptr->insert( face(h, this->tm_edges) );
    tmf_collected_faces_ptr->insert( face(fb->info(), this->tm_faces) );

    // throw if one of the faces are degenerated
    if (this->is_face_degenerated(fb->info()) ||
        this->are_edge_faces_degenerated(h))
    {
      throw Self_intersection_exception();
    }

    Base::operator()(fb, eb);
  }
  bool self_intersections_found()
  {
    return
      Polygon_mesh_processing::does_self_intersect(
        *tmf_collected_faces_ptr,
        this->tm_faces,
        CGAL::parameters::vertex_point_map(this->vpmap_tmf))
    ||
      Polygon_mesh_processing::does_self_intersect(
        *tme_collected_faces_ptr,
        this->tm_edges,
        CGAL::parameters::vertex_point_map(this->vpmap_tme));
  }
};

} } } // CGAL::Polygon_mesh_processing::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_CALLBACK_H
