// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_REFINE_POLYHEDRON_3_H
#define CGAL_POLYGON_MESH_PROCESSING_REFINE_POLYHEDRON_3_H

#include <cmath>
#include <map>
#include <set>
#include <list>

#include <CGAL/assertions.h>
#include <CGAL/trace.h>
#include <CGAL/Timer.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/properties.h>
#include <boost/property_map/property_map.hpp>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

template<class PolygonMesh, class VertexPointMap>
class Refine_Polyhedron_3 {
//// typedefs
  typedef typename boost::property_traits<VertexPointMap>::value_type     Point_3;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor      face_descriptor;

  typedef Halfedge_around_face_circulator<PolygonMesh>    Halfedge_around_facet_circulator;
  typedef Halfedge_around_target_circulator<PolygonMesh>  Halfedge_around_vertex_circulator;

private:
  PolygonMesh& pmesh;
  VertexPointMap vpmap;
  
  bool flippable(halfedge_descriptor h) {
    // this check is added so that edge flip does not break manifoldness
    // it might happen when there is an edge where flip_edge(h) will be placed (i.e. two edges collide after flip)
    vertex_descriptor v_tip_0 = target(next(h,pmesh),pmesh);
    vertex_descriptor v_tip_1 = target(next(opposite(h,pmesh),pmesh),pmesh);
    Halfedge_around_vertex_circulator v_cir(next(h,pmesh), pmesh), v_end(v_cir);
    do {
      if(target(opposite(*v_cir, pmesh),pmesh) == v_tip_1) { return false; }
    } while(++v_cir != v_end);
    
    // also eliminate collinear triangle generation
    if( CGAL::collinear(vpmap[v_tip_0], vpmap[v_tip_1], vpmap[target(h, pmesh)]) ||
        CGAL::collinear(vpmap[v_tip_0], vpmap[v_tip_1], vpmap[target(opposite(h, pmesh),pmesh)]) ) {
      return false;
    }

    return true;
  }

  bool relax(halfedge_descriptor h)
  {
    const Point_3& p = vpmap[target(h,pmesh)];
    const Point_3& q = vpmap[target(opposite(h,pmesh),pmesh)];
    const Point_3& r = vpmap[target(next(h,pmesh),pmesh)];
    const Point_3& s = vpmap[target(next(opposite(h,pmesh),pmesh),pmesh)];
    if( (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,r,s)) ||
        (CGAL::ON_UNBOUNDED_SIDE  != CGAL::side_of_bounded_sphere(p,q,s,r)) )
    {
      if(flippable(h)) {
        Euler::flip_edge(h,pmesh);
        return true;
      }
    }
    return false;
  }

  template<class VertexOutputIterator,
           class FaceOutputIterator,
           class FaceRange>
  bool subdivide(const FaceRange& faces,
                 const std::set<halfedge_descriptor>& border_edges,
                 std::map<vertex_descriptor, double>& scale_attribute, 
                 VertexOutputIterator& vertex_out,
                 FaceOutputIterator& facet_out,
                 std::vector<face_descriptor>& new_faces,
                 double alpha)
  {
    BOOST_FOREACH(face_descriptor fd, faces)
    {
      CGAL_assertion(fd  != boost::graph_traits<PolygonMesh>::null_face());

      vertex_descriptor vi = target(halfedge(fd,pmesh),pmesh);
      vertex_descriptor vj = target(next(halfedge(fd,pmesh),pmesh),pmesh);
      vertex_descriptor vk = target(prev(halfedge(fd,pmesh),pmesh),pmesh);
      Point_3 c = CGAL::centroid(vpmap[vi], vpmap[vj], vpmap[vk]);
      double sac  = (scale_attribute[vi] + scale_attribute[vj] + scale_attribute[vk])/3.0;
      double dist_c_vi = std::sqrt(CGAL::squared_distance(c,vpmap[vi]));
      double dist_c_vj = std::sqrt(CGAL::squared_distance(c,vpmap[vj]));
      double dist_c_vk = std::sqrt(CGAL::squared_distance(c,vpmap[vk]));
      if((alpha * dist_c_vi > sac) &&
         (alpha * dist_c_vj > sac) &&
         (alpha * dist_c_vk > sac) &&
         (alpha * dist_c_vi > scale_attribute[vi]) &&
         (alpha * dist_c_vj > scale_attribute[vj]) &&
         (alpha * dist_c_vk > scale_attribute[vk]))
      {
        halfedge_descriptor h = Euler::add_center_vertex(halfedge(fd,pmesh),pmesh);
        put(vpmap, target(h,pmesh), c);
          scale_attribute[target(h,pmesh)] = sac;
          *vertex_out++ = target(h,pmesh);

          // collect 2 new facets for next round 
          face_descriptor h1 = face(opposite(next(h,pmesh),pmesh),pmesh);
          face_descriptor h2 = face(opposite(h,pmesh),pmesh);
          new_faces.push_back(h1); new_faces.push_back(h2);
          *facet_out++ = h1;    *facet_out++ = h2;
          // relax edges of the  patching mesh 
          halfedge_descriptor e_ij = prev(h,pmesh);
          halfedge_descriptor e_ik = next(opposite(h,pmesh),pmesh);
          halfedge_descriptor e_jk = prev(opposite(next(h,pmesh),pmesh),pmesh);

          if(border_edges.find(e_ij) == border_edges.end()){
            relax(e_ij);
          }
          if(border_edges.find(e_ik) == border_edges.end()){
            relax(e_ik);
          }
          if(border_edges.find(e_jk) == border_edges.end()){
            relax(e_jk);
          }
      }
    }
    return !new_faces.empty();
  }

  template<typename FaceRange>
  bool relax(const FaceRange& faces,
             const std::vector<face_descriptor>& new_faces,
             const std::set<halfedge_descriptor>& border_edges)
  {
    int flips = 0;
    std::list<halfedge_descriptor> interior_edges;
    std::set<halfedge_descriptor> included_map; 

    collect_interior_edges(faces, border_edges, interior_edges, included_map);
    collect_interior_edges(new_faces, border_edges, interior_edges, included_map);

    CGAL_TRACE_STREAM << "Test " << interior_edges.size() << " edges " << std::endl;
    //do not just use std::set (included_map) for iteration, the order effects the output (we like to make it deterministic)
    BOOST_FOREACH(halfedge_descriptor h, interior_edges)
    {
      if (relax(h)) {
        ++flips;
      }
    }

    CGAL_TRACE_STREAM << "|flips| = " << flips << std::endl;
    return flips > 0;
  }

  template<typename FaceRange>
  void collect_interior_edges(const FaceRange& faces,
        const std::set<halfedge_descriptor>& border_edges,
        std::list<halfedge_descriptor>& interior_edges,
        std::set<halfedge_descriptor>& included_map)
  {
    BOOST_FOREACH(face_descriptor fd, faces)
    {
      Halfedge_around_face_circulator<PolygonMesh> circ(halfedge(fd, pmesh), pmesh), done(circ);
      do {
        halfedge_descriptor h = *circ;
        if (border_edges.find(h) == border_edges.end()){
          // do not remove included_map and use if(&*h < &*oh) { interior_edges.push_back(h) } 
          // which will change the order of edges from run to run
          halfedge_descriptor oh = opposite(h, pmesh);
          halfedge_descriptor h_rep = (h < oh) ? h : oh; // AF: was &*h < &*oh
          if (included_map.insert(h_rep).second) {
            interior_edges.push_back(h_rep);
          }
        }
      } while (++circ != done);
    }
  }

  double average_length(vertex_descriptor vh,
                        const std::set<face_descriptor>& interior_map, 
                        bool accept_internal_facets)
  {
    const Point_3& vp = vpmap[vh]; 
    Halfedge_around_target_circulator<PolygonMesh> circ(halfedge(vh,pmesh),pmesh), done(circ);
    int deg = 0;
    double sum = 0;
    do {
      face_descriptor f(face(*circ,pmesh)), f_op(face(opposite(*circ,pmesh),pmesh));

      if(!accept_internal_facets) {
        if(interior_map.find(f) != interior_map.end() && interior_map.find(f_op) != interior_map.end())
        { continue; } // which means current edge is an interior edge and should not be included in scale attribute calculation
      }

      const Point_3& vq = vpmap[target(opposite(*circ,pmesh),pmesh)];
      sum += std::sqrt(CGAL::squared_distance(vp, vq));
      ++deg;
    } while(++circ != done);

    CGAL_assertion(deg != 0); // this might happen when accept_internal_facets = false but there is
    return sum/deg;
  }

  template<typename FaceRange>
  void calculate_scale_attribute(const FaceRange& faces,
                                 const std::set<face_descriptor>& interior_map,
                                 std::map<vertex_descriptor, double>& scale_attribute,
                                 bool accept_internal_facets) 
  {
    BOOST_FOREACH(face_descriptor fd, faces)
    {
      Halfedge_around_face_circulator<PolygonMesh> circ(halfedge(fd,pmesh),pmesh), done(circ);
      do {
        vertex_descriptor v = target(*circ,pmesh);
        std::pair<typename std::map<vertex_descriptor, double>::iterator, bool> v_insert 
          = scale_attribute.insert(std::make_pair(v, 0));
        if(!v_insert.second) { continue; } // already calculated
        v_insert.first->second = average_length(v, interior_map, accept_internal_facets);
      } while(++circ != done);
    }
  }

  template<typename FaceRange>
  bool contain_internal_facets(const FaceRange& faces,
                               const std::set<face_descriptor>& interior_map) const
  {
    BOOST_FOREACH(face_descriptor fd, faces)
    {
      Halfedge_around_face_circulator<PolygonMesh> circ(halfedge(fd,pmesh),pmesh), done(circ);
      do {
        Halfedge_around_target_circulator<PolygonMesh> circ_v(*circ,pmesh), done_v(circ_v);
        bool internal_v = true;
        do {
          face_descriptor f(face(*circ,pmesh)), f_op(face(opposite(*circ_v,pmesh),pmesh));

          if(interior_map.find(f) == interior_map.end() || interior_map.find(f_op) == interior_map.end()) {
            internal_v = false;
            break;
          } 
        } while(++circ_v != done_v);

        if(internal_v) { return true; }
      } while(++circ != done);
    }
    return false;
  }

public:
  Refine_Polyhedron_3(PolygonMesh& pmesh
                    , VertexPointMap vpmap_)
    : pmesh(pmesh)
    , vpmap(vpmap_)
  {}

  template<class FaceRange, class FacetOutputIterator, class VertexOutputIterator>
  void refine(const FaceRange& faces,
              FacetOutputIterator& facet_out,
              VertexOutputIterator& vertex_out,
              double alpha)
  {
      // do not use just std::set, the order effects the output (for the same input we want to get same output)
    std::set<face_descriptor> interior_map(boost::begin(faces), boost::end(faces));

    // store boundary edges - to be used in relax 
    std::set<halfedge_descriptor> border_edges;
    BOOST_FOREACH(face_descriptor f, faces)
    {
      Halfedge_around_face_circulator<PolygonMesh> circ(halfedge(f,pmesh),pmesh), done(circ);
      do {
        if(interior_map.find(face(opposite(*circ,pmesh),pmesh)) == interior_map.end()) {
          border_edges.insert(*circ);
        }
      } while(++circ != done);
    }

    // check whether there is any need to accept internal facets
    bool accept_internal_facets = contain_internal_facets(faces, interior_map);
    std::map<vertex_descriptor, double> scale_attribute;
    calculate_scale_attribute(faces, interior_map, scale_attribute, accept_internal_facets);

    std::vector<face_descriptor> new_faces;
    CGAL::Timer total_timer; total_timer.start();
    for(int i = 0; i < 10; ++i)
    {
      CGAL::Timer timer; timer.start();
      bool is_subdivided = subdivide(faces, border_edges, scale_attribute, vertex_out, facet_out, new_faces, alpha);
      CGAL_TRACE_STREAM << "**Timer** subdivide() :" << timer.time() << std::endl; timer.reset();
      if(!is_subdivided) { 
        break; }

      bool is_relaxed = relax(faces, new_faces, border_edges);
      CGAL_TRACE_STREAM << "**Timer** relax() :" << timer.time() << std::endl;
      if(!is_relaxed) { 
        break; }
    }

    CGAL_TRACE_STREAM << "**Timer** TOTAL: " << total_timer.time() << std::endl;
  }

}; //end class Refine_Polyhedron_3

}//namespace internal

}//namespace Polygon_mesh_processing

}//namespace CGAL

#endif //CGAL_POLYGON_MESH_PROCESSING_REFINE_POLYHEDRON_3_H
