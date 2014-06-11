// Copyright (c) 2013 GeometryFactory (France).
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
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_POLYHEDRON_SLICER_3_H
#define CGAL_POLYHEDRON_SLICER_3_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/intersection_of_Polyhedra_3.h>

#include <CGAL/Timer.h>
#include <CGAL/Profile_counter.h>

#include <vector>
#include <map>
#include <deque>

#include <boost/tuple/tuple.hpp>
#include <boost/optional.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {

template<class Polyhedron, class Kernel>
class Polyhedron_slicer_3
{
private:
  typedef AABB_halfedge_graph_segment_primitive<Polyhedron>         AABB_primitive;
  typedef AABB_traits<Kernel, AABB_primitive>                       AABB_traits_;
  typedef AABB_tree<AABB_traits_>                                   AABB_tree_;

  typedef typename AABB_tree_::Object_and_primitive_id             Object_and_primitive_id;
  typedef typename AABB_tree_::Primitive_id                        Primitive_id;

  typedef typename Kernel::Plane_3    Plane;
  typedef typename Kernel::Segment_3  Segment;
  typedef typename Kernel::Point_3    Point;

  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor   Edge_const_handle;
  typedef typename boost::graph_traits<Polyhedron>::edge_iterator   Edge_const_iterator;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor Halfedge_const_handle;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor    Facet_const_handle;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor   Vertex_const_handle;
  typedef Halfedge_around_target_circulator<Polyhedron> Halfedge_around_vertex_const_circulator;
  typedef Halfedge_around_face_circulator<Polyhedron>  Halfedge_around_facet_const_circulator;

  // to unite halfedges under an "edge"
  struct Edge_comparator {
    bool operator()(Halfedge_const_handle h1, Halfedge_const_handle h2) const {
      return (std::min)(&*h1, &*(h1->opposite())) 
        < (std::min)(&*h2, &*(h2->opposite()));
    }
  };

  ////////////////////////////////////////////////
  // to represent a intersection point with plane
  struct Node {
    Point point;
  };

  struct Node_edge { 
    bool is_processed;
    Node_edge() : is_processed(false) 
    { }
  };

  // out_edges is setS for preventing parallel edges automatically (do not change it, or add custom check before add_edge)
  typedef boost::adjacency_list<
    boost::setS, boost::vecS, boost::undirectedS,
    Node, Node_edge>
  Node_graph;

  typedef typename Node_graph::vertex_descriptor Vertex_g;
  typedef typename Node_graph::vertex_iterator   Vertex_iterator_g;
  typedef typename Node_graph::edge_descriptor   Edge_g;
  typedef typename Node_graph::out_edge_iterator Out_edge_iterator_g;
  ////////////////////////////////////////////////

  // to store intersections of an edge
  // there is two nodes (intersection points) when an edge is coplanar to plane
  struct Node_pair {
    Vertex_g v1, v2;
    int vertex_count;

    Node_pair() : vertex_count(0){ }
    Node_pair(Vertex_g v1, Vertex_g v2) : v1(v1), v2(v2), vertex_count(2) { }

    void put(Vertex_g v) {
      CGAL_assertion(vertex_count <= 2);

      (vertex_count == 0 ? v1 : v2) = v;
      ++vertex_count;
    }
  };

  typedef std::map<Halfedge_const_handle, Node_pair, Edge_comparator> Edge_node_map;
  typedef typename Edge_node_map::iterator Edge_node_map_iterator;

  enum Intersection_type { BOUNDARY, INTERVAL, PLANAR };

  typedef std::map<Halfedge_const_handle, Intersection_type, Edge_comparator> Edge_intersection_map;
  typedef typename Edge_intersection_map::iterator Edge_intersection_map_iterator;

  // member variables //
  typename Kernel::Intersect_3 intersect_3_functor;
  AABB_tree_ tree;
  mutable Node_graph node_graph;
  Polyhedron& polyhedron;

  boost::tuple<Point, Intersection_type, Vertex_const_handle>
  halfedge_intersection(Halfedge_const_handle hf, const Plane& plane) const
  { 
    boost::tuple<Point, Intersection_type, Vertex_const_handle> ret;

    const Point& s = hf->vertex()->point();
    const Point& t = hf->opposite()->vertex()->point();
    Oriented_side s_os, t_os;
    s_os = plane.oriented_side(s);
    t_os = plane.oriented_side(t);

    // in case of planar just Intersection_type is not empty in tuple
    if(t_os == ON_ORIENTED_BOUNDARY && s_os == ON_ORIENTED_BOUNDARY) {
      ret.template get<1>() = PLANAR;
      return ret;
    }
    // in case of boundary Intersection_type and vertex are not empty in tuple
    if(t_os == ON_ORIENTED_BOUNDARY || s_os == ON_ORIENTED_BOUNDARY) { 
      ret.template get<1>() = BOUNDARY;
      ret.template get<2>() = t_os == ON_ORIENTED_BOUNDARY ? hf->opposite()->vertex() :
                                                    hf->vertex();
      return ret;
    }

    CGAL_assertion(t_os != s_os); // should one positive one negative
    // in case of interval Intersection_type and point are not empty in tuple
    ret.template get<1>() = INTERVAL;
    Object intersection = intersect_3_functor(plane, Segment(s,t));

    if(const Point* i_point  = object_cast<Point>(&intersection)) { 
      ret.template get<0>() = *i_point;
    }
    else if( object_cast<Segment>(&intersection) ) {
      // is it possible to predicate not-planar but construct segment ?
      CGAL_warning(!"on interval case - predicate not-planar but construct segment");
      ret.template get<1>() = PLANAR;
    }
    else {
      // prediction indicates intersection but construction returns nothing, returning closest point
      CGAL_warning(!"on interval case - no intersection found");
      ret.template get<0>() = squared_distance(plane, s) < squared_distance(plane, t) ? s : t; 
    }
    return ret;
  }

  void add_intersection_node_to_one_ring(Vertex_const_handle v, Edge_node_map& edge_node_map) const {
    Vertex_g v_g = boost::add_vertex(node_graph);
    node_graph[v_g].point = v->point();

    Halfedge_around_vertex_const_circulator around_vertex_c(v,polyhedron), done(around_vertex_c);
    do {
      edge_node_map[*around_vertex_c].put(v_g);
    } 
    while(++around_vertex_c != done);
  }
  
  void add_intersection_edge_to_facet_neighbors(Halfedge_const_handle hf, Edge_node_map& edge_node_map) const {
    Node_pair& hf_node_pair = edge_node_map.find(hf)->second;
    CGAL_assertion(hf_node_pair.vertex_count == 1);

    for(int i = 0; i < 2; ++i) // loop for hf and hf->opposite
    {
      if(i == 1) { hf = hf->opposite(); }
      if(hf->is_border()) { continue;} 

      for(int i = 0; i < 2; ++i) 
      {
        Halfedge_const_handle facet_hf = i == 0 ? hf->next() : hf->prev();
        Edge_node_map_iterator facet_hf_it = edge_node_map.find(facet_hf);

        if(facet_hf_it != edge_node_map.end()) {
          Node_pair& facet_hf_node_pair = facet_hf_it->second;
          CGAL_assertion(facet_hf_node_pair.vertex_count == 1); // == 2 if hf is planar and that cant happen 
          boost::add_edge(hf_node_pair.v1, facet_hf_node_pair.v1, node_graph); 
        }
      }
    } // for hf and hf->opposite
  }

  void add_intersection_edge_to_vertex_neighbors(Halfedge_const_handle hf, Vertex_const_handle v, Edge_node_map& edge_node_map) const {
    // do not worry about duplicate edges, they are not allowed (but performance might be a concern)
    Node_pair& hf_node_pair = edge_node_map.find(hf)->second;
    Vertex_g v_g = node_graph[hf_node_pair.v1].point == v->point() ? hf_node_pair.v1 : hf_node_pair.v2;// find node containing v

    Halfedge_around_vertex_const_circulator around_vertex_c(v, polyhedron), done(around_vertex_c);
    do {
      if((*around_vertex_c)->is_border()) { continue;} 
      Node_pair& around_vertex_node_pair = edge_node_map.find(*around_vertex_c)->second;
      Halfedge_around_facet_const_circulator around_facet_c(*around_vertex_c,polyhedron), done2(around_facet_c);
      do {
        CGAL_assertion(around_vertex_node_pair.vertex_count != 0);
        if(around_vertex_node_pair.v1 != v_g) {
          boost::add_edge(around_vertex_node_pair.v1, v_g, node_graph); 
        }
        if(around_vertex_node_pair.vertex_count == 2 && around_vertex_node_pair.v2 != v_g) {
          boost::add_edge(around_vertex_node_pair.v2, v_g, node_graph); 
        }
      } 
      while(++around_facet_c != done2);
    } 
    while(++around_vertex_c != done);
  }

  template<class OutputIterator>
  OutputIterator intersect_plane(const Plane& plane, OutputIterator out) const
  {
    node_graph.clear();

    // find out intersecting halfedges (note that tree contains edges only with custom comparator)
    std::vector<Edge_const_handle> intersected_edges;
    tree.all_intersected_primitives(plane, std::back_inserter(intersected_edges));

    // create node graph from segments
    // each node is associated with multiple edges
    Edge_node_map edge_node_map;
    Edge_intersection_map edge_intersection_map;
    for(typename std::vector<Edge_const_handle>::iterator it = intersected_edges.begin();
      it != intersected_edges.end(); ++it)
    {
      Halfedge_const_handle hf = halfedge(*it,polyhedron);
      Node_pair& assoc_nodes = edge_node_map[hf];
      CGAL_assertion(assoc_nodes.vertex_count < 3); // every Node_pair can at most contain 2 nodes

      boost::tuple<Point, Intersection_type, Vertex_const_handle> intersection = halfedge_intersection(hf, plane);
      edge_intersection_map[hf] = intersection.template get<1>(); // save intersection type

      if(intersection.template get<1>() == INTERVAL) {
        CGAL_PROFILER(" number of INTERVAL intersections.");
        CGAL_assertion(assoc_nodes.vertex_count == 0); 

        Vertex_g v = boost::add_vertex(node_graph);
        node_graph[v].point = intersection.template get<0>();
        assoc_nodes.put(v);
        // no related hf to edge node, done
      }
      else if(intersection.template get<1>() == BOUNDARY) {
        CGAL_PROFILER(" number of BOUNDARY intersections.");
        CGAL_assertion(assoc_nodes.vertex_count < 2); 

        if(assoc_nodes.vertex_count == 0) {
          Vertex_const_handle v_h = intersection.template get<2>(); // boundary vertex
          // update edge_node_map for neighbor halfedges of v_h
          add_intersection_node_to_one_ring(v_h, edge_node_map);
        } // else node is already added by other hf
      }
      else {// intersection.get<1>() == PLANAR
        CGAL_PROFILER(" number of PLANAR intersections.");
        CGAL_assertion(intersection.template get<1>() == PLANAR);

        if(assoc_nodes.vertex_count != 2) {
          if(assoc_nodes.vertex_count == 1) { // there is one intersection node that we need to add
            if(node_graph[assoc_nodes.v1].point == hf->vertex()->point()) {
              add_intersection_node_to_one_ring(hf->opposite()->vertex(), edge_node_map);
            }
            else {
              CGAL_assertion(node_graph[assoc_nodes.v1].point == hf->opposite()->vertex()->point());
              add_intersection_node_to_one_ring(hf->vertex(), edge_node_map);
            }
          }
          else { // assoc_nodes.vertex_count == 0
            // update edge_node_map for neighbor halfedges 
            add_intersection_node_to_one_ring(hf->vertex(), edge_node_map);
            add_intersection_node_to_one_ring(hf->opposite()->vertex(), edge_node_map);
          }
        } // else both nodes are already added by other hfs
      }
    } // for(typename std::vector<Halfedge_const_handle>::iterator it = intersected_edges.begin()...
    
    // introduce node connectivity
    for(typename std::vector<Edge_const_handle>::iterator it = intersected_edges.begin();
      it != intersected_edges.end(); ++it)
    {
      Halfedge_const_handle hf = halfedge(*it,polyhedron);
      Edge_intersection_map_iterator intersection_it = edge_intersection_map.find(hf);
      CGAL_assertion(intersection_it != edge_intersection_map.end());

      Intersection_type intersection = intersection_it->second;
      if(intersection == INTERVAL) {
        add_intersection_edge_to_facet_neighbors(hf, edge_node_map);
      }
      else if(intersection == BOUNDARY) {
        Node_pair& node_pair = edge_node_map[hf];
        Vertex_const_handle v = hf->vertex()->point() == node_graph[node_pair.v1].point ?
          hf->vertex() : hf->opposite()->vertex();
        add_intersection_edge_to_vertex_neighbors(hf, v, edge_node_map);
      }
      else {
        add_intersection_edge_to_vertex_neighbors(hf, hf->vertex(), edge_node_map);
        add_intersection_edge_to_vertex_neighbors(hf, hf->opposite()->vertex(), edge_node_map);
      }
    }
    
    // use node_graph to construct polylines
    return construct_polylines(out);
  }

  // find a new node to advance, if no proper node exists return v
  std::pair<Vertex_g, bool> find_next(Vertex_g v) const {
    std::pair<Vertex_g, bool> ret;
    Out_edge_iterator_g ei_b, ei_e;
    for(boost::tie(ei_b, ei_e) = boost::out_edges(v, node_graph); ei_b != ei_e; ++ei_b) {
      if(node_graph[*ei_b].is_processed) { continue; } // do not go over same edge twice

      ret.first = boost::target(*ei_b, node_graph);
      ret.second = true;
      node_graph[*ei_b].is_processed = true;
      
      CGAL_assertion(ret.first != v);
      return ret;
    }

    ret.second = false;
    return ret;
  }

  template<class OutputIterator>
  OutputIterator construct_polylines(OutputIterator out) const  {
    //std::cout << boost::num_vertices(node_graph) << std::endl;
    //std::cout << boost::num_edges(node_graph) << std::endl;
    //Vertex_iterator_g v_b, v_e;
    //for(boost::tie(v_b, v_e) = boost::vertices(node_graph); v_b != v_e; ++v_b)
    //{
    //  Out_edge_iterator_g e_b, e_e;
    //  boost::tie(e_b, e_e) = boost::out_edges(*v_b, node_graph);
    //  for( ;e_b != e_e; ++e_b) // keep current vertex active while there is 
    //  {
    //    std::cout << "source " << boost::source(*e_b, node_graph) 
    //      << " target " << boost::target(*e_b, node_graph) << std::endl;
    //  }
    //}

    Vertex_iterator_g v_b, v_e;
    for(boost::tie(v_b, v_e) = boost::vertices(node_graph); v_b != v_e; ++v_b)
    {
      Vertex_g v = *v_b;
      Out_edge_iterator_g e_b, e_e;
      boost::tie(e_b, e_e) = boost::out_edges(v, node_graph);
      // isolated intersection, construct a polyline with one point
      if(e_b == e_e) {
        std::vector<Point> polyline;
        polyline.push_back(node_graph[v].point);
        *out++ = polyline;
        continue;
      }

      Vertex_g v_head = v;
      // there are edges, construct one or more polylines
      for( ;e_b != e_e; ++e_b) // keep current vertex active while there is an non-processed edge to go
      {
        if(node_graph[*e_b].is_processed) { continue; } // we previously put it to polylines

        // construct new polyline
        std::vector<Point> polyline;

        bool first_border = true;
        while(true) {
          polyline.push_back(node_graph[v].point);
          bool found;
          boost::tie(v, found) = find_next(v);

          if(!found) { // intersection at boundary
            if(!first_border) { break; } // completed non-loop polyline

            first_border = false;
            boost::tie(v, found) = find_next(v_head); // try to go other direction
            if(!found) { break; }
            std::reverse(polyline.begin(), polyline.end());
            continue;
          }

          if(v == v_head) { // loop is completed
            polyline.push_back(node_graph[v_head].point); // put first point again
            break;
          }
        } // while( true )...
        *out++ = polyline;

      } // for( ;e_b != e_e;...
    } // for(boost::tie(v_b, v_e) = boost::vertices(node_graph)...
    return out;
  }

public:

  /**
  * Constructor. `polyhedron` must be valid polyhedron as long as this functor is used.
  * @param polyhedron the polyhedron to be cut
  * @param kernel the kernel
  */
  Polyhedron_slicer_3(const Polyhedron& polyhedron, const Kernel& kernel = Kernel())
  : intersect_3_functor(kernel.intersect_3_object()),
    tree( edges(polyhedron).first,
          edges(polyhedron).second,
          polyhedron),
    polyhedron(const_cast<Polyhedron&>(polyhedron))
  { }

  /**
   * @tparam OutputIterator an output iterator accepting polylines. A polyline is considered to be a `std::vector<Kernel::Point_3>`. A polyline is closed if its first and last points are identical.
   * @param plane the plane to intersect the polyhedron with
   * @out output iterator of polylines
   * computes the intersection polylines of the polyhedron passed in the constructor with `plane` and puts each of them in `out`
   */
  template <class OutputIterator>
  OutputIterator operator() (const typename Kernel::Plane_3& plane, OutputIterator out) const {
    CGAL_precondition(!plane.is_degenerate());
    return intersect_plane(plane, out);
  }
};

}// end of namespace CGAL
#endif //CGAL_POLYHEDRON_SLICER_3_H
