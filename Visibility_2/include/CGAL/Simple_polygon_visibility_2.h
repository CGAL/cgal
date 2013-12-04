// Copyright (c) 2013 Technical University Braunschweig (Germany).
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
// $URL$`
// $Id$
//
//
// Author(s):  Francisc Bungiu <fbungiu@gmail.com>
//             Michael Hemmer <michael.hemmer@cgal.org>

#ifndef CGAL_SIMPLE_POLYGON_VISIBILITY_2_H
#define CGAL_SIMPLE_POLYGON_VISIBILITY_2_H

#include <CGAL/Arrangement_2.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/Visibility_2/visibility_utils.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <stack>
#include <map>

namespace CGAL {

template<class Arrangement_2_, class Visibility_arrangement_2_ = Arrangement_2_, class RegularizationTag = CGAL::Tag_true> 
class Simple_polygon_visibility_2 {

public:
  // Currently only consider with same type for both
  typedef Arrangement_2_                                Arrangement_2;
  typedef Visibility_arrangement_2_                     Visibility_arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            K;

  typedef typename Arrangement_2::Halfedge_const_handle       
                                                        Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
                                        Halfedge_around_vertex_const_circulator;

  typedef typename Geometry_traits_2::Point_2           Point_2;
  typedef typename Geometry_traits_2::Ray_2             Ray_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Line_2            Line_2;
  typedef typename Geometry_traits_2::Vector_2          Vector_2;
  typedef typename Geometry_traits_2::Direction_2       Direction_2;
  typedef typename Geometry_traits_2::FT                Number_type;
  typedef typename Geometry_traits_2::Object_2          Object_2;

  typedef RegularizationTag                       Regularization_tag;
  typedef CGAL::Tag_false                         Supports_general_polygon_tag;
  typedef CGAL::Tag_true                          Supports_simple_polygon_tag;                                        

  Simple_polygon_visibility_2() : p_arr(NULL), geom_traits(NULL) {};

  /*! Constructor given an arrangement and the Regularization tag. */
  Simple_polygon_visibility_2(const Arrangement_2& arr): 
    p_arr(&arr) {
    geom_traits = p_arr->geometry_traits();
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
  }

  /*! Method to check if the visibility object is attached or not to
      an arrangement*/
  bool is_attached() {
    return (p_arr != NULL);
  }

  /*! Attaches the visibility object to the 'arr' arrangement */
  void attach(const Arrangement_2& arr) {
    p_arr = &arr;
    geom_traits = p_arr->geometry_traits();
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
  }

  /*! Detaches the visibility object from the arrangement it is
      attached to*/
  void detach() {
    p_arr = NULL;
    geom_traits = NULL;
    vertices.clear();
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;
    p_cdt = boost::shared_ptr<CDT>();
  }

  /*! Getter method for the input arrangement*/
  const Arrangement_2& arr() {
    return *p_arr;
  }

  /*! Computes the visibility object from the query point 'q' in the face
      'face' and constructs the output in 'out_arr'*/
  typename Visibility_arrangement_2::Face_handle 
  compute_visibility(const Point_2& q, 
      const Face_const_handle face,
      Visibility_arrangement_2& out_arr) {
    
    assert(query_pt_is_vertex == false);
    assert(query_pt_is_on_halfedge == false);

    // Now retrieve the circulator to first visible vertex from triangulation
    Ccb_halfedge_const_circulator circ = find_visible_start(face, q);
    Ccb_halfedge_const_circulator curr = circ;
    Halfedge_const_handle he;

    do {
      he = curr;
      vertices.push_back(he->source()->point());
    } while(++curr != circ);

    vertices.push_back(vertices[0]);

    visibility_region_impl(q);

    typename std::vector<Point_2> points;
    while (!s.empty()) {
      Point_2 curr_point = s.top();
      points.push_back(curr_point);
      s.pop();
    }

    std::reverse(points.begin(), points.end());

    CGAL::Visibility_2::report_while_handling_needles
                              <Simple_polygon_visibility_2>(geom_traits, 
                                                            q, 
                                                            points,                                 
                                                            out_arr);  

    CGAL_precondition(out_arr.number_of_isolated_vertices() == 0);
    CGAL_precondition(s.size() == 0);
    conditional_regularize(out_arr, Regularization_tag());
    vertices.clear();
    if (out_arr.faces_begin()->is_unbounded()) {
      return ++out_arr.faces_begin();
    }
    else {
      return out_arr.faces_begin();
    }
  }

  /*! Computes the visibility region of the query point 'q' located on the
      halfedge 'he' and constructs the output in 'out_arr'*/
  typename Visibility_arrangement_2::Face_handle 
  compute_visibility(
      const Point_2& q, 
      const Halfedge_const_handle he,
      Visibility_arrangement_2& out_arr ) 
  {
    query_pt_is_vertex = false;
    query_pt_is_on_halfedge = false;

    if (q != he->source()->point()) {
      if (q != he->target()->point()) {
        vertices.push_back(q);
        vertices.push_back(he->target()->point());
        query_pt_is_on_halfedge = true;
      }
      else {
        vertices.push_back(q);
        query_pt_is_vertex = true;
      }
    }

    Ccb_halfedge_const_circulator circ = he;
    circ++;
    Ccb_halfedge_const_circulator curr = circ;

    do {
      Halfedge_const_handle he_handle = curr;
      Point_2 curr_vertex = he_handle->target()->point();
      vertices.push_back(curr_vertex);
    } while (++curr != circ);

    vertices.pop_back();
    vertices.push_back(vertices[0]);

    visibility_region_impl(q);
    
    typename std::vector<Point_2> points;
    if (!s.empty()) {
      Point_2 prev_pt = s.top();
      if (prev_pt != q) {
        points.push_back(prev_pt);
      }
      else if (query_pt_is_vertex) {
        points.push_back(prev_pt); 
      }
      if (!s.empty()) {
        s.pop();
      }
      while(!s.empty()) {
        Point_2 curr_pt = s.top();
        if (curr_pt != q) {
          points.push_back(curr_pt);
        }
        else if (query_pt_is_vertex) {
          points.push_back(curr_pt); 
        }
        s.pop();
      }
    }

    std::reverse(points.begin(), points.end());

    CGAL::Visibility_2::report_while_handling_needles
                              <Simple_polygon_visibility_2>(geom_traits, 
                                                            q, 
                                                            points,                                
                                                            out_arr);
    CGAL_precondition(out_arr.number_of_isolated_vertices() == 0);
    CGAL_precondition(s.size() == 0);
    conditional_regularize(out_arr, Regularization_tag());
    vertices.clear();
    if (out_arr.faces_begin()->is_unbounded()) {
      return ++out_arr.faces_begin();
    }
    else {
      return out_arr.faces_begin();
    }
  }

private:
  typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::No_intersection_tag                                Itag;
  typedef CGAL::Constrained_triangulation_2<K, TDS, Itag> CDT;

private:
  const Arrangement_2 *p_arr;
  /*! Boost pointer to the constrained Delaunay triangulation object*/
  boost::shared_ptr<CDT> p_cdt;
  /*! Mapping of the vertices of the input to the corresponding circulator
      needed for finding the first visible vertex in case of face queries*/
  std::map<Point_2, typename Arrangement_2::Ccb_halfedge_const_circulator>
                                                                  point_itr_map;
  const Geometry_traits_2 *geom_traits;
  /*! Stack of visibile points; manipulated when going through the sequence
      of input vertices; contains the vertices of the visibility region after 
      the run of the algorithm*/
  std::stack<Point_2> s;
  /*! Sequence of input vertices*/
  std::vector<Point_2> vertices;
  /*! State of visibility region algorithm*/
  enum {LEFT, RIGHT, SCANA, SCANC, FINISH} upcase;
  bool query_pt_is_vertex;
  bool query_pt_is_on_halfedge;

  /*! Regularize output if flag is set to true*/
  void conditional_regularize(Visibility_arrangement_2& out_arr, CGAL::Tag_true) {
    regularize_output(out_arr);
  }
  /*! No need to regularize output if flag is set to false*/
  void conditional_regularize(Visibility_arrangement_2& out_arr, CGAL::Tag_false) {
    //do nothing
  }

  /*! Regularizes the output - removes edges that have the same face on both
      sides */
  void regularize_output(Visibility_arrangement_2& out_arr) {
    typename Visibility_arrangement_2::Edge_iterator e_itr;
    for (e_itr = out_arr.edges_begin() ; 
         e_itr != out_arr.edges_end() ; e_itr++) {

      typename Visibility_arrangement_2::Halfedge_handle he = e_itr;
      typename Visibility_arrangement_2::Halfedge_handle he_twin = he->twin();
      if (he->face() == he_twin->face()) {
        out_arr.remove_edge(he);
      }
    }
  }

  /*! Initialized the constrained Delaunay triangulation using the edges of
      the outer boundary of 'face' */
  void init_cdt(const Face_const_handle &face) { 

    std::vector<std::pair<Point_2,Point_2> > constraints; 
    typename Arrangement_2::Ccb_halfedge_const_circulator circ = 
                                                            face->outer_ccb();
    typename Arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    typename Arrangement_2::Halfedge_const_handle he;

    do {
      he = curr;
      Point_2 source = he->source()->point();
      Point_2 target = he->target()->point();
      point_itr_map.insert(std::make_pair(source, curr));
      constraints.push_back(std::make_pair(source,target));
    } while(++curr != circ);

    p_cdt = boost::shared_ptr<CDT>(new CDT(constraints.begin(),constraints.end()));
  }


  /*! Finds a visible vertex from the query point 'q' in 'face' 
      to start the algorithm from*/
  Ccb_halfedge_const_circulator find_visible_start(Face_const_handle face, const Point_2 &q) {
    init_cdt(face);
    typename CDT::Face_handle fh = p_cdt->locate(q);
    Point_2 start_point = fh->vertex(0)->point();

    // Now retrieve the circulator to first visible vertex from triangulation
    Ccb_halfedge_const_circulator circ = point_itr_map[start_point];
    Halfedge_const_handle he_curr = circ;

    Halfedge_around_vertex_const_circulator incident_circ = he_curr->source()->incident_halfedges();
    Halfedge_around_vertex_const_circulator incident_curr = incident_circ;

    do {
      Ccb_halfedge_const_circulator curr_inc = incident_curr;
      Halfedge_const_handle he_curr_inc = curr_inc;

      if (he_curr_inc->face() == face) {
        Ccb_halfedge_const_circulator incident_next = incident_curr;
        incident_next++;
        Halfedge_const_handle he_next_inc = incident_next;

        if (CGAL::Visibility_2::orientation_2(geom_traits,
                                          he_curr_inc->source()->point(),
                                          he_curr_inc->target()->point(),
                                          q) == CGAL::LEFT_TURN
         || CGAL::Visibility_2::orientation_2(geom_traits,
                                          he_next_inc->source()->point(),
                                          he_next_inc->target()->point(),
                                          q) == CGAL::LEFT_TURN) {
          Ccb_halfedge_const_circulator result_circ = incident_next;
          Halfedge_const_handle he_print = result_circ;
          return result_circ;
        }
      }
    } while (++incident_curr != incident_circ);
  }

  /*! Main method of the algorithm - initializes the stack and variables
      and calles the corresponding methods acc. to the algorithm's state;
      'q' - query point;
      'i' - current vertex' index
      'w' - endpoint of ray shot from query point */
  void visibility_region_impl(const Point_2& q) {

    int i = 0;
    Point_2 w;
    if (query_pt_is_vertex) {
      if (CGAL::Visibility_2::orientation_2(geom_traits, 
                                            q,
                                            vertices[1],
                                            vertices[2]) == CGAL::LEFT_TURN) {
        upcase = LEFT;
        i = 1;
        w = vertices[1];
        vertices.pop_back();
        s.push(vertices[vertices.size()-1]);
        s.push(vertices[0]);
        s.push(vertices[1]);
      }
      else {
        upcase = SCANA;
        i = 1;
        w = vertices[1];
        vertices.pop_back();
        s.push(vertices[vertices.size()-1]);
        s.push(vertices[0]);
        s.push(vertices[1]);
      }
    }
    else {
      CGAL::Orientation orient = CGAL::Visibility_2::orientation_2(geom_traits, 
                                                                   q, 
                                                                   vertices[0], 
                                                                   vertices[1]);
      if (orient == CGAL::LEFT_TURN || orient == CGAL::COLLINEAR) {

        upcase = LEFT;
        i = 1;
        w = vertices[1];
        s.push(vertices[0]);
        s.push(vertices[1]);
      }
      else {
        upcase = SCANA;
        i = 1;
        w = vertices[0];
        s.push(vertices[0]);
      }
    }
    do {
      switch(upcase) {
        case LEFT: 
          left(i, w, q);
          break;
        case RIGHT:
          right(i, w, q);
          break;
        case SCANA:
          scana(i, w, q);
          break;
        case SCANC:
          scanc(i, w, q);
          break;
      }
    } while(upcase != FINISH);
  }

  /*! Method that handles the left turns in the vertex algorithm */
  void left(int& i, Point_2& w, const Point_2& query_pt) {
    if (i >= vertices.size() - 1) {
      upcase = FINISH;
    }
    else {
        if (w == query_pt) {
          Point_2 s_t = s.top();
          s.pop();
          Point_2 s_t_prev = s.top();
          if (CGAL::Visibility_2::orientation_2(geom_traits,                                                                                                            
                                              s_t_prev, 
                                              vertices[i], 
                                              vertices[i+1]) == CGAL::RIGHT_TURN) {
            upcase = SCANA;
            w = vertices[i];
          }
          else {
            upcase = RIGHT;
            w = vertices[i+1];
            i++;
          }
          s.push(s_t);
        }
        else {
          CGAL::Orientation orient = CGAL::Visibility_2::orientation_2(geom_traits,
                                                                         query_pt, 
                                                                         w, 
                                                                         vertices[i+1]);
          if (orient == CGAL::LEFT_TURN || orient == CGAL::COLLINEAR) {
            s.push(vertices[i+1]);
            upcase = LEFT;
            w = vertices[i+1];
            i++;
          }
          else {
            Point_2 s_t = s.top();
            s.pop();
            Point_2 s_t_prev = s.top();
            if (CGAL::Visibility_2::orientation_2(geom_traits,                                                                                                            
                                                s_t_prev, 
                                                vertices[i], 
                                                vertices[i+1]) == CGAL::RIGHT_TURN) {
              upcase = SCANA;
              w = vertices[i];
            }
            else {
              upcase = RIGHT;
              w = vertices[i+1];
              i++;
            }
            s.push(s_t);
          }
        }
    }
  }

  /*! Scans the stack such that all vertices that were pushed before to the 
      stack and are now not visible anymore. */
  void right(int& i, Point_2& w, const Point_2& query_pt) {

    bool found = false;
    while(!found && upcase == RIGHT) {
      assert(!s.empty());
      Point_2 s_j = s.top();
      s.pop();
      assert(!s.empty());
      Point_2 s_j_prev = s.top();

      if (vertices[i-1] != s_j && CGAL::Visibility_2::do_intersect_2
          <Geometry_traits_2, Segment_2, Segment_2>(geom_traits,
                                                    Segment_2(s_j, s_j_prev), 
                                                    Segment_2(vertices[i-1], 
                                                    vertices[i]))) {
        Segment_2 seg(vertices[i-1], vertices[i]);
        Ray_2 ray(s_j, s_j_prev);
        Object_2 result = CGAL::Visibility_2::intersect_2
                   <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, seg, ray);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          if (*ipoint != s_j_prev) {
            scanb(i, s_j_prev, *ipoint, query_pt, w);
            found = true;
          }
        }
      }
      else {
        CGAL::Orientation orient = 
                                CGAL::Visibility_2::orientation_2(geom_traits, 
                                                                  query_pt,
                                                                  w, 
                                                                  s_j_prev);

        if (orient == CGAL::RIGHT_TURN || orient == CGAL::COLLINEAR) {
          found = true;
          if (i+1 >= vertices.size()) {
            Segment_2 s1(s_j_prev, s_j);
            Ray_2 s2(query_pt, w);
            Object_2 result = CGAL::Visibility_2::intersect_2
                     <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, s1, s2);
            if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
              if (*ipoint != s_j_prev) {
                s.push(*ipoint);
              }
            }
            upcase = FINISH;
            break;
          }
          CGAL::Orientation qwv_orient = 
                              CGAL::Visibility_2::orientation_2(geom_traits, 
                                                                query_pt,
                                                                w, 
                                                                vertices[i+1]);
          if (qwv_orient == CGAL::RIGHT_TURN) {

            upcase = RIGHT;
            w = vertices[i+1];
            s.push(s_j);
            i++;
          }
          else if ((qwv_orient == CGAL::LEFT_TURN 
                 || qwv_orient == CGAL::COLLINEAR) 
              && (CGAL::Visibility_2::orientation_2(geom_traits, 
                                          vertices[i-1],
                                          vertices[i], 
                                          vertices[i+1]) == CGAL::RIGHT_TURN)) {

            Segment_2 s1(s_j_prev, s_j);
            Ray_2 s2(query_pt, w);
            Object_2 result = CGAL::Visibility_2::intersect_2
                     <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, s1, s2);
            if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
              if (i < vertices.size()-1) {
                
                upcase = LEFT;
                if (*ipoint != s_j_prev) {
                  s.push(*ipoint);
                }
                s.push(w);
                s.push(vertices[i+1]); 
                w = vertices[i+1];
                i++;                
              }
              else {
                if (query_pt_is_vertex && *ipoint != s_j_prev) {
                  s.push(*ipoint);
                }
                upcase = FINISH;
              }
            }
          }
          else if (CGAL::Visibility_2::orientation_2(geom_traits,
                                            vertices[i-1],
                                            vertices[i],
                                            vertices[i+1]) == CGAL::LEFT_TURN) {

            upcase = SCANC;
            w = vertices[i];
            s.push(s_j);
            i++;        
          }
          else {
            Segment_2 s1(s_j_prev, s_j);
            Ray_2 s2(query_pt, w);
            Object_2 result = CGAL::Visibility_2::intersect_2
                     <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, s1, s2);
            if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
              upcase = LEFT;
              s.push(*ipoint);
              s.push(w);
              s.push(vertices[i+1]);  
              w = vertices[i+1];
              i++;
            }
          }
        }
      }
    }
  }

  /*! Scans the vertices starting from index 'i' for the first visible vertex
      out of the back hidden window */
  void scana(int& i, Point_2& w, const Point_2& query_pt) {
    // Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
    int k = i;
    while (k+1 < vertices.size()) {
      assert(query_pt != vertices[i+1]);
      CGAL::Orientation qwv_orient = CGAL::Visibility_2::orientation_2(geom_traits,
                                                                       query_pt,
                                                                       w,
                                                                       vertices[k+1]);
      if (qwv_orient == CGAL::LEFT_TURN) {
        Ray_2 s2(query_pt, s.top());
        Segment_2 s1(vertices[k], vertices[k+1]);
        Object_2 result = CGAL::Visibility_2::intersect_2
                       <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, s1, s2);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          s.push(*ipoint);
          s.push(vertices[k+1]);
          w = vertices[k+1];
          i = k+1;
          upcase = LEFT;
          break;
        }
      }
      else if (qwv_orient == CGAL::COLLINEAR) {
        if ((!query_pt_is_vertex && !query_pt_is_on_halfedge) 
            || ((query_pt_is_vertex || query_pt_is_on_halfedge)
          && CGAL::Visibility_2::collinear_are_ordered_along_line_2
                <Geometry_traits_2>(geom_traits, query_pt, w, vertices[k+1]))) {
    
          if (CGAL::Visibility_2::orientation_2(geom_traits, 
                                                query_pt, 
                                                w, 
                                                vertices[k+2]) == RIGHT_TURN) {

            CGAL::Orientation v_orient = 
                          CGAL::Visibility_2::orientation_2(geom_traits, 
                                                            vertices[k], 
                                                            vertices[k+1], 
                                                            vertices[k+2]);
            if (v_orient == LEFT_TURN) {
              upcase = SCANA;
              i = k+2;
            }
            else if (v_orient == RIGHT_TURN) {
              upcase = SCANA;
              s.push(vertices[k+1]);
              w = vertices[k+1];
              i = k+1;
            }
            else {
              upcase = LEFT;
              s.push(vertices[k+1]);
              w = vertices[k+1];
              i = k+1;
            }
          }
          else {
            s.push(vertices[k+1]);
            upcase = LEFT;
            w = vertices[k+1];
            i = k+1;
            break;
          }
          break;
        }
      }
      k++;
    }
  }

  /*! Finds the exit from a front hidden window for the special case from
  right state, when (s[j], s[j-1]) intersects (vertices[i-1], vertices[i])
  by finding the first edge to intersect the segment (st, ipt)*/
  void scanb(int &i, const Point_2 &st, const Point_2 &ipt, 
             const Point_2 &query_pt, Point_2 &w) {
    
    int k = i;
    bool found = false;
    while (!found && k+1 < vertices.size()) {
      Segment_2 seg1(st, ipt);
      Segment_2 seg2(vertices[k], vertices[k+1]);
      if (CGAL::Visibility_2::do_intersect_2
          <Geometry_traits_2, Segment_2, Segment_2>(geom_traits, seg1, seg2)) {

        Object_2 result = CGAL::Visibility_2::intersect_2
             <Geometry_traits_2, Segment_2, Segment_2>(geom_traits, seg1, seg2);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
          if (*ipoint != vertices[k+1]) {
            s.push(*ipoint);
          }
          w = *ipoint;
          i = k;
          upcase = LEFT;
          break;
        }
      }
      k++;
    }
  }

  /*! Finds the exit from a general front hidden window by finding the first
      vertex to the right of the ray defined by the query_point and w*/
  void scanc(int& i, Point_2& w, const Point_2& query_pt) {
    assert(i != vertices.size()-1);
    Point_2 s_t = s.top();
    int k = i;
    Point_2 intersection_pt;
    while (k < vertices.size()) {
      assert(query_pt != vertices[k]);
      CGAL::Orientation qwv_orient = CGAL::Visibility_2::orientation_2(geom_traits, query_pt, w, vertices[k]);
      if (qwv_orient == CGAL::RIGHT_TURN || qwv_orient == CGAL::COLLINEAR) {
        break;
      }
      k++;
    }
    w = vertices[k];
    i = k;
    upcase = RIGHT;
  }
};

} // namespace CGAL
#endif
