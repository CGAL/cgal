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
// $URL$
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
#include <stack>

namespace CGAL {

template<class Arrangement_2, class RegularizationTag> 
class Simple_polygon_visibility_2 {

public:
  // Currently only consider with same type for both
  typedef Arrangement_2                                 Input_arrangement_2;
  typedef Arrangement_2                                 Output_arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;

  typedef typename Arrangement_2::Halfedge_const_handle       
                                                        Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

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
  Simple_polygon_visibility_2(const Input_arrangement_2& arr): 
    p_arr(&arr) {
    geom_traits = p_arr->geometry_traits();
    query_pt_is_vertex = false;
  }

  bool is_attached() {
    return (p_arr != NULL);
  }

  void attach(const Input_arrangement_2& arr) {
    p_arr = &arr;
    geom_traits = p_arr->geometry_traits();
    query_pt_is_vertex = false;
  }

  void detach() {
    p_arr = NULL;
    geom_traits = NULL;
    vertices.clear();
    query_pt_is_vertex = false;
  }

  const Input_arrangement_2& arr() {
    return *p_arr;
  }

  Face_handle compute_visibility(const Point_2& q, const Face_const_handle face,
                         Output_arrangement_2& out_arr) {

    assert(query_pt_is_vertex == false);
   
    typename Input_arrangement_2::Ccb_halfedge_const_circulator circ = 
                                                            face->outer_ccb();
    typename Input_arrangement_2::Ccb_halfedge_const_circulator curr = circ;
    typename Input_arrangement_2::Halfedge_const_handle he = curr;

    std::vector<Point_2> temp_vertices;
    Point_2 min_intersect_pt;
    bool intersect_on_endpoint = false;
    Segment_2 curr_edge(he->source()->point(), he->target()->point());
    Segment_2 curr_min_edge(he->source()->point(), he->target()->point());
    Point_2 curr_vertex = he->target()->point();
    temp_vertices.push_back(curr_vertex);
    Number_type min_dist = CGAL::Visibility_2::Compute_squared_distance_2
             <Geometry_traits_2, Point_2, Segment_2>(geom_traits, q, curr_edge);
    int min_dist_index = 0;
    int index = 1;
    curr++;
    // Push all vertices and determine edge minimum in terms 
    // of squared distance to query point
    do {
      he = curr;          
      curr_edge = Segment_2(he->source()->point(), he->target()->point());
      Number_type curr_dist = CGAL::Visibility_2::Compute_squared_distance_2
             <Geometry_traits_2, Point_2, Segment_2>(geom_traits, q, curr_edge);
        
      if (curr_dist < min_dist) {
        min_dist = curr_dist;
        min_dist_index = index;
        curr_min_edge = curr_edge;
      }
      temp_vertices.push_back(he->target()->point());
      index++;
    } while (++curr != circ);

    // Only now compute the intersection point
    min_intersect_pt = CGAL::Visibility_2::Construct_projected_point_2
         <Geometry_traits_2, Segment_2, Point_2>(geom_traits, curr_min_edge, q);

    bool intersect_pt_on_seg_endpoint = false;
    if (min_intersect_pt != curr_min_edge.source() && 
        min_intersect_pt != curr_min_edge.target()) {
      vertices.push_back(min_intersect_pt);
    }
    else {
      intersect_pt_on_seg_endpoint = true;
    }
    // Now create vector so that first vertex v0 is visible
    for (unsigned int k = min_dist_index ; k < temp_vertices.size() ; k++) {
      vertices.push_back(temp_vertices[k]);
    }
    for (unsigned int k = 0 ; k < min_dist_index ; k++) {
      vertices.push_back(temp_vertices[k]);
    }

    // Push first vertex again to fulfill algo precondition
    if (min_intersect_pt != curr_min_edge.source() && 
        min_intersect_pt != curr_min_edge.target()) {
      vertices.push_back(min_intersect_pt);
    }
    else {
      vertices.push_back(vertices[0]);
    }
    std::cout << "VERTICES\n";
    for(unsigned int i = 0 ; i < vertices.size() ; i++) {
      std::cout << vertices[i] << std::endl;
    }
    std::cout << "END VERTICES\n";

    visibility_region_impl(q);

    typename std::vector<Point_2> points;
    if (!s.empty()) {
      Point_2 prev_pt = s.top();
      if (prev_pt == min_intersect_pt) {
        if (intersect_pt_on_seg_endpoint) {
          points.push_back(prev_pt);
        }
        s.pop();
        if (!s.empty()) {
          prev_pt = s.top();
          points.push_back(prev_pt);
        }
      }
      if (!s.empty()) {
        s.pop();
      }
      while(!s.empty()) {
        Point_2 curr_pt = s.top();
        if (curr_pt == min_intersect_pt) {
          if (intersect_pt_on_seg_endpoint) {
            points.push_back(curr_pt);
          }
          s.pop();
        }
        else {
          points.push_back(curr_pt);
          prev_pt = curr_pt;
          s.pop();
        }
      }
    }

    std::cout << "POINTS\n";
    for (unsigned int i = 0 ; i < points.size() ; i++) {
      std::cout << points[i]<<std::endl;
    }
    std::cout << "END POINTS\n";

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

  Face_handle compute_visibility(const Point_2& q, const Halfedge_const_handle he,
                           Output_arrangement_2& out_arr ) {

    query_pt_is_vertex = false;

    if (q != he->source()->point()) {
      if (q != he->target()->point()) {
        vertices.push_back(q);
        vertices.push_back(he->target()->point());
      }
      else {
        vertices.push_back(q);
        query_pt_is_vertex = true;
      }
    }

    typename Input_arrangement_2::Face_const_handle face = he->face();
    typename Input_arrangement_2::Ccb_halfedge_const_circulator circ = 
                                                              face->outer_ccb();
    typename Input_arrangement_2::Ccb_halfedge_const_circulator curr;
    typename Input_arrangement_2::Halfedge_const_handle he_handle = circ;

    while (he_handle != he) {
      circ++;
      he_handle = circ;
    }
    circ++;
    curr = circ;
    do {
      he_handle = curr;
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
  const Input_arrangement_2 *p_arr;
  const Geometry_traits_2 *geom_traits;
  std::stack<Point_2> s;
  std::vector<Point_2> vertices;
  enum {LEFT, RIGHT, SCANA, SCANC, FINISH} upcase;
  bool query_pt_is_vertex;

  void conditional_regularize(Output_arrangement_2& out_arr, CGAL::Tag_true) {
    regularize_output(out_arr);
  }

  void conditional_regularize(Output_arrangement_2& out_arr, CGAL::Tag_false) {
    //do nothing
  }

  void regularize_output(Output_arrangement_2& out_arr) {
    typename Output_arrangement_2::Edge_iterator e_itr;
    for (e_itr = out_arr.edges_begin() ; 
         e_itr != out_arr.edges_end() ; e_itr++) {

      Halfedge_handle he = e_itr;
      Halfedge_handle he_twin = he->twin();
      if (he->face() == he_twin->face()) {
        out_arr.remove_edge(he);
      }
    }
  }

  void visibility_region_impl(const Point_2& q) {

    int i = 0;
    Point_2 w;
    if (query_pt_is_vertex) {
      upcase = LEFT;
      i = 1;
      w = vertices[1];
      vertices.pop_back();
      s.push(vertices[vertices.size()-1]);
      s.push(vertices[0]);
      s.push(vertices[1]);
    }
    else if (CGAL::Visibility_2::Orientation_2(geom_traits, 
                                          q, 
                                          vertices[0], 
                                          vertices[1]) == CGAL::LEFT_TURN
      || CGAL::Visibility_2::Orientation_2(geom_traits, 
                                          q, 
                                          vertices[0], 
                                          vertices[1]) == CGAL::COLLINEAR) {

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

  void left(int& i, Point_2& w, const Point_2& query_pt) {
    std::cout << "in left\n";
    if (i == vertices.size() - 1) {
      upcase = FINISH;
    }
    else if (CGAL::Visibility_2::Orientation_2(geom_traits,
                                               query_pt, 
                                               w, 
                                               vertices[i+1]) == CGAL::LEFT_TURN
            || CGAL::Visibility_2::Orientation_2(geom_traits,
                                            query_pt, 
                                            w, 
                                            vertices[i+1]) == CGAL::COLLINEAR) {

      upcase = LEFT;
      s.push(vertices[i+1]);
      w = vertices[i+1];
      i++;
    }
    else {

      Point_2 s_t = s.top();
      s.pop();
      Point_2 s_t_prev = s.top();                                                                                            
      if (CGAL::Visibility_2::Orientation_2(geom_traits,                                                                                                            
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

  void right(int& i, Point_2& w, const Point_2& query_pt) {
    // Scan s_t, s_t-1, ..., s_1, s_0 for the first edge (s_j, s_j-1) such that
    // (z, s_j, v_i) is a right turn and (z, s_j-1, v_i) is a left turn, or
    bool found = false;
    std::cout << "right w = " << w << std::endl;
    std::cout << "R: v[i+1] = " << vertices[i+1] << std::endl;
    std::cout << "R: v[i-1] = " << vertices[i-1] << std::endl;
    std::cout << "R: v[i] = " << vertices[i] << std::endl;
    std::cout << "R: s top = " << s.top() << std::endl;
    while(!found && upcase == RIGHT) {
      assert(!s.empty());
      Point_2 s_j = s.top();
      std::cout << "R s.top = " << s.top();
      s.pop();
      assert(!s.empty());
      Point_2 s_j_prev = s.top();

      if (vertices[i-1] != s_j && CGAL::do_intersect(Segment_2(s_j, s_j_prev), Segment_2(vertices[i-1], vertices[i]))) {
        std::cout << "R switch to scana\n";
        upcase = SCANA;
        found = true;
        w = s.top();
      }
      else {

        std::cout << "R: s t-1 = " << s_j_prev << std::endl;
        assert(CGAL::Visibility_2::Orientation_2(geom_traits, 
                                                 query_pt, 
                                                 w, 
                                                 s_j) == CGAL::LEFT_TURN 
            || CGAL::Visibility_2::Orientation_2(geom_traits, 
                                                 query_pt, 
                                                 w, 
                                                 s_j) == CGAL::COLLINEAR);

        if (CGAL::Visibility_2::Orientation_2(geom_traits, 
                                              query_pt,
                                              w, 
                                              s_j_prev) == CGAL::RIGHT_TURN
          || CGAL::Visibility_2::Orientation_2(geom_traits, 
                                              query_pt,
                                              w, 
                                              s_j_prev) == CGAL::COLLINEAR) {
          found = true;
          std::cout << "R found\n";

          if (CGAL::Visibility_2::Orientation_2(geom_traits, 
                                                query_pt,
                                                w, 
                                                vertices[i+1]) == CGAL::RIGHT_TURN) {
          // SEE TEST CASE 
            upcase = RIGHT;
            w = vertices[i+1];
            s.push(s_j);
            i++;
            std::cout << "R continues with w = " << w << std::endl;
          }
          else if ((CGAL::Visibility_2::Orientation_2(geom_traits, 
                                                      query_pt,
                                                      w, 
                                                      vertices[i+1]) == CGAL::LEFT_TURN
                  || CGAL::Visibility_2::Orientation_2(geom_traits, 
                                                      query_pt,
                                                      w, 
                                                      vertices[i+1]) == CGAL::COLLINEAR) &&
                    (CGAL::Visibility_2::Orientation_2(geom_traits, 
                                                       vertices[i-1],
                                                       vertices[i], 
                                                       vertices[i+1]) == CGAL::RIGHT_TURN)) {
            std::cout << "in jere\n";
            Segment_2 s1(s_j_prev, s_j);
            Ray_2 s2(query_pt, w);
            Object_2 result = CGAL::Visibility_2::Intersect_2
                        <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, s1, s2);
            if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
              std::cout << "ipt = " << *ipoint << std::endl;
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
                std::cout << "bla\n";
                upcase = FINISH;
              }
            }
          }
          else if (CGAL::Visibility_2::Orientation_2(geom_traits,
                                                     vertices[i-1],
                                                     vertices[i],
                                                     vertices[i+1]) == CGAL::LEFT_TURN) {

            std::cout << "R: scanc\n";
            upcase = SCANC;
            w = vertices[i];
      //      s.push(s_j);
            i++;        
          }
          else {
            std::cout << "R : left\n";
            Segment_2 s1(s_j_prev, s_j);
            Ray_2 s2(query_pt, w);
            Object_2 result = CGAL::Visibility_2::Intersect_2
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

  void scana(int& i, Point_2& w, const Point_2& query_pt) {
    // Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
    int k = i;
    std::cout << "scana w = " << w << std::endl;
    while (k+1 < vertices.size()-1) {

      if (CGAL::Visibility_2::Orientation_2(geom_traits,
                                            query_pt,
                                            w,
                                            vertices[k+1]) == CGAL::LEFT_TURN) {
        Segment_2 s1(vertices[k], vertices[k+1]);
        Ray_2 s2(query_pt, s.top());

        Object_2 result = CGAL::Visibility_2::Intersect_2
                       <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, s1, s2);
        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
          s.push(*ipoint);
          std::cout << "scana ipt = " << *ipoint << std::endl;
          s.push(vertices[k+1]);
          w = vertices[k+1];
          i = k+1;
          std::cout << "scana i = " << i << std::endl;
          upcase = LEFT;
          break;
        }
      }
      else if (CGAL::Visibility_2::Orientation_2(geom_traits,
                                            query_pt,
                                            w,
                                            vertices[k+1]) == CGAL::COLLINEAR) {

        s.push(vertices[k+1]);
        w = vertices[k+1];
        i = k+1;
        std::cout << "scana i = " << i << std::endl;
        upcase = LEFT;
        break;
      }
      k++;
    }
  }

  void scanc(int& i, Point_2& w, const Point_2& query_pt) {
    // Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, w)
    assert(i != vertices.size()-1);
    Point_2 s_t = s.top();
    int k = i;
    Point_2 intersection_pt;
    while (k < vertices.size()-1) {
      if (CGAL::Visibility_2::Orientation_2(geom_traits, query_pt, w, vertices[k]) == CGAL::RIGHT_TURN
        ||CGAL::Visibility_2::Orientation_2(geom_traits, query_pt, w, vertices[k]) == CGAL::COLLINEAR) {
        break;
      }
      k++;
    }
    w = vertices[k];
    std::cout << "scanc w = " << w << std::endl;
    i = k;
    upcase = RIGHT;
  }
};

} // namespace CGAL

#endif
