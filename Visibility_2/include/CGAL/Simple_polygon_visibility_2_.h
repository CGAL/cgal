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

#ifndef CGAL_SIMPLE_POLYGON_VISIBILITY_2__H
#define CGAL_SIMPLE_POLYGON_VISIBILITY_2__H

#include <CGAL/Arrangement_2.h>
#include <CGAL/rational_rotation.h>
#include <CGAL/tags.h>
#include <CGAL/enum.h>
#include <CGAL/Visibility_2/visibility_utils.h>
#include <stack>
#include <map>

namespace CGAL {

template<class Arrangement_2, class RegularizationTag> 
class Simple_polygon_visibility_2_ {

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

  Simple_polygon_visibility_2_() : p_arr(NULL), geom_traits(NULL) {};

  /*! Constructor given an arrangement and the Regularization tag. */
  Simple_polygon_visibility_2_(const Input_arrangement_2 &arr): 
    p_arr(&arr) {
    geom_traits = p_arr->geometry_traits();
  };

  bool is_attached() {
    return (p_arr != NULL);
  }

  void attach(const Input_arrangement_2 &arr) {
    p_arr = &arr;
    geom_traits = p_arr->geometry_traits();
  }

  void detach() {
    p_arr = NULL;
    geom_traits = NULL;
    vertices.clear();
  }

  const Input_arrangement_2& arr() {
    return *p_arr;
  }

  Face_handle visibility_region(Point_2 &q, const Face_const_handle face,
                         Output_arrangement_2 &out_arr) {

    CGAL::Visibility_2::print_arrangement_by_face<Input_arrangement_2>(*p_arr);
   
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
    std::cout << "******VERTICES***************\n";
    for (unsigned int i = 0 ; i < vertices.size() ; i++) {
      std::cout << vertices[i] << std::endl;
    }
    std::cout << "*********************\n";

    compute_angular_displacement(q);
    print_angular_displacement();
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

    std::reverse(points.begin(), points.end());
    std::cout << "POINTS\n";
    for (unsigned int k = 0 ; k < points.size() ; k++) {
      std::cout << points[k] << std::endl;
    }
    std::cout << "END POINTS\n";
    CGAL::Visibility_2::report_while_handling_needles
                              <Simple_polygon_visibility_2_>(geom_traits, 
                                                            q, 
                                                            points,                                 
                                                            out_arr);  
    std::cout << "OUTPUT\n";
    CGAL::Visibility_2::print_arrangement_by_face<Output_arrangement_2>(out_arr);                              
    std::cout << "END OUTPUT\n";
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

  Face_handle visibility_region(const Point_2 &q, const Halfedge_const_handle he,
                           Output_arrangement_2 &out_arr ) {
/*
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
    std::cout << "******VERTICES***************\n";
    for (unsigned int i = 0 ; i < vertices.size() ; i++) {
      std::cout << vertices[i] << std::endl;
    }
    std::cout << "*********************\n";

    visibility_region_impl(q);

    std::cout << "STACK\n";
    while (!s.empty()) {
      std::cout << s.top() << std::endl;
      s.pop();
    }
    std::cout << "END STACK\n";

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

    std::cout << "POINTS\n";
    for (unsigned int i = 0 ; i < points.size() ; i++) {
      std::cout << points[i] << std::endl;
    }
    std::cout << "*****************\n";
    CGAL::Visibility_2::report_while_handling_needles
                              <Simple_polygon_visibility_2>(geom_traits, 
                                                            q, 
                                                            points,                                
                                                            out_arr);
    CGAL_precondition(out_arr.number_of_isolated_vertices() == 0);
    CGAL_precondition(s.size() == 0);
    conditional_regularize(out_arr, Regularization_tag());
    vertices.clear();
//    CGAL::Visibility_2::print_arrangement_by_face<Output_arrangement_2>(out_arr);
    if (out_arr.faces_begin()->is_unbounded()) {
      return ++out_arr.faces_begin();
    }
    else {
      return out_arr.faces_begin();
    }*/
  }

private:
  const Input_arrangement_2 *p_arr;
  const Geometry_traits_2 *geom_traits;
  std::stack<Point_2> s;
  std::vector<Point_2> vertices;
  std::map<Point_2, double> angular_displacement;
  std::pair<Point_2, double> angular_displacement_vn;
  enum {ADVANCE, SCAN, RETARD, FINISH} upcase;
  enum {RAY, SEGMENT} polar_mode;

  bool do_overlap(const Point_2 &a, const Point_2 &b, const Point_2 &c) {
    if (CGAL::Visibility_2::Collinear(geom_traits, a, b, c)) {
      Segment_2 s1(a, b);
      Segment_2 s2(a, c);
      const Segment_2 *seg_overlap;
      Object_2 result = CGAL::Visibility_2::Intersect_2
                 <Geometry_traits_2, Segment_2, Segment_2>(geom_traits, s1, s2);
      if (seg_overlap = CGAL::object_cast<Segment_2>(&result)) { 
        return true;
      }
    }
    return false;
  }

  void conditional_regularize(Output_arrangement_2 &out_arr, CGAL::Tag_true) {
    regularize_output(out_arr);
  }

  void conditional_regularize(Output_arrangement_2 &out_arr, CGAL::Tag_false) {
    //do nothing
  }

  double angle(const Point_2 &r, const Point_2 &p, const Point_2 &q) {

  }

  void compute_angular_displacement(const Point_2 &q) {

    Point_2 v0 = vertices[0];
    angular_displacement.insert(std::pair<Point_2, double>(v0, 0));

    for (unsigned int k = 1 ; k < vertices.size() ; k++) {
      Vector_2 vec1(vertices[k-1].x() - q.x(), vertices[k-1].y() - q.y());
      Vector_2 vec2(vertices[k].x() - q.x(), vertices[k].y() - q.y());
      Number_type scalar_prod = vec1*vec2;
      Number_type len1 = vec1.squared_length();
      Number_type len2 = vec2.squared_length();
      double angle = acos(CGAL::to_double(scalar_prod)/
                                          sqrt(CGAL::to_double(len1*len2)));
      switch(CGAL::Visibility_2::Orientation_2(geom_traits, 
                                               q, 
                                               vertices[k-1], 
                                               vertices[k])) {
        case CGAL::LEFT_TURN:
          if (k == vertices.size() - 1) {
            angular_displacement_vn = std::make_pair(vertices[k],
                          angular_displacement[vertices[k-1]] + angle);
          }
          else {
            angular_displacement.insert(std::pair<Point_2, double>(vertices[k], 
                          angular_displacement[vertices[k-1]] + angle));
          }
          break;
        case CGAL::RIGHT_TURN:
          if (k == vertices.size() - 1) {
            angular_displacement_vn = std::make_pair(vertices[k],
                          angular_displacement[vertices[k-1]] - angle);
          }
          else {
            angular_displacement.insert(std::pair<Point_2, double>(vertices[k], 
                          angular_displacement[vertices[k-1]] - angle));
          }
          break;
        case CGAL::COLLINEAR:
          if (k == vertices.size() - 1) {
            angular_displacement_vn = std::make_pair(vertices[k],
                          angular_displacement[vertices[k-1]]);
          }
          else {
            angular_displacement.insert(std::pair<Point_2, double>(vertices[k], 
                          angular_displacement[vertices[k-1]]));
          }
          break;
      }
    }
  }

  double get_angular_displacement(const Point_2 &pt, const int i) {
    if (i == vertices.size() - 1) {
      return angular_displacement_vn.second;
    }
    else {
      return angular_displacement[pt];
    }
  }

  void print_angular_displacement() {
    typename std::map<Point_2, double>::iterator it = angular_displacement.begin();
    for (it ; it != angular_displacement.end() ; it++) {
      std::cout << it->first << " - " << it->second << std::endl;
    }
    std::cout << angular_displacement_vn.first << " - " 
              << angular_displacement_vn.second << std::endl;
  }

  void regularize_output(Output_arrangement_2 &out_arr) {
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

  void visibility_region_impl(const Point_2 &q) {

    int i = 0;
    Point_2 w;

    bool ccw = false;
    s.push(vertices[0]);
    if (angular_displacement[vertices[1]] >= angular_displacement[vertices[0]]){
      upcase = ADVANCE;
    }
    else {
      upcase = SCAN;
      ccw = true;
      w = Point_2(vertices[0]);
      polar_mode = RAY;
      std::cout << "in here\n";
    }
    while (upcase != FINISH) {
      switch(upcase) {
        case ADVANCE:
          advance(q, i, ccw, w);
          break;
        case RETARD:
          retard(q, i, ccw, w);
          break;
        case SCAN:
          scan(q, i, ccw, w);
          break;
      }
    }
  }

  void advance(const Point_2 &q, int &i, bool &ccw, Point_2 &w) {

    while (upcase == ADVANCE) {

      if (get_angular_displacement(vertices[i+1], i+1) <= 2*CGAL_PI) {
        i++;
        s.push(vertices[i]);

        if (i == vertices.size()-1) {
          upcase = FINISH;
        }
        else if (get_angular_displacement(vertices[i+1], i+1) < get_angular_displacement(vertices[i], i)
                && CGAL::Visibility_2::Orientation_2(geom_traits,
                                                     vertices[i-1], 
                                                     vertices[i],
                                                     vertices[i+1]) == CGAL::RIGHT_TURN) {
          upcase = SCAN;
          ccw = true;
          w = Point_2(vertices[i]);
          polar_mode = RAY;
        }
        else if (get_angular_displacement(vertices[i+1], i+1) < get_angular_displacement(vertices[i], i)
                && CGAL::Visibility_2::Orientation_2(geom_traits,
                                                     vertices[i-1], 
                                                     vertices[i],
                                                     vertices[i+1]) == CGAL::LEFT_TURN) {
          upcase = RETARD;
        }
      }
      else {
        if (angular_displacement[s.top()] < 2*CGAL_PI) {
          // Compute intersection of v[i]v[i+1] and line qv[0]
          Segment_2 seg(vertices[i], vertices[i+1]);
          Ray_2 ray(q, vertices[0]);
          Object_2 result = CGAL::Visibility_2::Intersect_2
                       <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, seg, ray);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            s.push(*ipoint);
            angular_displacement.insert(std::make_pair(*ipoint, 2*CGAL_PI));
          }
        }
        upcase = SCAN;
        ccw = false;
        w = vertices[0];
        polar_mode  = SEGMENT;
      }
    }
  }

  void retard(const Point_2 &q, int &i, bool &ccw, Point_2 &w) {

    while (upcase == RETARD && !s.empty()) {

      Point_2 s_j_prev = s.top();

      if (!s.empty()) {
        s.pop();
        Point_2 s_j = s.top();
       
        if (angular_displacement[s_j] < get_angular_displacement(vertices[i+1], i+1) 
            && get_angular_displacement(vertices[i+1], i+1) > angular_displacement[s_j]) {

          i++;
          // Compute intersection of s[j]s[j+1] and line qv[i]
          Segment_2 seg(s_j, s_j_prev);
          Ray_2 ray(q, vertices[i]);
          Object_2 result = CGAL::Visibility_2::Intersect_2
                       <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, seg, ray);

          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {

            angular_displacement.insert(std::pair<Point_2, double>(*ipoint, CGAL::to_double(2*CGAL_PI)));
            s.push(*ipoint);
          }
          s.push(vertices[i]);
          if (i == vertices.size() - 1) {
            upcase = FINISH;
          }
          else if (get_angular_displacement(vertices[i], i) <= get_angular_displacement(vertices[i+1], i+1)
                  && CGAL::Visibility_2::Orientation_2(geom_traits,
                                                       vertices[i-1], 
                                                       vertices[i],
                                                       vertices[i+1]) == CGAL::RIGHT_TURN) {

            upcase = ADVANCE;
          }
          else if (get_angular_displacement(vertices[i], i) < get_angular_displacement(vertices[i+1], i+1)
                && CGAL::Visibility_2::Orientation_2(geom_traits,
                                                     vertices[i-1], 
                                                     vertices[i],
                                                     vertices[i+1]) == CGAL::LEFT_TURN) {
 
            upcase = SCAN;
            ccw = false;
            w = vertices[i];
            polar_mode = SEGMENT;
            s.pop();
          }
          else {
            s.pop();
          }
        }
        else {

          if (get_angular_displacement(vertices[i+1], i+1) == angular_displacement[s_j]
                && get_angular_displacement(vertices[i+2], i+2) > get_angular_displacement(vertices[i+1], i+1)
                && CGAL::Visibility_2::Orientation_2(geom_traits,
                                                     vertices[i], 
                                                     vertices[i+1],
                                                     vertices[i+2]) == CGAL::RIGHT_TURN) {
            upcase = ADVANCE;
            i++;
            s.push(vertices[i]);
          }
          else {
            Segment_2 seg_fst(vertices[i], vertices[i+1]);
            Segment_2 seg_snd(s_j, s_j_prev);
            Object_2 result = CGAL::Visibility_2::Intersect_2
                         <Geometry_traits_2, Segment_2, Segment_2>(geom_traits, seg_fst, seg_snd);

            if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
              w = *ipoint;
              polar_mode = SEGMENT;
              upcase = SCAN;
              ccw = true;
            }
          }
        }
      }
    }
  }

  void scan(const Point_2 &q, int &i, bool &ccw, Point_2 &w) {

    while (upcase == SCAN && i < vertices.size()-2) {
      i++;

      if (ccw && get_angular_displacement(vertices[i+1], i+1) > angular_displacement[s.top()]
              && angular_displacement[s.top()] >= get_angular_displacement(vertices[i], i)) {

        Segment_2 seg(vertices[i], vertices[i+1]);
        if (polar_mode == RAY) {
          Ray_2 w_ray(q, w);
          Object_2 result = CGAL::Visibility_2::Intersect_2
                         <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, seg, w_ray);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            upcase = ADVANCE;
            s.push(*ipoint);
          }
        }
        else {
          Segment_2 w_seg(s.top(), w);
          Object_2 result = CGAL::Visibility_2::Intersect_2
                         <Geometry_traits_2, Segment_2, Segment_2>(geom_traits, seg, w_seg);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            upcase = ADVANCE;
            s.push(*ipoint);
          }
        }
      }
      else if (!ccw && get_angular_displacement(vertices[i+1], i+1) <= angular_displacement[s.top()]
                    && angular_displacement[s.top()] < get_angular_displacement(vertices[i], i)) {

        Segment_2 seg(vertices[i], vertices[i+1]);
        if (polar_mode == RAY) {
          Ray_2 w_ray(q, w);
          Object_2 result = CGAL::Visibility_2::Intersect_2
                         <Geometry_traits_2, Segment_2, Ray_2>(geom_traits, seg, w_ray);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            upcase = RETARD;
          }
        }
        else {
          Segment_2 w_seg(s.top(), w);
          Object_2 result = CGAL::Visibility_2::Intersect_2
                         <Geometry_traits_2, Segment_2, Segment_2>(geom_traits, seg, w_seg);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            upcase = RETARD;
          }
        }
      }
    }
  }
};
} // namespace CGAL

#endif
