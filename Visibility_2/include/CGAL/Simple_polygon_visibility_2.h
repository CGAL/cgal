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
#include <stack>

namespace CGAL {

template<class Arrangement_2, class RegularizationTag> 
class Simple_polygon_visibility_2 {

public:
  // Currently only consider with same type for both
  typedef Arrangement_2                                 Input_arrangement_2;
  typedef Arrangement_2                                 Output_arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;

  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
                                                  Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Face_handle     Face_handle;

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
  Simple_polygon_visibility_2(const Input_arrangement_2 &arr): 
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

  Input_arrangement_2 arr() {
    return *p_arr;
  }

  Face_handle visibility_region(Point_2 &q, const Face_handle face,
                         Output_arrangement_2 &out_arr) {

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
    min_intersect_pt = 
                Construct_projected_point_2(curr_min_edge.supporting_line(), q);

    temp_vertices.push_back(curr_vertex);
    Number_type min_dist = Compute_squared_distance_2(q, curr_edge);

    int min_dist_index = 0;
    int index = 1;

    curr++;
    // Push all vertices and determine edge minimum in terms 
    // of squared distance to query point
    do {
      he = curr;          
      curr_edge = Segment_2(he->source()->point(), he->target()->point());
      Number_type curr_dist = Compute_squared_distance_2(q, curr_edge);
        
      if (curr_dist < min_dist) {
        min_dist = curr_dist;
        min_dist_index = index;
        curr_min_edge = curr_edge;
      }
      temp_vertices.push_back(he->target()->point());
      index++;
    } while (++curr != circ);

    // Only now compute the intersection point
    min_intersect_pt = 
                Construct_projected_point_2(curr_min_edge.supporting_line(), q);

    if (min_intersect_pt != curr_min_edge.source() && 
        min_intersect_pt != curr_min_edge.target()) {
      vertices.push_back(min_intersect_pt);
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

    visibility_region_impl(q);

    typename std::vector<Point_2> points;
    if (!s.empty()) {
      Point_2 prev_pt = s.top();
      if (prev_pt == min_intersect_pt) {
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
    std::vector<Segment_2> segments;
    treat_needles(q, points, segments);
    CGAL::insert_non_intersecting_curves(out_arr, 
                                         segments.begin(), 
                                         segments.end());
    CGAL_precondition(out_arr.number_of_isolated_vertices() == 0);
    CGAL_precondition(s.size() == 0);
    conditional_regularize(out_arr, Regularization_tag());
    vertices.clear();

    if (out_arr.faces_begin()->is_unbounded())
      return ++out_arr.faces_begin();
    else
      return out_arr.faces_begin();
  }

  Face_handle visibility_region(const Point_2 &q, const Halfedge_handle he,
                           Output_arrangement_2 &out_arr ) {

    if (q != he->source()->point()) {
      if (q != he->target()->point()) {
        vertices.push_back(q);
        vertices.push_back(he->target()->point());
      }
      else {
        vertices.push_back(q);
      }
    }
    else {
      vertices.push_back(he->target()->point());
    }

    typename Input_arrangement_2::Face_const_handle face = he->face();
    typename Input_arrangement_2::Ccb_halfedge_const_circulator circ = 
                                                              face->outer_ccb();
    typename Input_arrangement_2::Ccb_halfedge_const_circulator curr;
    typename Input_arrangement_2::Halfedge_const_handle he_handle = circ;

    while (he_handle != he) {
      he_handle = circ;
      circ++;
    }

    curr = circ;
    curr++;
    typename Input_arrangement_2::Ccb_halfedge_const_circulator curr_next = curr;
    curr_next++;

    he_handle = curr;
    vertices.push_back(Point_2(he_handle->source()->point()));

    while (curr_next != circ) {
      he_handle = curr;
      Point_2 curr_vertex = he_handle->target()->point();
      vertices.push_back(curr_vertex);
      curr++;
      curr_next++;
    }
    vertices.push_back(vertices[0]);

    visibility_region_impl(q);

    typename std::vector<Point_2> points;
    if (!s.empty()) {
      Point_2 prev_pt = s.top();
      if (prev_pt != q) {
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
        s.pop();
      }
    }

    std::reverse(points.begin(), points.end());
    std::vector<Segment_2> segments;
    treat_needles(q, points, segments);
    CGAL::insert_non_intersecting_curves(out_arr, 
                                         segments.begin(), 
                                         segments.end());
    CGAL_precondition(out_arr.number_of_isolated_vertices() == 0);
    CGAL_precondition(s.size() == 0);
    conditional_regularize(out_arr, Regularization_tag());
    vertices.clear();

    if (out_arr.faces_begin()->is_unbounded())
      return ++out_arr.faces_begin();
    else
      return out_arr.faces_begin();
  }

  void print_arrangement(const Arrangement_2 &arr) {
    typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
    Edge_const_iterator eit;
    std::cout << arr.number_of_edges() << " edges:" << std::endl;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
      std::cout << "[" << eit->curve() << "]" << std::endl;
  }   

private:
  const Input_arrangement_2 *p_arr;
  const Geometry_traits_2  *geom_traits;
  std::stack<Point_2> s;
  std::vector<Point_2> vertices;
  enum {LEFT, RIGHT, SCANA, SCANB, SCANC, SCAND, FINISH} upcase;

  bool LessDistanceToPoint_2(const Point_2 &p, const Point_2 &q, 
                                               const Point_2 &r) const {
    typename Geometry_traits_2::Less_distance_to_point_2 less_dist = 
                                geom_traits->less_distance_to_point_2_object();
    return less_dist(p, q, r);
  }

  bool Collinear(const Point_2 &p, const Point_2 &q,
                                   const Point_2 &r) const {
    typename Geometry_traits_2::Collinear_2 collinear_fnct = 
                                geom_traits->collinear_2_object();
    return collinear_fnct(p, q, r);
  }

  template < class _Curve_first, class _Curve_second >
  Object_2 Intersect_2(const _Curve_first &s1, const _Curve_second &s2) {
    typedef typename Geometry_traits_2::Kernel Kernel;
    const Kernel *kernel = static_cast<const Kernel*> (geom_traits);
    typename Kernel::Intersect_2 intersect_fnct = 
                                            kernel->intersect_2_object();
    return intersect_fnct(s1, s2);
  }

  Orientation Orientation_2(const Point_2 &p, const Point_2 &q, 
                                              const Point_2 &r) {
    typename Geometry_traits_2::Orientation_2 orient = 
                                          geom_traits->orientation_2_object();
    return orient(p, q, r);
  }

  Point_2 Construct_projected_point_2(const Line_2 &l, const Point_2 &p) {
    typename Geometry_traits_2::Construct_projected_point_2 construct_proj =
                            geom_traits->construct_projected_point_2_object();
    return construct_proj(l, p);
  }

  Number_type Compute_squared_distance_2(const Point_2 &p, 
                                         const Segment_2 &seg) {
    typename Geometry_traits_2::Compute_squared_distance_2 compute_dist = 
                              geom_traits->compute_squared_distance_2_object();
    return compute_dist(p, seg);
  }

  bool do_overlap(const Point_2 &a, const Point_2 &b, const Point_2 &c) {
    if (collinear(a, b, c)) {
      Segment_2 s1(a, b);
      Segment_2 s2(a, c);
      const Segment_2 *seg_overlap;
      Object_2 result = Intersect_2<Segment_2, Segment_2>(s1, s2);
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

  void treat_needles(const Point_2 &q, typename std::vector<Point_2> &points, 
                     typename std::vector<Segment_2> &segments) {

    typename std::vector<Point_2>::size_type i = 0;

    while (CGAL::collinear(points[i], points[points.size()-1],
                                     points[points.size()-2]) ||
           CGAL::collinear(points[i], points[i+1], points[points.size()-1])) {

      points.push_back(points[i]);
      i++;
    }

    points.push_back(points[i]);

    std::vector<Point_2> forward_needle;
    std::vector<Point_2> backward_needle;

    while (i+1 < points.size()) {
      if ((i+2 < points.size()) &&
          (Orientation_2(points[i], 
                         points[i+1], 
                         points[i+2]) == CGAL::COLLINEAR)) {
                
        Point_2 needle_start = points[i];
        Direction_2 forward_dir(Segment_2(points[i], points[i+1]));
        forward_needle.push_back(points[i]);
        forward_needle.push_back(points[i+1]);

        while ((i+2 < points.size()) && 
              (Orientation_2(points[i], 
                             points[i+1], 
                             points[i+2]) == CGAL::COLLINEAR)) {

          Direction_2 check_dir(Segment_2(points[i+1], points[i+2]));
          if (forward_dir == check_dir) {
            forward_needle.push_back(points[i+2]);
          }
          else if (check_dir == -forward_dir) {
            backward_needle.push_back(points[i+2]);
          }
          i++;
        }
        std::reverse(backward_needle.begin(), backward_needle.end());

        std::vector<Point_2> merged_needle;

        // Now merge the two vectors
        unsigned int itr_fst = 0, itr_snd = 0;
        while (itr_fst < forward_needle.size() && 
               itr_snd < backward_needle.size()) {

          if (LessDistanceToPoint_2(q, forward_needle[itr_fst], 
                                       backward_needle[itr_snd])) {
              merged_needle.push_back(forward_needle[itr_fst]);
              itr_fst++;
          }
          else {
            merged_needle.push_back(backward_needle[itr_snd]);
            itr_snd++;
          }
        }
        while (itr_fst < forward_needle.size()) {
          merged_needle.push_back(forward_needle[itr_fst]);
          itr_fst++;
        }
        while (itr_snd < backward_needle.size()) {
          merged_needle.push_back(backward_needle[itr_snd]);
          itr_snd++;
        }
        for (unsigned int p = 0 ; p+1 < merged_needle.size() ; p++) {
          segments.push_back(Segment_2(merged_needle[p], merged_needle[p+1]));
        }
      }
      else {
        segments.push_back(Segment_2(points[i], points[i+1]));
      }
      i++;
    }
  }

  void visibility_region_impl(const Point_2 &q) {

    int i = 0;
    Point_2 w;

    if (Orientation_2(q, vertices[0], vertices[1]) == CGAL::LEFT_TURN
        || Orientation_2(q, vertices[0], vertices[1]) == CGAL::COLLINEAR) {
      upcase = LEFT;
      i = 1;
      w = vertices[1];
      s.push(vertices[0]);
      s.push(vertices[1]);
    }
    else {
      upcase = SCANA;
      i = 1;
      w = vertices[1];
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
        case SCANB:
          scanb(i, w, q);
          break;
        case SCANC:
          scanc(i, w, q);
          break;
        case SCAND:
          scand(i, w, q);
          break;
      }

      if (upcase == LEFT) {
        // Check if (s_t-1, s_t) intersects (q, vn) 
        Point_2 s_t = s.top();
        s.pop();
        Point_2 s_t_prev = s.top();
        Segment_2 s1(s_t_prev, s_t);
        Segment_2 s2(q, vertices[vertices.size()-1]);
        Object_2 result = Intersect_2<Segment_2, Segment_2>(s1, s2);

        if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
          Segment_2 s3(s_t_prev, vertices[i]);
          Object_2 result2 = Intersect_2<Segment_2, Segment_2>(s3, s2);
          if (const Point_2 *vertex_new = CGAL::object_cast<Point_2>(&result2)){
            if ((*vertex_new) != (s_t_prev) && (*vertex_new != s_t)) {
              upcase = SCANB;
              s.push(*vertex_new);
            }
            else { // Do not alter stack if it doesn't intersect - push back s_t
              s.push(s_t);
            }
          }
          else {
            s.push(s_t);
          }
        }
        else {
          s.push(s_t);
        }
      }
    } while(upcase != FINISH);
  }

  void left(int &i, Point_2 &w, const Point_2 &query_pt) {
   
    if (i == vertices.size() - 1) {
      upcase = FINISH;
    }
    else if (Orientation_2(query_pt, 
                           vertices[i], 
                           vertices[i+1]) == CGAL::LEFT_TURN
            || Orientation_2(query_pt, 
                           vertices[i], 
                           vertices[i+1]) == CGAL::COLLINEAR) {
      upcase = LEFT;
      s.push(vertices[i+1]);
      w = vertices[i+1];
      i++;
    }
    else if (Orientation_2(query_pt, 
                           vertices[i], 
                           vertices[i+1]) == CGAL::RIGHT_TURN) {
      Point_2 s_t = s.top();
      s.pop();
      Point_2 s_t_prev = s.top();
      s.pop();
      if (Orientation_2(s_t_prev, 
                        vertices[i], 
                        vertices[i+1]) == CGAL::RIGHT_TURN) {
        upcase = SCANA;
        w = vertices[i+1];
        i++;
      } // Both conditions have to be met to move on. Thus same else branch as below
      else {
        upcase = RIGHT;
        w = vertices[i];
        i++;
      }
        s.push(s_t_prev);
        s.push(s_t);
    }
    else {
      upcase = RIGHT;
      i++;
      w = vertices[i];
    }
  }

  void right(int &i, Point_2 &w, const Point_2 &query_pt) {
    // Scan s_t, s_t-1, ..., s_1, s_0 for the first edge (s_j, s_j-1) such that
    // (a) (z, s_j, v_i) is a right turn and (z, s_j-1, v_i) is a left turn, or
    // (b) (z, s_j-1, s_j) is a forward move and (v_i-1, v_i) intersects (s_j-1, s_j)
    bool found = false;
    while(!found && !s.empty()) {
      Point_2 s_j = s.top();
      s.pop();
      if (!s.empty()) {
        Point_2 s_j_prev = s.top();
        // Check condition (a)
        if ((Orientation_2(query_pt, 
                           s_j, 
                           vertices[i]) == CGAL::RIGHT_TURN) &&
            (Orientation_2(query_pt,
                           s_j_prev, 
                           vertices[i]) == CGAL::LEFT_TURN)) {
          found = true;
          Segment_2 s1(s_j_prev, s_j);
          Ray_2 s2(query_pt, vertices[i]);
          Object_2 result = Intersect_2<Segment_2, Ray_2>(s1, s2);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            s_j = *ipoint;
          }

          if (Orientation_2(query_pt,
                            vertices[i], 
                            vertices[i+1]) == CGAL::RIGHT_TURN) {
            upcase = RIGHT;
            s.push(s_j);
            w = vertices[i];
            i++;
          }
          else if ((Orientation_2(query_pt,
                                  vertices[i], 
                                  vertices[i+1]) == CGAL::LEFT_TURN) &&
                   (Orientation_2(vertices[i-1],
                                  vertices[i], 
                                  vertices[i+1]) == CGAL::RIGHT_TURN)) {
            upcase = LEFT;
            s.push(s_j);
            s.push(vertices[i]);
            s.push(vertices[i+1]);
            w = vertices[i+1];
            i++;
          }
          else {
            upcase = SCANC;
            s.push(s_j);
            w = vertices[i];
            i++;        
          }
        }
        else if (do_overlap(query_pt, s_j_prev, s_j)) { // Case (b)
          // Check if v_i-1, v_i intersects (s_j-1, s_j)
          Segment_2 s1(s_j_prev, s_j);
          Segment_2 s2(vertices[i-1], vertices[i]);
          Object_2 result = Intersect_2<Segment_2, Segment_2>(s1, s2);
          if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
            // Keep s_j off the stack
            found = true;
            upcase = SCAND;
            w = *ipoint;
          }
        }
        else if ((Orientation_2(query_pt,
                                s_j, 
                                vertices[i]) == CGAL::RIGHT_TURN) &&
                 (Orientation_2(query_pt,
                                s_j_prev, 
                                vertices[i]) == CGAL::COLLINEAR)) {
          found = true;
          upcase = LEFT;
          s.push(vertices[i]);
          s.push(vertices[i+1]);
          w = vertices[i+1];
          i++;
        }
      }
    }
  }

  void scana(int &i, Point_2 &w, const Point_2 &query_pt) {
    // Scan v_i, v_i+1, ..., v_n for the first edge to intersect (z, s_t)
    bool found = false;
    int k = i;
    Point_2 intersection_pt;
    while (k+1 < vertices.size()) {
      Segment_2 s1(vertices[k], vertices[k+1]);
      Ray_2 s2(query_pt, s.top());
      Object_2 result = Intersect_2<Segment_2, Ray_2>(s1, s2);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
        found = true;
        intersection_pt = *ipoint;
        break;
      }
      k++;
    }
    if (found) {
      if ((Orientation_2(query_pt,
                         vertices[k], 
                         vertices[k+1]) == CGAL::RIGHT_TURN) &&
         (!do_overlap(query_pt, s.top(), intersection_pt))) {
                
        upcase = RIGHT;
        i = k+1;
        w = intersection_pt;
      }
      else if ((Orientation_2(query_pt, 
                              vertices[k], 
                              vertices[k+1]) == CGAL::RIGHT_TURN) &&
               (do_overlap(query_pt, s.top(), intersection_pt))) {

        upcase = SCAND;
        i = k+1;
        w = intersection_pt;
      }
      else if ((Orientation_2(query_pt, 
                              vertices[k], 
                              vertices[k+1]) == CGAL::LEFT_TURN) &&
               (do_overlap(query_pt, s.top(), intersection_pt))) {

        upcase = LEFT;
        i = k+1;
        s.push(intersection_pt);
        if (intersection_pt != vertices[k+1]) {
          s.push(vertices[k+1]);
        }
          w = vertices[k+1];
      }
      else {
          // This case never occurs
      }
    }
  }

  void scanb(int &i, Point_2 &w, const Point_2 &query_pt) {
    // Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, v_n]
    Point_2 s_t = s.top();
    int k = i;
    bool found = false;
    Point_2 intersection_pt;
    while (k+1 < vertices.size()) {
      Segment_2 s1(vertices[k], vertices[k+1]);
      Segment_2 s2(s_t, vertices[vertices.size()-1]);
      Object_2 result = Intersect_2<Segment_2, Segment_2>(s1, s2);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) { 
        if (*ipoint != s_t) {
          intersection_pt = *ipoint;
          found = true;
          break;
        }
      }
      k++;
    }
    if (found) {
      if ((intersection_pt == vertices[k+1]) && 
          (intersection_pt == vertices[vertices.size()-1])) {

        upcase = FINISH;
        w = vertices[vertices.size()-1];
        s.push(vertices[vertices.size()-1]);
      }
      else {
        upcase = RIGHT;
        i = k+1;
        w = intersection_pt;
      }
    }
    else {
      upcase = LEFT;
      i++;
    }
  }

  void scanc(int &i,Point_2 &w, const Point_2 &query_pt) {
    // Scan v_i, v_i+1, ..., v_n-1, v_n for the first edge to intersect (s_t, w)
    Point_2 s_t = s.top();
    int k = i;
    bool found = false;
    Point_2 intersection_pt;
    while (k+1 < vertices.size()) {
      Segment_2 s1(vertices[k], vertices[k+1]);
      Segment_2 s2(s_t, w);
      Object_2 result = Intersect_2<Segment_2, Segment_2>(s1, s2);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        found = true;
        intersection_pt = *ipoint;
        break;
      }
      k++;
    }
    if (found) {
      upcase = RIGHT;
      i = k+1;
      w = intersection_pt;
    }
  }

  void scand(int &i, Point_2 &w, const Point_2 &query_pt) {
    // Scan v_i, v_i+1, v_n-1, v_n for the fist edge to intersect (s_t, w)
    Point_2 s_t = s.top();
    int k = i;
    bool found = false;
    Point_2 intersection_pt;
    while (k+1 < vertices.size()) {
      Segment_2 s1(vertices[k], vertices[k+1]);
      Segment_2 s2(s_t, w);
      Object_2 result = Intersect_2<Segment_2, Segment_2>(s1, s2);
      if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result)) {
        found = true;
        intersection_pt = *ipoint;
        break;
      }
      k++;
    }
    if (found) {
      upcase = LEFT;
      i = k+1;
      s.push(intersection_pt);
      s.push(vertices[k+1]);
      w = vertices[k+1];
    }
  }
};

} // namespace CGAL

#endif
