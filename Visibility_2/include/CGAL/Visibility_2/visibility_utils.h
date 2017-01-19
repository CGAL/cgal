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

#ifndef CGAL_VISIBILITY_UTILS_H
#define CGAL_VISIBILITY_UTILS_H

#include <CGAL/license/Visibility_2.h>


#include <iostream>
#include <vector>
#include <CGAL/tags.h>
#include <CGAL/enum.h>

namespace CGAL {
namespace Visibility_2 {

  template <class Arrangement_2>
int count_edges_in_face(typename Arrangement_2::Face_const_handle fch) {
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
                                                  Ccb_halfedge_const_circulator;

  Ccb_halfedge_const_circulator circ = fch->outer_ccb();
  Ccb_halfedge_const_circulator curr = circ;

  int edge_cnt(0);
  do {
    edge_cnt++;
  } while (++curr != circ);
  return edge_cnt;
}

template <class Edge_const_iterator>
void print_edge(Edge_const_iterator eit) {
  std::cout << "[" << eit->curve() << "]" << std::endl;
}
template <class Face_const_handle, class Ccb_halfedge_const_circulator>
void print_simple_face(Face_const_handle fh) {
  Ccb_halfedge_const_circulator  cir = fh->outer_ccb();
  Ccb_halfedge_const_circulator  curr = cir;
  do {
    std::cout << "[" << curr->curve() << "]" << std::endl;
  } while (++ curr != cir);
}

template <class Arrangement_2> 
void print_arrangement(const Arrangement_2& arr) {
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  Edge_const_iterator eit;
  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
    print_edge(eit);
}

template <class Arrangement_2>
void print_arrangement_by_face(const Arrangement_2& arr) {
  typedef typename Arrangement_2::Face_const_iterator     Face_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator 
                                                 Ccb_halfedge_const_circulator;
  Face_const_iterator f;
  for (f = arr.faces_begin() ; f != arr.faces_end() ; f++) {
    if (!f->is_unbounded()) {
      std::cout << "FACE\n";
      print_simple_face<Face_const_iterator, Ccb_halfedge_const_circulator>(f);
      std::cout << "END FACE\n";
    }
  }
}

template <class Geometry_traits_2> 
Orientation orientation_2(const Geometry_traits_2 *geom_traits,
                          const typename Geometry_traits_2::Point_2& p, 
                          const typename Geometry_traits_2::Point_2& q, 
                          const typename Geometry_traits_2::Point_2& r) {

  typename Geometry_traits_2::Orientation_2 orient = 
                                          geom_traits->orientation_2_object();
  return orient(p, q, r);
}

template <class Geometry_traits_2>
bool less_distance_to_point_2(const Geometry_traits_2 *geom_traits,
                           const typename Geometry_traits_2::Point_2& p, 
                           const typename Geometry_traits_2::Point_2& q, 
                           const typename Geometry_traits_2::Point_2& r) {

  typename Geometry_traits_2::Less_distance_to_point_2 less_dist = 
                                geom_traits->less_distance_to_point_2_object();
  return less_dist(p, q, r);
}

template <class Geometry_traits_2>
bool collinear(const Geometry_traits_2 *geom_traits,
               const typename Geometry_traits_2::Point_2& p,
               const typename Geometry_traits_2::Point_2& q,
               const typename Geometry_traits_2::Point_2& r) {

    typename Geometry_traits_2::Collinear_2 collinear_fnct = 
                                geom_traits->collinear_2_object();
    return collinear_fnct(p, q, r);
}

template <class Geometry_traits_2, class _Curve_first, class _Curve_second >
typename Geometry_traits_2::Object_2 intersect_2(
        const Geometry_traits_2 *geom_traits,
        const _Curve_first& s1,
        const _Curve_second& s2)
{

  typedef typename Geometry_traits_2::Kernel Kernel;
  const Kernel *kernel = static_cast<const Kernel*> (geom_traits);
  typename Kernel::Intersect_2 intersect_fnct = kernel->intersect_2_object();
  return intersect_fnct(s1, s2);
}

template <class Geometry_traits_2>
CGAL::Comparison_result compare_xy_2(
        const Geometry_traits_2 *geom_traits,
        const typename Geometry_traits_2::Point_2 &p,
        const typename Geometry_traits_2::Point_2 &q)
{

  typename Geometry_traits_2::Compare_xy_2 cmp = 
                            geom_traits->compare_xy_2_object();
  return cmp(p, q);
}

template <class Geometry_traits_2, class Type1, class Type2>
typename Geometry_traits_2::FT compute_squared_distance_2(  
                            const Geometry_traits_2 *geom_traits,
                            const Type1& p, 
                            const Type2& seg) {

  typename Geometry_traits_2::Compute_squared_distance_2 compute_dist = 
                              geom_traits->compute_squared_distance_2_object();
  return compute_dist(p, seg);
}

template <class Geometry_traits_2, class Type1, class Type2>
bool do_intersect_2(const Geometry_traits_2 *geom_traits,
                    const Type1& c1, 
                    const Type2& c2) {

  typename Geometry_traits_2::Do_intersect_2 intersect = 
                              geom_traits->do_intersect_2_object();
  return intersect(c1, c2);
}

template <class Geometry_traits_2> 
bool collinear_are_ordered_along_line_2(
        const Geometry_traits_2 *geom_traits,
        const typename Geometry_traits_2::Point_2 &p,
        const typename Geometry_traits_2::Point_2 &q,
        const typename Geometry_traits_2::Point_2 &r)
{

  typename Geometry_traits_2::Collinear_are_ordered_along_line_2 coll = 
                      geom_traits->collinear_are_ordered_along_line_2_object();
  return coll(p, q, r);
}

template <class Geometry_traits_2, class Type1, class Type2>
typename Geometry_traits_2::Point_2 construct_projected_point_2(
                                const Geometry_traits_2 *geom_traits,
                                const Type1& s, 
                                const Type2& p) {

  typedef typename Geometry_traits_2::Point_2         Point_2;
  typedef typename Geometry_traits_2::FT              Number_type;
  typename Geometry_traits_2::Construct_projected_point_2 construct_proj =
                              geom_traits->construct_projected_point_2_object();
  Point_2 proj_pt = construct_proj(s.supporting_line(), p);
  if (s.has_on(proj_pt)) {
    return proj_pt;
  }
  else {
    Number_type d_to_src = compute_squared_distance_2
        <Geometry_traits_2, Point_2, Point_2>(geom_traits, proj_pt, s.source());
    Number_type d_to_trg = compute_squared_distance_2
        <Geometry_traits_2, Point_2, Point_2>(geom_traits, proj_pt, s.target());              
    if (d_to_src < d_to_trg) {
      return s.source();
    }
    else {
      return s.target();
    }
  }
}

// construct an arrangement of visibility region from a vector of
// circular ordered vertices with respect to the query point
template <class Visibility_2, class Visibility_arrangement_2>
void report_while_handling_needles(
  const typename Visibility_2::Arrangement_2::Geometry_traits_2 *traits,
  const typename Visibility_2::Arrangement_2::Point_2& q,
  std::vector<typename Visibility_2::Arrangement_2::Point_2>& points,
  Visibility_arrangement_2& arr_out) {

  typedef typename Visibility_2::Arrangement_2              Arrangement_2;
  typedef typename Arrangement_2::Point_2                   Point_2;
  typedef typename Arrangement_2::Geometry_traits_2         Geometry_traits_2;
  typedef typename Geometry_traits_2::Segment_2             Segment_2;

  typedef typename Visibility_arrangement_2::Halfedge_handle    Halfedge_handle;
  typedef typename Visibility_arrangement_2::Vertex_handle      Vertex_handle;



  typename std::vector<Segment_2>::size_type i = 0;

  if (points.front() == points.back()) {
    points.pop_back();
  }

  while (collinear(traits, q, points[i], points.back())) {
    points.push_back(points[i]);
    i++;
  }

  points.push_back(points[i]);

  Halfedge_handle he_handle;

  //the handle of vertex where the next segment is inserted
  Vertex_handle v_trg;

  //the handle of vertex inserted first
  Vertex_handle v_fst;

  //the handle of vertex of the end of a needle
  Vertex_handle v_needle_end;


  v_trg = v_fst = arr_out.insert_in_face_interior(points[i],
                                                  arr_out.unbounded_face());

  //find a point that is right after a needle
  while (i+1 < points.size()) {
    if ( collinear(traits, points[i], points[i+1], q)) {
      Vertex_handle v_needle_begin = v_trg;

      std::vector<Point_2> forward_needle; //vertices of the needle that are not
                                           //between q and v_needle_begin;
                                           //their direction is leaving q;

      std::vector<Point_2> backward_needle;//vertices of the needle that are not
                                           //between q and v_needle_begin;
                                           //their direction is towards q;

      std::vector<Point_2> part_in_q_side; //vertices of the needle that are
                                           //between q and v_needle_begin
      part_in_q_side.push_back(points[i]);
      forward_needle.push_back((points[i]));

      bool same_side_of_q = (compare_xy_2(traits, points[i], q) ==
                             compare_xy_2(traits, points[i], points[i+1]));

      if (same_side_of_q)
        part_in_q_side.push_back(points[i+1]);
      else
        forward_needle.push_back(points[i+1]);
      i++;
      while (i+1 < points.size() &&
             orientation_2(traits, points[i], points[i+1], q ) ==
             CGAL::COLLINEAR)
      {
        if (same_side_of_q) {
          part_in_q_side.push_back(points[i+1]);
        }
        else {
          if (compare_xy_2(traits, part_in_q_side.front(), q) ==
              compare_xy_2(traits, part_in_q_side.front(), points[i+1]))
          {
            same_side_of_q = true;
            part_in_q_side.push_back(points[i+1]);
          }
          else {
            if (less_distance_to_point_2(traits, q, points[i], points[i+1]))
              forward_needle.push_back(points[i+1]);
            else
              backward_needle.push_back(points[i+1]);
          }
        }
        i++;
      }
      //obtain the end point of a needle
      Point_2 end_of_needle;
      if (same_side_of_q)
        end_of_needle = part_in_q_side.back();
      else {
        if (backward_needle.empty()) {
          end_of_needle = forward_needle.back();
        }
        else {
          end_of_needle = backward_needle.back();
        }
      }

      std::reverse(backward_needle.begin(), backward_needle.end());
      std::vector<Point_2> merged_needle;

      // merge the forward_needle and backward_needle
      unsigned int itr_fst = 0, itr_snd = 0;
      while (itr_fst < forward_needle.size() &&
             itr_snd < backward_needle.size()) {

        if (less_distance_to_point_2(traits,
                                     q,
                                     forward_needle[itr_fst],
                                     backward_needle[itr_snd]))
        {
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
        if (compare_xy_2<Geometry_traits_2>(
                    traits, merged_needle[p], merged_needle[p+1]) == SMALLER)
        {
          he_handle = arr_out.insert_from_left_vertex(
                      Segment_2(merged_needle[p], merged_needle[p+1]), v_trg);
        }
        else {
          he_handle = arr_out.insert_from_right_vertex(
                      Segment_2(merged_needle[p], merged_needle[p+1]), v_trg);
        }
        v_trg = he_handle->target();
        if (merged_needle[p+1] == end_of_needle) {
          v_needle_end = v_trg;
        }
      }
      if (same_side_of_q) {
        //insert the part of needle between q and v_needle_begin
        v_trg = v_needle_begin;
        for (unsigned int p = 0 ; p+1 < part_in_q_side.size() ; p++) {
          if (compare_xy_2<Geometry_traits_2>(
                     traits, part_in_q_side[p], part_in_q_side[p+1]) == SMALLER)
          {
            he_handle = arr_out.insert_from_left_vertex(
                      Segment_2(part_in_q_side[p], part_in_q_side[p+1]), v_trg);
          }
          else {
            he_handle = arr_out.insert_from_right_vertex(
                      Segment_2(part_in_q_side[p], part_in_q_side[p+1]), v_trg);
          }
          v_trg = he_handle->target();
        }
      }
      else
        v_trg = v_needle_end;
    }
    else {
      if (compare_xy_2<Geometry_traits_2>(
                  traits, v_trg->point(), points[i+1]) == SMALLER)
      {
        he_handle = arr_out.insert_from_left_vertex(
                    Segment_2(points[i], points[i+1]), v_trg);
      }
      else {
        he_handle = arr_out.insert_from_right_vertex(
                    Segment_2(points[i], points[i+1]), v_trg);
      }
      v_trg = he_handle->target();
      i++;
    }
    if (i+2 == points.size()) {
      //close the boundary
      v_trg = he_handle->target();
      arr_out.insert_at_vertices(Segment_2(points[points.size()-2],
                                           points[points.size()-1]),
                                 v_trg, v_fst);
      break;
    }
  }
}

template <typename VARR>
void regularize_output(VARR& arr_out) {
  typename VARR::Edge_iterator it = arr_out.edges_begin();

  while(it != arr_out.edges_end()) {
    if (it->face() == it->twin()->face()) {
      typename VARR::Halfedge_handle he = it;
      ++it;
      arr_out.remove_edge(he);
    }
    else {
      ++it;
    }
  }
}

template <typename VARR>
void conditional_regularize(VARR& arr_out, CGAL::Tag_true) {
  regularize_output(arr_out);
}

template <typename VARR>
void conditional_regularize(VARR&, CGAL::Tag_false) {
  //do nothing
}



} // end namespace Visibility_2
} // end namespace CGAL

#endif
