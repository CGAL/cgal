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

#include <vector>
#include <CGAL/tags.h>
#include <CGAL/enum.h>

namespace CGAL {
namespace Visibility_2 {

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
  Face_const_iterator fit;
  for (fit = arr.faces_begin() ; fit != arr.faces_end() ; fit++) {
    if (!fit->is_unbounded()) {
      print_simple_face<Face_const_iterator, Ccb_halfedge_const_circulator>(fit);
    }
    else {
      std::cout << "unbounded\n";
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
typename Geometry_traits_2::Object_2 intersect_2(const Geometry_traits_2 *geom_traits,
                     const _Curve_first& s1, 
                     const _Curve_second& s2) {

  typedef typename Geometry_traits_2::Kernel Kernel;
  const Kernel *kernel = static_cast<const Kernel*> (geom_traits);
  typename Kernel::Intersect_2 intersect_fnct = kernel->intersect_2_object();
  return intersect_fnct(s1, s2);
}

template <class Geometry_traits_2>
CGAL::Comparison_result compare_xy_2(const Geometry_traits_2 *geom_traits,
                                    const typename Geometry_traits_2::Point_2 &p,
                                    const typename Geometry_traits_2::Point_2 &q) {

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
bool collinear_are_ordered_along_line_2(const Geometry_traits_2 *geom_traits,
                                        const typename Geometry_traits_2::Point_2 &p,
                                        const typename Geometry_traits_2::Point_2 &q,
                                        const typename Geometry_traits_2::Point_2 &r) {

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

template <class Visibility_2>
void report_while_handling_needles(
  const typename Visibility_2::Input_arrangement_2::Geometry_traits_2 *geom_traits,
  const typename Visibility_2::Input_arrangement_2::Point_2& q,
  std::vector<typename Visibility_2::Input_arrangement_2::Point_2>& points,
  typename Visibility_2::Output_arrangement_2& arr_out) {

  typedef typename Visibility_2::Input_arrangement_2      Input_arrangement_2;
  typedef typename Input_arrangement_2::Point_2           Point_2;
  typedef typename Input_arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Input_arrangement_2::Halfedge_handle   Halfedge_handle;
  typedef typename Input_arrangement_2::Vertex_handle     Vertex_handle;
  typedef typename Geometry_traits_2::Segment_2           Segment_2;
  typedef typename Geometry_traits_2::Direction_2         Direction_2;

  typename std::vector<Segment_2>::size_type i = 0;
  typename std::vector<Segment_2>::size_type start_idx;

  if (points.front() == points.back()) {
    points.pop_back();
  }

  while (collinear(geom_traits,
                   points[i],
                   points[points.size()-1],
                   points[points.size()-2])
      || collinear(geom_traits,
                   points[i],
                   points[i+1],
                   points[points.size()-1])) {

    points.push_back(points[i]);
    i++;
  }

  points.push_back(points[i]);
  start_idx = i;

  Halfedge_handle he_handle;
  Vertex_handle v_trg;
  Vertex_handle v_fst;
  Vertex_handle v_needle_end;
/*
  std::cout << "\nPOINTS\n";
  for(unsigned int k = 0 ; k < points.size() ; k++) {
    std::cout << points[k] << std::endl;
  }
  std::cout << "END\n";
*/
  while (i+1 < points.size()) {
    bool had_needle = false;
    if (i == start_idx) {
      he_handle = arr_out.insert_in_face_interior(Segment_2(points[i], points[i+1]), arr_out.unbounded_face());
 //     std::cout << "inserting in face interior seg " << Segment_2(points[i], points[i+1]) << std::endl;
      v_trg = he_handle->target();
      if (v_trg->point() != points[i+1]) {
        v_fst = he_handle->target();
      }
      else {
        v_fst = he_handle->source();
      }
      i++;
    }
    if (v_trg->point() != points[i]) {
      v_trg = he_handle->source();
    }

    if ((i+2 < points.size()) &&
        (orientation_2(geom_traits, 
                       points[i], 
                       points[i+1], 
                       points[i+2]) == CGAL::COLLINEAR)) {

 //     std::cout << "needle\n";
      had_needle = true;
      std::vector<Point_2> forward_needle;
      std::vector<Point_2> backward_needle;          
      Direction_2 forward_dir(Segment_2(points[i], points[i+1]));
      forward_needle.push_back(points[i]);
      forward_needle.push_back(points[i+1]);

      while ((i+2 < points.size()) && 
            (orientation_2(geom_traits, 
                           points[i], 
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
      Point_2 end_of_needle;
      if (backward_needle.size() != 0) {
        end_of_needle = backward_needle[backward_needle.size()-1];
      }
      else {
        end_of_needle = forward_needle[forward_needle.size()-1];
      }
  //    std::cout << "end of needle set to " << end_of_needle << std::endl;
      std::reverse(backward_needle.begin(), backward_needle.end());
      std::vector<Point_2> merged_needle;

      // Now merge the two vectors
      unsigned int itr_fst = 0, itr_snd = 0;
      while (itr_fst < forward_needle.size() && 
             itr_snd < backward_needle.size()) {

        if (less_distance_to_point_2(geom_traits, 
                                  q, 
                                  forward_needle[itr_fst], 
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
        if (v_trg->point() != merged_needle[p]) {
          v_trg = he_handle->source();
        }
        if (CGAL::Visibility_2::compare_xy_2<Geometry_traits_2>(geom_traits, v_trg->point(), merged_needle[p+1]) == CGAL::SMALLER) {
 //         std::cout << "insertingN seg " << Segment_2(merged_needle[p], merged_needle[p+1]) << " from left of " << v_trg->point() << std::endl;
          he_handle = arr_out.insert_from_left_vertex(Segment_2(merged_needle[p], merged_needle[p+1]), v_trg);  
        }
        else {
 //         std::cout << "insertingN seg " << Segment_2(merged_needle[p], merged_needle[p+1]) << " from right of " << v_trg->point() << std::endl;
          he_handle = arr_out.insert_from_right_vertex(Segment_2(merged_needle[p], merged_needle[p+1]), v_trg);   
        }
        v_trg = he_handle->target();
        if (merged_needle[p+1] == end_of_needle) {
          v_needle_end = v_trg;
  //        std::cout << "set needle end to " << v_needle_end->point();
        }
      }
    }
    else {
      if (had_needle) {
        v_trg = v_needle_end;
  //      std::cout << "set v_trg to " << v_trg->point() << std::endl;
      }
      if (v_trg->point() != points[i]) {
        v_trg = he_handle->source();
      }
      if (CGAL::Visibility_2::compare_xy_2<Geometry_traits_2>(geom_traits, v_trg->point(), points[i+1]) == CGAL::SMALLER) {
  //      std::cout << "inserting seg " << Segment_2(points[i], points[i+1]) << " from left of " << v_trg->point() << std::endl;
        he_handle = arr_out.insert_from_left_vertex(Segment_2(points[i], points[i+1]), v_trg);  
      }
      else {
  //      std::cout << "inserting seg " << Segment_2(points[i], points[i+1]) << " from right of " << v_trg->point() << std::endl;
        he_handle = arr_out.insert_from_right_vertex(Segment_2(points[i], points[i+1]), v_trg);   
      }
      v_trg = he_handle->target();
    }
    i++;

    if (i+2 == points.size()) {
      v_trg = he_handle->target();
      arr_out.insert_at_vertices(Segment_2(points[points.size()-2], points[points.size()-1]), v_trg, v_fst);
      break;
    }
  }
}

template <class Visibility_2>
void report_while_handling_needles_(
  const typename Visibility_2::Input_arrangement_2::Geometry_traits_2 *geom_traits,
  const typename Visibility_2::Input_arrangement_2::Point_2& q,
  std::vector<typename Visibility_2::Input_arrangement_2::Point_2>& points,
  typename Visibility_2::Output_arrangement_2& arr_out) {

  typedef typename Visibility_2::Input_arrangement_2      Input_arrangement_2;
  typedef typename Input_arrangement_2::Point_2           Point_2;
  typedef typename Input_arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Input_arrangement_2::Halfedge_handle   Halfedge_handle;
  typedef typename Input_arrangement_2::Vertex_handle     Vertex_handle;
  typedef typename Geometry_traits_2::Segment_2           Segment_2;
  typedef typename Geometry_traits_2::Direction_2         Direction_2;

  typename std::vector<Segment_2>::size_type i = 0;

  if (points.front() == points.back()) {
    points.pop_back();
  }

  while (collinear(geom_traits,
                   q,
                   points[i],
                   points.back())) {

    points.push_back(points[i]);
    i++;
  }

  points.push_back(points[i]);

  Halfedge_handle he_handle;
  Vertex_handle v_trg;
  Vertex_handle v_fst;
  Vertex_handle v_needle_end;

  v_trg = v_fst = arr_out.insert_in_face_interior(points[i], arr_out.unbounded_face());

//  std::cout << "\nPOINTS\n";
//  for(unsigned int k = 0 ; k < points.size() ; k++) {
//    std::cout << points[k] << std::endl;
//  }
//  std::cout << "END\n";


  while (i+1 < points.size()) {
    if ( collinear(geom_traits,
                   points[i],
                   points[i+1],
                   q)) {
      Vertex_handle v_needle_begin = v_trg;
      std::vector<Point_2> forward_needle;
      std::vector<Point_2> backward_needle;
      std::vector<Point_2> part_in_q_side;
      part_in_q_side.push_back(points[i]);
      forward_needle.push_back((points[i]));
//      bool same_side_of_q = less_distance_to_point_2(geom_traits,
//                                                     q,
//                                                     points[i+1],
//                                                     points[i]);
      bool same_side_of_q = (compare_xy_2(geom_traits, points[i], q)==compare_xy_2(geom_traits, points[i], points[i+1]));
      if (same_side_of_q)
        part_in_q_side.push_back(points[i+1]);
      else
        forward_needle.push_back(points[i+1]);
      i++;
      while (i+1< points.size() && orientation_2(geom_traits,
                           points[i],
                           points[i+1],
                           q ) == CGAL::COLLINEAR) {
        if (same_side_of_q) {
          part_in_q_side.push_back(points[i+1]);
        }
        else {
          if (compare_xy_2(geom_traits, part_in_q_side.front(), q)
              == compare_xy_2(geom_traits, part_in_q_side.front(), points[i+1])) {
            same_side_of_q = true;
            part_in_q_side.push_back(points[i+1]);
          }
          else {
            if (less_distance_to_point_2(geom_traits, q, points[i], points[i+1]))
              forward_needle.push_back(points[i+1]);
            else
              backward_needle.push_back(points[i+1]);
          }
        }
        i++;
      }

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
  //    std::cout << "end of needle set to " << end_of_needle << std::endl;
      std::reverse(backward_needle.begin(), backward_needle.end());
      std::vector<Point_2> merged_needle;

      // Now merge the two vectors
      unsigned int itr_fst = 0, itr_snd = 0;
      while (itr_fst < forward_needle.size() &&
             itr_snd < backward_needle.size()) {

        if (less_distance_to_point_2(geom_traits,
                                  q,
                                  forward_needle[itr_fst],
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
        if (CGAL::Visibility_2::compare_xy_2<Geometry_traits_2>(geom_traits, merged_needle[p], merged_needle[p+1]) == CGAL::SMALLER) {
//          std::cout << "insertingN seg " << Segment_2(merged_needle[p], merged_needle[p+1]) << " from left of " << v_trg->point() << std::endl;
          he_handle = arr_out.insert_from_left_vertex(Segment_2(merged_needle[p], merged_needle[p+1]), v_trg);
        }
        else {
//          std::cout << "insertingN seg " << Segment_2(merged_needle[p], merged_needle[p+1]) << " from right of " << v_trg->point() << std::endl;
          he_handle = arr_out.insert_from_right_vertex(Segment_2(merged_needle[p], merged_needle[p+1]), v_trg);
        }
        v_trg = he_handle->target();
        if (merged_needle[p+1] == end_of_needle) {
          v_needle_end = v_trg;
//         std::cout << "set needle end to " << v_needle_end->point();
        }
      }
      if (same_side_of_q) {
        v_trg = v_needle_begin;
        for (unsigned int p = 0 ; p+1 < part_in_q_side.size() ; p++) {
          if (CGAL::Visibility_2::compare_xy_2<Geometry_traits_2>(geom_traits, part_in_q_side[p], part_in_q_side[p+1]) == CGAL::SMALLER) {
//         std::cout << "insertingN seg " << Segment_2(part_in_q_side[p], part_in_q_side[p+1]) << " from left of " << v_trg->point() << std::endl;
            he_handle = arr_out.insert_from_left_vertex(Segment_2(part_in_q_side[p], part_in_q_side[p+1]), v_trg);
          }
          else {
//        std::cout << "insertingN seg " << Segment_2(part_in_q_side[p], part_in_q_side[p+1]) << " from right of " << v_trg->point() << std::endl;
            he_handle = arr_out.insert_from_right_vertex(Segment_2(part_in_q_side[p], part_in_q_side[p+1]), v_trg);
          }
          v_trg = he_handle->target();
        }
      }
      else
        v_trg = v_needle_end;
    }
    else {
      if (CGAL::Visibility_2::compare_xy_2<Geometry_traits_2>(geom_traits, v_trg->point(), points[i+1]) == CGAL::SMALLER) {
//      std::cout << "inserting seg " << Segment_2(points[i], points[i+1]) << " from left of " << v_trg->point() << std::endl;
        he_handle = arr_out.insert_from_left_vertex(Segment_2(points[i], points[i+1]), v_trg);
      }
      else {
//      std::cout << "inserting seg " << Segment_2(points[i], points[i+1]) << " from right of " << v_trg->point() << std::endl;
        he_handle = arr_out.insert_from_right_vertex(Segment_2(points[i], points[i+1]), v_trg);
      }
      v_trg = he_handle->target();
      i++;
    }
//    std::cout<<"next insertion: "<<v_trg->point()<<std::endl;
    if (i+2 == points.size()) {
      v_trg = he_handle->target();
//      std::cout<<v_fst->point()<<std::endl;
//      std::cout<<v_trg->point()<<std::endl;
      arr_out.insert_at_vertices(Segment_2(points[points.size()-2], points[points.size()-1]), v_trg, v_fst);
      break;
    }
  }
}

} // end namespace Visibility_2
} // end namespace CGAL

#endif
