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
void print_arrangement(const Arrangement_2 &arr) {
  typedef typename Arrangement_2::Edge_const_iterator Edge_const_iterator;
  Edge_const_iterator eit;
  std::cout << arr.number_of_edges() << " edges:" << std::endl;
  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
    print_edge(eit);
}   

template <class Geometry_traits_2> 
Orientation Orientation_2(const Geometry_traits_2 *geom_traits,
                          const typename Geometry_traits_2::Point_2 &p, 
                          const typename Geometry_traits_2::Point_2 &q, 
                          const typename Geometry_traits_2::Point_2 &r) {

  typename Geometry_traits_2::Orientation_2 orient = 
                                          geom_traits->orientation_2_object();
  return orient(p, q, r);
}

template <class Geometry_traits_2>
bool LessDistanceToPoint_2(const Geometry_traits_2 *geom_traits,
                           const typename Geometry_traits_2::Point_2 &p, 
                           const typename Geometry_traits_2::Point_2 &q, 
                           const typename Geometry_traits_2::Point_2 &r) {

  typename Geometry_traits_2::Less_distance_to_point_2 less_dist = 
                                geom_traits->less_distance_to_point_2_object();
  return less_dist(p, q, r);
}

template <class Geometry_traits_2>
bool Collinear(const Geometry_traits_2 *geom_traits,
               const typename Geometry_traits_2::Point_2 &p,
               const typename Geometry_traits_2::Point_2 &q,
               const typename Geometry_traits_2::Point_2 &r) {

    typename Geometry_traits_2::Collinear_2 collinear_fnct = 
                                geom_traits->collinear_2_object();
    return collinear_fnct(p, q, r);
}

template <class Geometry_traits_2, class _Curve_first, class _Curve_second >
typename Geometry_traits_2::Object_2 Intersect_2(const Geometry_traits_2 *geom_traits,
                     const _Curve_first &s1, 
                     const _Curve_second &s2) {

  typedef typename Geometry_traits_2::Kernel Kernel;
  const Kernel *kernel = static_cast<const Kernel*> (geom_traits);
  typename Kernel::Intersect_2 intersect_fnct = kernel->intersect_2_object();
  return intersect_fnct(s1, s2);
}

template <class Geometry_traits_2>
typename Geometry_traits_2::Point_2 Construct_projected_point_2(
                                const Geometry_traits_2 *geom_traits,
                                const typename Geometry_traits_2::Line_2 &l, 
                                const typename Geometry_traits_2::Point_2 &p) {

  typename Geometry_traits_2::Construct_projected_point_2 construct_proj =
                              geom_traits->construct_projected_point_2_object();
  return construct_proj(l, p);
}

template <class Geometry_traits_2>
typename Geometry_traits_2::FT Compute_squared_distance_2(  
                            const Geometry_traits_2 *geom_traits,
                            const typename Geometry_traits_2::Point_2 &p, 
                            const typename Geometry_traits_2::Segment_2 &seg) {

  typename Geometry_traits_2::Compute_squared_distance_2 compute_dist = 
                              geom_traits->compute_squared_distance_2_object();
  return compute_dist(p, seg);
}

template <class Visibility_2>
void report_while_handling_needles(
  const typename Visibility_2::Input_arrangement_2::Geometry_traits_2 *geom_traits,
  const typename Visibility_2::Input_arrangement_2::Point_2 &q,
  std::vector<typename Visibility_2::Input_arrangement_2::Point_2> &points,
  typename Visibility_2::Output_arrangement_2 &arr_out) {

  typedef typename Visibility_2::Input_arrangement_2      Input_arrangement_2;
  typedef typename Input_arrangement_2::Point_2           Point_2;
  typedef typename Input_arrangement_2::Geometry_traits_2 Geometry_traits_2;
  typedef typename Geometry_traits_2::Segment_2           Segment_2;
  typedef typename Geometry_traits_2::Direction_2         Direction_2;

  typename std::vector<Segment_2>::size_type i = 0;

  if (points[0] == points[points.size()-1]) {
    points.pop_back();
  }

  while (Collinear(geom_traits, 
                   points[i], 
                   points[points.size()-1],
                   points[points.size()-2]) 
      || Collinear(geom_traits, 
                   points[i], 
                   points[i+1], 
                   points[points.size()-1])) {

    points.push_back(points[i]);
    i++;
  }

  points.push_back(points[i]);

  std::vector<Point_2> forward_needle;
  std::vector<Point_2> backward_needle;
  std::vector<Segment_2> segments;

  while (i+1 < points.size()) {
    if ((i+2 < points.size()) &&
        (Orientation_2(geom_traits, 
                       points[i], 
                       points[i+1], 
                       points[i+2]) == CGAL::COLLINEAR)) {
                
      Point_2 needle_start = points[i];
      Direction_2 forward_dir(Segment_2(points[i], points[i+1]));
      forward_needle.push_back(points[i]);
      forward_needle.push_back(points[i+1]);

      while ((i+2 < points.size()) && 
            (Orientation_2(geom_traits, 
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
      std::reverse(backward_needle.begin(), backward_needle.end());
      std::vector<Point_2> merged_needle;

      // Now merge the two vectors
      unsigned int itr_fst = 0, itr_snd = 0;
      while (itr_fst < forward_needle.size() && 
             itr_snd < backward_needle.size()) {

        if (LessDistanceToPoint_2(geom_traits, 
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
        segments.push_back(Segment_2(merged_needle[p], merged_needle[p+1]));
      }
    }
    else {
      segments.push_back(Segment_2(points[i], points[i+1]));
    }
    i++;
  }
 
  CGAL::insert_non_intersecting_curves(arr_out, 
                                       segments.begin(), 
                                       segments.end());
}

} // end namespace Visibility_2
} // end namespace CGAL

#endif
