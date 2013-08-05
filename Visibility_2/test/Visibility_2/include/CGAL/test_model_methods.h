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

#ifndef CGAL_TEST_MODEL_METHODS_H
#define CGAL_TEST_MODEL_METHODS_H

#include <CGAL/basic.h>
#include <CGAL/test_utils.h>
#include <cassert>
#include <vector>

namespace CGAL {

template <class Visibility_2>
void test_model_methods_for_arr(
              const typename Visibility_2::Input_arrangement_2 &arr) {

  typedef typename Visibility_2::Input_arrangement_2       Input_arrangement_2;
  typedef typename Visibility_2::Output_arrangement_2      Output_arrangement_2;
  typedef typename Input_arrangement_2::Point_2            Point_2;
  typedef typename Input_arrangement_2::Geometry_traits_2::Segment_2 
                                                           Segment_2;

  Visibility_2 visibility;
  assert(false == visibility.is_attached());
  visibility.attach(arr);
  assert(true == visibility.is_attached());
  assert(true == (CGAL::test_are_equal<Input_arrangement_2>(arr, 
                                                            visibility.arr())));

  Output_arrangement_2 arr_out;
  Output_arrangement_2 arr_out_check;

  typename Input_arrangement_2::Face_const_iterator fit;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (!fit->is_unbounded()) {
      break;
    }
  }
  // First consider query point in the unbounded face
  Point_2 query_pt(1, 1);
  visibility.visibility_region(query_pt, fit, arr_out);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr, arr_out));
  visibility.detach();
  assert(false == visibility.is_attached());
  visibility.attach(arr);

  visibility.visibility_region(query_pt, fit, arr_out_check);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_out_check, arr));
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_out, arr_out_check));
  arr_out.clear();
  visibility.visibility_region(query_pt, fit, arr_out);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_out, arr_out_check));

  // Now consider the query point on a halfedge
  query_pt = Point_2(0, 4);
  arr_out.clear();
  typename Input_arrangement_2::Halfedge_const_iterator hit;

  for (hit = arr.halfedges_begin(); 
       hit != arr.halfedges_end(); ++hit) {

    Segment_2 curr_seg(hit->source()->point(), hit->target()->point());
    if (curr_seg.has_on(query_pt)) {
      visibility.visibility_region(query_pt, hit, arr_out);
      break;
    }
  }
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_out, arr));
  arr_out_check.clear();
  visibility.visibility_region(query_pt, hit, arr_out_check);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_out, arr_out_check));
}

template <class Visibility_2>
void test_model_methods() {

  // Check concept obediance
	typedef typename Visibility_2::Input_arrangement_2       Input_arrangement_2;
	typedef typename Visibility_2::Output_arrangement_2      Output_arrangement_2;
  typedef typename Input_arrangement_2::Point_2            Point_2;
  typedef typename Input_arrangement_2::Face_handle        Face_handle;
  typedef typename Input_arrangement_2::Halfedge_handle    Halfedge_handle;
  typedef typename Visibility_2::Regularization_tag        Regularization_tag;
  typedef typename Visibility_2::Supports_general_polygon_tag
                                                  Supports_general_polygon_tag;
  typedef typename Visibility_2::Supports_simple_polygon_tag 
                                                  Supports_simple_polygon_tag;
  typedef typename Input_arrangement_2::Geometry_traits_2::Segment_2 
                                                           Segment_2;
                                         
  Point_2 p1(0, 0), p2(8, 0), p3(8, 8), p4(0, 8);
  std::vector<Segment_2> seg_sq;
  seg_sq.push_back(Segment_2(p1, p2));
  seg_sq.push_back(Segment_2(p2, p3));
  seg_sq.push_back(Segment_2(p3, p4));
  seg_sq.push_back(Segment_2(p4, p1));
  Input_arrangement_2 arr_square;
  CGAL::insert(arr_square, seg_sq.begin(), seg_sq.end());

  test_model_methods_for_arr<Visibility_2>(arr_square);

  std::vector<Segment_2> seg_tri;
  seg_tri.push_back(Segment_2(p1, p2));
  seg_tri.push_back(Segment_2(p2, p4));
  seg_tri.push_back(Segment_2(p4, p1));
  Input_arrangement_2 arr_triangle;
  CGAL::insert(arr_triangle, seg_tri.begin(), seg_tri.end());

  test_model_methods_for_arr<Visibility_2>(arr_triangle);
}

} // end CGAL namespace
#endif