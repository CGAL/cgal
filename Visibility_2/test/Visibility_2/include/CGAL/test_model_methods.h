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
              
  // Test with simple square and query point at its center                                    
  Point_2 p1(0, 0), p2(8, 0), p3(8, 8), p4(0, 8);
  Point_2 query_pt(4, 4);
  std::vector<Segment_2> seg_sq;
  seg_sq.push_back(Segment_2(p1, p2));
  seg_sq.push_back(Segment_2(p2, p3));
  seg_sq.push_back(Segment_2(p3, p4));
  seg_sq.push_back(Segment_2(p4, p1));
  Input_arrangement_2 arr_square;
  CGAL::insert(arr_square, seg_sq.begin(), seg_sq.end());

  Visibility_2 visibility;
  assert(false == visibility.is_attached());
  visibility.attach(arr_square);
  assert(true == visibility.is_attached());
  assert(true == (CGAL::test_are_equal<Input_arrangement_2>(arr_square, 
                                                            visibility.arr())));

  Output_arrangement_2 arr_square_out;
  typename Input_arrangement_2::Face_const_iterator fit;

  for (fit = arr_square.faces_begin(); fit != arr_square.faces_end(); ++fit) {
    if (!fit->is_unbounded()) {
      break;
    }
  }
  visibility.visibility_region(query_pt, fit, arr_square_out);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_square, arr_square_out));
  visibility.detach();
  assert(false == visibility.is_attached());
  visibility.attach(arr_square);
  Output_arrangement_2 arr_square_out_check;

  visibility.visibility_region(query_pt, fit, arr_square_out_check);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_square_out_check, arr_square));
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_square_out, arr_square_out_check));
  arr_square_out.clear();
  visibility.visibility_region(query_pt, fit, arr_square_out);
  assert(true == test_are_equal<Output_arrangement_2>
                                        (arr_square_out, arr_square_out_check));

  std::vector<Segment_2> seg_tri;
  seg_tri.push_back(Segment_2(p1, p2));
  seg_tri.push_back(Segment_2(p2, p4));
  seg_tri.push_back(Segment_2(p4, p1));
  Input_arrangement_2 arr_triangle;
  CGAL::insert(arr_triangle, seg_tri.begin(), seg_tri.end());
}

} // end CGAL namespace
#endif