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

template <class Visibility_2, class Visibility_arrangement_2>
void test_model_methods_for_arr(
              typename Visibility_2::Arrangement_2 &arr) {

  typedef typename Visibility_2::Arrangement_2              Arrangement_2;
  typedef typename Arrangement_2::Point_2                   Point_2;


  typedef typename Visibility_arrangement_2::Face_handle    VFH;


  Visibility_2 visibility;
  assert(false == visibility.is_attached());
  visibility.attach(arr);
  assert(true == visibility.is_attached());
  assert(true == (CGAL::test_are_equal(arr,visibility.arrangement_2())));

  Visibility_arrangement_2 arr_out;
  Visibility_arrangement_2 arr_out_check;

  typename Arrangement_2::Face_const_iterator fit;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if (!fit->is_unbounded()) {
      break;
    }
  }
  // First consider query point in the unbounded face
  const Point_2 query_pt(1, 1);
  // Check returned face_handle
  VFH face_check = visibility.compute_visibility(query_pt, fit, arr_out);
  VFH face;
  if (arr_out.faces_begin()->is_unbounded()) {
    face = ++arr_out.faces_begin();
  }
  else {
    face = arr_out.faces_begin();
  }
  assert(face_check == face);

  assert((true == test_are_equal<Arrangement_2, Visibility_arrangement_2>
          (arr, arr_out)));
  visibility.detach();
  assert(false == visibility.is_attached());
  visibility.attach(arr);

  visibility.compute_visibility(query_pt, fit, arr_out_check);
  assert((true == test_are_equal<Visibility_arrangement_2, Arrangement_2>
          (arr_out_check, arr)));
  assert((true == test_are_equal<Visibility_arrangement_2, Visibility_arrangement_2>
          (arr_out, arr_out_check)));
  arr_out.clear();
  visibility.compute_visibility(query_pt, fit, arr_out);
  assert((true == test_are_equal<Visibility_arrangement_2, Visibility_arrangement_2>
          (arr_out, arr_out_check)));

  // Now consider the query point on a halfedge
  const Point_2 query_pt2 = Point_2(0, 4);
  arr_out.clear();
  typename Arrangement_2::Halfedge_const_iterator hit;
  VFH face_check_he;
  for (hit = arr.halfedges_begin(); 
    hit != arr.halfedges_end(); ++hit) {

    if (hit->source()->point() == Point_2(0, 8) && hit->target()->point() == Point_2(0, 0)) {
      face_check_he = visibility.compute_visibility(query_pt2, hit, arr_out);
      break;
    }
  }
  if (arr_out.faces_begin()->is_unbounded()) {
    face = ++arr_out.faces_begin();
  }
  else {
    face = arr_out.faces_begin();
  }
  assert(face_check_he == face);
  assert((true == test_are_equal(arr_out, arr)));
  arr_out_check.clear();
  visibility.compute_visibility(query_pt2, hit, arr_out_check);
  assert((true == test_are_equal<Visibility_arrangement_2, Visibility_arrangement_2>
          (arr_out, arr_out_check)));  

  // Now consider the query point as the target of a halfedge
  typename Arrangement_2::Halfedge_const_iterator hit_snd;
  for (hit_snd = arr.halfedges_begin(); hit_snd != arr.halfedges_end(); ++hit_snd) {
    if(!hit_snd->face()->is_unbounded()){
      arr_out.clear();
      VFH face_check_he_snd = visibility.compute_visibility(hit_snd->target()->point(), hit_snd, arr_out);
      assert(!face_check_he_snd->is_unbounded());
      if (arr_out.faces_begin()->is_unbounded()) {
        face = ++arr_out.faces_begin();
      }
      else {
        face = arr_out.faces_begin();
      }
      assert(face_check_he_snd == face);
      if (! test_are_equal(arr_out, arr)) {
        assert(false);
      }
    }
  }   
}

template <class Visibility_2, class Visibility_arrangement_2>
void test_model_methods() {

  typedef typename Visibility_2::Arrangement_2                  Arrangement_2;
  typedef typename Arrangement_2::Point_2                       Point_2;
  typedef typename Arrangement_2::Geometry_traits_2::Segment_2  Segment_2;


  Point_2 p1(0, 0), p2(8, 0), p3(8, 8), p4(0, 8);
  std::vector<Segment_2> seg_sq;
  seg_sq.push_back(Segment_2(p1, p2));
  seg_sq.push_back(Segment_2(p2, p3));
  seg_sq.push_back(Segment_2(p3, p4));
  seg_sq.push_back(Segment_2(p4, p1));
  Arrangement_2 arr_square;
  CGAL::insert(arr_square, seg_sq.begin(), seg_sq.end());

  test_model_methods_for_arr<Visibility_2,Visibility_arrangement_2>(arr_square);

  std::vector<Segment_2> seg_tri;
  seg_tri.push_back(Segment_2(p1, p2));
  seg_tri.push_back(Segment_2(p2, p4));
  seg_tri.push_back(Segment_2(p4, p1));
  Arrangement_2 arr_triangle;
  CGAL::insert(arr_triangle, seg_tri.begin(), seg_tri.end());

  test_model_methods_for_arr<Visibility_2,Visibility_arrangement_2>(arr_triangle);
}

} // end CGAL namespace
#endif
