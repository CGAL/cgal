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

#ifndef CGAL_TEST_SIMPLE_POLYGONS_H
#define CGAL_TEST_SIMPLE_POLYGONS_H

namespace CGAL {

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_face_test_case(std::ifstream &input, std::ifstream &correct_output) {

    typedef _Visibility_2                                   Visibility_2;
    typedef _Arrangement_2                                  Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2       Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2             Point_2;

    Visibility_2 visibility;

    // First read arrangement 
    Arrangement_2 arr, out_arr, correct_out_arr;

    CGAL::create_arrangement_from_file<Arrangement_2>(arr, input);

    // Read query point from file
    double x, y;
    input >> x >> y;
    Point_2 query_pt(x, y);

    CGAL::create_arrangement_from_file<Arrangement_2>(correct_out_arr, correct_output);

    typename Arrangement_2::Face_handle fit;

    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (!fit->is_unbounded()) {
            visibility.visibility_region(query_pt, fit, out_arr);
        }
    }

    return CGAL::test_are_equal<Arrangement_2>(correct_out_arr, out_arr);    
}

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_halfedge_test_case(std::ifstream &input, std::ifstream &correct_output) {

    typedef _Visibility_2                                   Visibility_2;
    typedef _Arrangement_2                                  Arrangement_2;
    typedef typename Arrangement_2::Geometry_traits_2       Geometry_traits_2;
    typedef typename Geometry_traits_2::Point_2             Point_2;
    typedef typename Geometry_traits_2::Segment_2           Segment_2;


    // First read arrangement 
    Arrangement_2 arr, out_arr, correct_out_arr;

    CGAL::create_arrangement_from_file<Arrangement_2>(arr, input);
    Visibility_2 visibility(arr);

    // Read query point from file
    double x, y;
    input >> x >> y;
    Point_2 query_pt(x, y);

    CGAL::create_arrangement_from_file<Arrangement_2>(correct_out_arr, correct_output);

    typename Arrangement_2::Halfedge_iterator hit;

    for (hit = arr.halfedges_begin(); hit != arr.halfedges_end(); ++hit) {
        Segment_2 curr_seg(hit->source()->point(), hit->target()->point());
        if (curr_seg.has_on(query_pt)) {
            visibility.visibility_region(query_pt, hit, out_arr);
            break;
        }
    }
    return CGAL::test_are_equal<Arrangement_2>(correct_out_arr, out_arr);    
}

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_test_case_1() {

    std::ifstream input("./data/simple_polygon_test_case_1.in");
    std::ifstream correct_output("./data/simple_polygon_test_case_1.out");
    return simple_polygon_face_test_case<_Visibility_2, _Arrangement_2>(input, correct_output);
}

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_test_case_2() {

    std::ifstream input("./data/simple_polygon_test_case_2.in");
    std::ifstream correct_output("./data/simple_polygon_test_case_2.out");
    return simple_polygon_face_test_case<_Visibility_2, _Arrangement_2>(input, correct_output);
}

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_test_case_3() {
    std::ifstream input("./data/simple_polygon_test_case_3.in");
    std::ifstream correct_output("./data/simple_polygon_test_case_3.out");
    return simple_polygon_face_test_case<_Visibility_2, _Arrangement_2>(input, correct_output);
}

template < class _Visibility_2, class _Arrangement_2 >
bool simple_polygon_test_case_4() {
    std::ifstream input("./data/simple_polygon_test_case_4.in");
    std::ifstream correct_output("./data/simple_polygon_test_case_4.out");
    return simple_polygon_halfedge_test_case<_Visibility_2, _Arrangement_2>(input, correct_output);
}



} // end namespace CGAL

#endif
