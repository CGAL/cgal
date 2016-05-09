// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <fstream>
#include <iostream>

#include <CGAL/Timer.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include "AABB_test_util.h"

enum Query_type {RAY_QUERY,
                 SEGMENT_QUERY,
                 LINE_QUERY};

template <class Tree, class K>
void test_speed_for_query(const Tree& tree,
                          const Query_type query_type,
                          const char *query_name,
                          const double duration)
{
    typedef typename K::Ray_3 Ray;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::Segment_3 Segment;

    CGAL::Timer timer;
    unsigned int nb = 0;
    timer.start();
    while(timer.time() < duration)
    {
        switch(query_type)
        {
            case RAY_QUERY:
                {
                    Point source = random_point_in<K>(tree.bbox());
                    Vector vec = random_vector<K>();
                    Ray ray(source, vec);
                    tree.do_intersect(ray);
                    break;
                }
            case SEGMENT_QUERY:
                {
                    Point a = random_point_in<K>(tree.bbox());
                    Point b = random_point_in<K>(tree.bbox());
                    tree.do_intersect(Segment(a,b));
                    break;
                }
                break;
            case LINE_QUERY:
                {
                    Point a = random_point_in<K>(tree.bbox());
                    Point b = random_point_in<K>(tree.bbox());
                    tree.do_intersect(Line(a,b));
                    break;
                }
        }
        nb++;
    }
    unsigned int speed = (unsigned int)(nb / timer.time());
    std::cout.precision(10);
    std::cout.width(15);
    std::cout << speed << " intersections/s with " << query_name << std::endl;
    timer.stop();
}

template <class Tree, class K>
void test_speed(Tree& tree,
                const double duration)
{
    std::cout << "Test for speed" << std::endl;
    test_speed_for_query<Tree,K>(tree,RAY_QUERY,"ray",duration);
    test_speed_for_query<Tree,K>(tree,LINE_QUERY,"line",duration);
    test_speed_for_query<Tree,K>(tree,SEGMENT_QUERY,"segment",duration);
}

template<class K, class Tree, class Polyhedron, Primitive_type Type>
void test_impl(Tree& tree, Polyhedron&, const double duration)
{
  test_all_intersection_query_types<Tree,K>(tree);
  test_speed<Tree,K>(tree,duration);
}

int main()
{
    const double duration = 0.1; // duration of each test
    test_kernels<TRIANGLE>("./data/cube.off",duration);
    test_kernels<TRIANGLE>("./data/coverrear.off",duration);
    test_kernels<TRIANGLE>("./data/finger.off",duration);
    test_kernels<TRIANGLE>("./data/pinion.off",duration);
    return EXIT_SUCCESS;
}
