// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#include <CGAL/AABB_intersections.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include "AABB_test_util.h"

template <class Tree, class K>
void test_all_query_types(Tree& tree)
{
    std::cout << "Test all query types" << std::endl;

    typedef typename K::FT FT;
    typedef typename K::Ray_3 Ray;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::Segment_3 Segment;
    typedef typename Tree::Primitive Primitive;
    typedef typename Tree::Point_and_primitive Point_and_primitive;
    typedef typename Tree::Object_and_primitive Object_and_primitive;

    Point p((FT)-0.5, (FT)-0.5, (FT)-0.5);
    Point q((FT) 0.5, (FT) 0.5, (FT) 0.5);
    Ray ray(p,q);
    Ray line(p,q);
    Ray segment(p,q);
    bool success = false;

    // do_intersect
    success = tree.do_intersect(ray);
    success = tree.do_intersect(line);
    success = tree.do_intersect(segment);

    // number_of_intersected_primitives
    tree.number_of_intersected_primitives(ray);
    tree.number_of_intersected_primitives(line);
    tree.number_of_intersected_primitives(segment);

    // all_intersected_primitives
    std::list<Primitive> primitives;
    tree.all_intersected_primitives(ray,std::back_inserter(primitives));
    tree.all_intersected_primitives(line,std::back_inserter(primitives));
    tree.all_intersected_primitives(segment,std::back_inserter(primitives));

    // any_intersection
    boost::optional<Object_and_primitive> optional_object_and_primitive;
    optional_object_and_primitive = tree.any_intersection(ray);
    optional_object_and_primitive = tree.any_intersection(line);
    optional_object_and_primitive = tree.any_intersection(segment);

    // any_intersected_primitive
    boost::optional<Primitive> optional_primitive;
    optional_primitive = tree.any_intersected_primitive(ray);
    optional_primitive = tree.any_intersected_primitive(line);
    optional_primitive = tree.any_intersected_primitive(segment);

    // all_intersections
    std::list<Object_and_primitive> intersections;
    tree.all_intersections(ray,std::back_inserter(intersections));
    tree.all_intersections(line,std::back_inserter(intersections));
    tree.all_intersections(segment,std::back_inserter(intersections));
}


enum Query_type {RAY_QUERY,
                 SEGMENT_QUERY,
                 LINE_QUERY};

template <class Tree, class K>
void test_speed_for_query(const Tree& tree,
                          const Query_type query_type,
                          const char *query_name)
{
    typedef typename K::FT FT;
    typedef typename K::Ray_3 Ray;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::Segment_3 Segment;

    CGAL::Timer timer;
    unsigned int nb = 0;
    timer.start();
    while(timer.time() < 1.0)
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
void test_speed(Tree& tree)
{
    std::cout << "Test for speed" << std::endl;
    test_speed_for_query<Tree,K>(tree,RAY_QUERY,"ray");
    test_speed_for_query<Tree,K>(tree,LINE_QUERY,"line");
    test_speed_for_query<Tree,K>(tree,SEGMENT_QUERY,"segment");
}

template <class K>
void test(const char *filename)
{
    typedef typename K::FT FT;
    typedef typename K::Ray_3 Ray;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    // loads triangle polyhedral surface
    Polyhedron polyhedron;
    std::ifstream ifs(filename);
    ifs >> polyhedron;

    // constructs tree 
    std::cout << "construct tree...";
    CGAL::Timer timer;
    timer.start();
    Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());
    timer.stop();
    std::cout << "done (" << timer.time() << " s)" << std::endl;

    // tests rebuilding
    tree.rebuild(polyhedron.facets_begin(),polyhedron.facets_end());

    // calls all tests
    test_all_query_types<Tree,K>(tree);
    test_speed<Tree,K>(tree);
}

void test_several_kernels(const char *filename)
{
    std::cout << std::endl;
    std::cout << "Polyhedron " << filename << std::endl;
    std::cout << "============================" << std::endl;

    std::cout << std::endl;
    std::cout << "Simple cartesian float kernel" << std::endl;
    test<CGAL::Simple_cartesian<float> >(filename);

    std::cout << std::endl;
    std::cout << "Cartesian float kernel" << std::endl;
    test<CGAL::Cartesian<float> >(filename);

    std::cout << std::endl;
    std::cout << "Simple cartesian double kernel" << std::endl;
    test<CGAL::Simple_cartesian<double> >(filename);

    std::cout << std::endl;
    std::cout << "Cartesian double kernel" << std::endl;
    test<CGAL::Cartesian<double> >(filename);

    std::cout << std::endl;
    std::cout << "Epic kernel" << std::endl;
    test<CGAL::Exact_predicates_inexact_constructions_kernel>(filename);
}

int main(void)
{
    std::cout << "AABB intersection tests" << std::endl;
    test_several_kernels("./data/cube.off");
    test_several_kernels("./data/coverrear.off");
    test_several_kernels("./data/nested_spheres.off");
    return 0;
}
