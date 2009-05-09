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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>


double random_in(const double a,
                 const double b)
{
    double r = rand() / (double)RAND_MAX;
    return a + (b - a) * r;
}

template <class K>
typename K::Point_3 random_point_in(const CGAL::Bbox_3& bbox)
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(bbox.xmin(),bbox.xmax());
    FT y = (FT)random_in(bbox.ymin(),bbox.ymax());
    FT z = (FT)random_in(bbox.zmin(),bbox.zmax());
    return typename K::Point_3(x,y,z);
}

template <class K>
typename K::Vector_3 random_vector()
{
    typedef typename K::FT FT;
    FT x = (FT)random_in(0.0,1.0);
    FT y = (FT)random_in(0.0,1.0);
    FT z = (FT)random_in(0.0,1.0);
    return typename K::Vector_3(x,y,z);
}

template <class Tree, class K>
void test_all_intersection_query_types(Tree& tree)
{
    std::cout << "Test all query types" << std::endl;

    typedef typename K::FT FT;
    typedef typename K::Ray_3 Ray;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typedef typename K::Vector_3 Vector;
    typedef typename K::Segment_3 Segment;
    typedef typename Tree::Primitive Primitive;
    typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;
    typedef typename Tree::Object_and_primitive_id Object_and_primitive_id;

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
    std::list<typename Primitive::Id> primitives;
    tree.all_intersected_primitives(ray,std::back_inserter(primitives));
    tree.all_intersected_primitives(line,std::back_inserter(primitives));
    tree.all_intersected_primitives(segment,std::back_inserter(primitives));

    // any_intersection
    boost::optional<Object_and_primitive_id> optional_object_and_primitive;
    optional_object_and_primitive = tree.any_intersection(ray);
    optional_object_and_primitive = tree.any_intersection(line);
    optional_object_and_primitive = tree.any_intersection(segment);

    // any_intersected_primitive
    boost::optional<typename Primitive::Id> optional_primitive;
    optional_primitive = tree.any_intersected_primitive(ray);
    optional_primitive = tree.any_intersected_primitive(line);
    optional_primitive = tree.any_intersected_primitive(segment);

    // all_intersections
    std::list<Object_and_primitive_id> intersections;
    tree.all_intersections(ray,std::back_inserter(intersections));
    tree.all_intersections(line,std::back_inserter(intersections));
    tree.all_intersections(segment,std::back_inserter(intersections));
}


template <class Tree, class K>
void test_all_distance_query_types(Tree& tree)
{
        typedef typename K::FT FT;
        typedef typename K::Ray_3 Ray;
        typedef typename K::Point_3 Point;
        typedef typename K::Vector_3 Vector;
        typedef typename Tree::Primitive Primitive;
        typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

        Point query = random_point_in<K>(tree.bbox());
        Point_and_primitive_id hint = tree.any_reference_point_and_id();

        FT sqd1 = tree.squared_distance(query);
        FT sqd2 = tree.squared_distance(query,hint.first);

        Point p1 = tree.closest_point(query);
        Point p2 = tree.closest_point(query,hint.first);

        Point_and_primitive_id pp1 = tree.closest_point_and_primitive(query);
        Point_and_primitive_id pp2 = tree.closest_point_and_primitive(query,hint);
}

template <class Tree, class K>
void test_distance_speed(Tree& tree)
{
        typedef typename K::FT FT;
        typedef typename K::Ray_3 Ray;
        typedef typename K::Point_3 Point;
        typedef typename K::Vector_3 Vector;

        CGAL::Timer timer;
        timer.start();
        unsigned int nb = 0;
        while(timer.time() < 1.0)
        {
                // picks a random point in the tree bbox
                Point query = random_point_in<K>(tree.bbox());
                Point closest = tree.closest_point(query);
                nb++;
        }
        double speed = (double)nb / timer.time();
        std::cout << speed << " distance queries/s" << std::endl;
        timer.stop();
}



