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

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

template <class Tree, class K>
void test_all_query_types(Tree& tree)
{
    std::cout << "Test all query types" << std::endl;

    typedef K::FT FT;
    typedef K::Ray_3 Ray;
    typedef K::Line_3 Line;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;
    typedef K::Segment_3 Segment;
    typedef Tree::Primitive Primitive;
    typedef Tree::Intersection Intersection;

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

    // number_of_intersections
    tree.number_of_intersections(ray); 
    tree.number_of_intersections(line); 
    tree.number_of_intersections(segment); 

    // all_intersected_primitives
    std::list<Primitive> primitives;
    tree.all_intersected_primitives(ray,std::back_inserter(primitives));
    tree.all_intersected_primitives(line,std::back_inserter(primitives));
    tree.all_intersected_primitives(segment,std::back_inserter(primitives));

    // any_intersection
    Intersection intersection;
    success = tree.any_intersection(ray,intersection);
    success = tree.any_intersection(line,intersection);
    success = tree.any_intersection(segment,intersection);

    // all_intersections
    std::list<Intersection> intersections;
    tree.all_intersections(ray,std::back_inserter(intersections));
    tree.all_intersections(line,std::back_inserter(intersections));
    tree.all_intersections(segment,std::back_inserter(intersections));
}

template <class Tree, class Polyhedron, class K>
void test_speed(Tree& tree, 
                Polyhedron& polyhedron)
{
    std::cout << "Test for speed" << std::endl;
    typedef K::FT FT;
    typedef K::Ray_3 Ray;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;

    CGAL::Timer timer;
    unsigned int nb = 0;
    timer.start();
    Point source((FT)0.0, (FT)0.0, (FT)0.0);
    Vector vec((FT)0.1, (FT)0.2, (FT)0.3);
    Ray ray(source, vec);
    while(timer.time() < 1.0)
    {
        tree.do_intersect(ray); 
        nb++;
    }
    double speed = (double)nb / timer.time();
    std::cout << speed << " intersections/s" << std::endl;
    timer.stop();
}

template <class K>
void test(const char *filename)
{
    typedef K::FT FT;
    typedef K::Ray_3 Ray;
    typedef K::Point_3 Point;
    typedef K::Vector_3 Vector;
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<K, Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;

    // load (triangle) polyhedral surface
    Polyhedron polyhedron;
    std::ifstream ifs(filename);
    ifs >> polyhedron;

    // construct tree (without internal KD-tree as we do not query any projection).
    Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());

    // call tests
    test_all_query_types<Tree,K>(tree);
    test_speed<Tree,Polyhedron,K>(tree,polyhedron);
}

void test_several_kernels(const char *filename)
{
    std::cout << std::endl;
    std::cout << "Polyhedron " << filename << std::endl;

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
    test_several_kernels("../data/cube.off");
    test_several_kernels("../data/coverrear.off");
    test_several_kernels("../data/nested_spheres.off");
    return 0;
}