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
// Author(s)     : Camille Wormser, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <fstream>
#include <CGAL/Timer.h>

#include <CGAL/AABB_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include <CGAL/spatial_sort.h>

#include "AABB_test_util.h"

#define NBQ 100000

template<class Value>
size_t check_outputs(const std::vector<Value>& a, const std::vector<Value>& b) {
        size_t counter = 0;
        for(size_t i = 0; i < a.size(); ++i) {
                if(a[i] != b[i])
                        ++counter;
        }
        return counter;
}

template <class Tree, class K>
void test_hint_strategies(Tree& tree, CGAL::Polyhedron_3<K>& polyhedron)
{
        typedef typename K::Point_3 Point;
        typedef typename Tree::Primitive::Id Id;
        typedef typename Tree::Point_and_primitive_id Point_and_primitive_id;

        std::vector<Point> queries;
        std::vector<Id> outputs1, outputs2, outputs3;
        
        queries.reserve(NBQ);
        outputs1.reserve(NBQ);
        outputs2.reserve(NBQ);
        outputs3.reserve(NBQ);
        
        size_t common_min = NBQ;
        size_t counter;
        
        for(size_t i = 0; i < NBQ; ++i)
                queries.push_back(random_point_in<K>(tree.bbox()));
        
        CGAL::spatial_sort(queries.begin(), queries.end());
        
        CGAL::Timer timer;
        timer.start();
        counter = 0;
        while(timer.time() < 1. && counter < NBQ) {
                outputs1.push_back(tree.closest_point_and_primitive(queries[counter]).second);
                ++counter;
        }
        timer.stop();
        double speed = static_cast<double>(counter)/(counter == NBQ? timer.time(): 1.); 
        std::cout << "without hint:      " << speed << " queries/s" << std::endl;
        timer.reset();
        
        Point_and_primitive_id hint = tree.any_reference_point_and_id();
        
        timer.start();
        counter = 0;
        while(timer.time() < 1. && counter < NBQ) {
                outputs2.push_back((hint = tree.closest_point_and_primitive(queries[counter], hint)).second);
                ++counter;
        }
        timer.stop();
        speed = static_cast<double>(counter)/(counter == NBQ? timer.time(): 1.); 
        std::cout << "with spatial sort: " << speed << " queries/s" << std::endl;
        timer.reset();        
      
        tree.accelerate_distance_queries(polyhedron.points_begin(),polyhedron.points_end());
        
        timer.start();   
        counter = 0;
        while(timer.time() < 1. && counter < NBQ) {
                outputs3.push_back(tree.closest_point_and_primitive(queries[counter]).second);
                ++counter;
        }
        timer.stop();
        speed = static_cast<double>(counter)/(counter == NBQ? timer.time(): 1.); 
        std::cout << "with KD-tree:      " << speed << " queries/s" << std::endl << std::endl;
        timer.stop();               
        std::cout << "Consistency:" << std::endl;
        if((counter = check_outputs(outputs1, outputs2)) == 0)
                std::cout << "         without hint and spatial sort are consistent" << std::endl;
        else
                std::cout << "WARNING, without hint and spatial sort have " << counter << " inconsistencies" << std::endl;
                
        if((counter = check_outputs(outputs1, outputs3)) == 0)
                std::cout << "         without hint and with KD-tree are consistent" << std::endl;
        else
                std::cout << "WARNING, without hint and with KD-tree have " << counter << " inconsistencies (hopefully just because the hint is returned)" << std::endl;
                
        std::cout << std::endl;
}


template <class Tree, class K>
void test_speed(Tree& tree)
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
                Point query = random_point_in<K>(tree.bbox());
                Point closest = tree.closest_point(query);
                nb++;
        }
        double speed = (double)nb / timer.time();
        std::cout << speed << " distance queries/s" << std::endl;
        timer.stop();
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

        Polyhedron polyhedron;
        std::ifstream ifs(filename);
        ifs >> polyhedron;

        // constructs AABB tree and internal search KD-tree with 
        // the points of the polyhedron
        Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());

        // tests hint strategies
        test_hint_strategies<Tree,K>(tree, polyhedron);
}

void test_kernels(const char *filename)
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
        std::cout << "AABB hint strategies tests" << std::endl;
        test_kernels("./data/cube.off");
        test_kernels("./data/coverrear.off");
        test_kernels("./data/nested_spheres.off");
        return 0;
}
