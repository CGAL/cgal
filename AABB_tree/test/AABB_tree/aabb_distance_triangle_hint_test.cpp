// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/spatial_sort.h>

#include "AABB_test_util.h"

#define NBQ 100000

template<class Value>
size_t check_outputs(const std::vector<Value>& a,
                     const std::vector<Value>& b,
                     Value deflt)
{
    typedef typename std::vector<Value>::const_iterator Iterator;

    size_t counter = 0;

    Iterator it_a = a.begin();
    Iterator it_b = b.begin();

    while ( it_a != a.end() && it_b != b.end() )
    {
      if ( *it_a != *it_b && *it_b != deflt && *it_a != deflt )
        ++counter;

      ++it_a;
      ++it_b;
    }

    return counter;
}

template <class Tree, class K>
void test_hint_strategies(Tree& tree,
                          CGAL::Polyhedron_3<K>& polyhedron,
                          const double duration)
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

        size_t counter;

        for(size_t i = 0; i < NBQ; ++i)
                queries.push_back(random_point_in<K>(tree.bbox()));

        CGAL::spatial_sort(queries.begin(), queries.end());

        CGAL::Timer timer;
        timer.start();
        counter = 0;
        while(timer.time() < duration && counter < NBQ) {
                outputs1.push_back(tree.closest_point_and_primitive(queries[counter]).second);
                ++counter;
        }
        timer.stop();
        double speed = static_cast<double>(counter)/(counter == NBQ? timer.time(): duration);
        std::cout << "without hint:      " << speed << " queries/s" << std::endl;
        timer.reset();

        Point_and_primitive_id hint = tree.any_reference_point_and_id();

        timer.start();
        counter = 0;
        while(timer.time() < duration && counter < NBQ) {
                outputs2.push_back((hint = tree.closest_point_and_primitive(queries[counter], hint)).second);
                ++counter;
        }
        timer.stop();
        speed = static_cast<double>(counter)/(counter == NBQ? timer.time(): duration);
        std::cout << "with spatial sort: " << speed << " queries/s" << std::endl;
        timer.reset();

        tree.accelerate_distance_queries(polyhedron.points_begin(),polyhedron.points_end());

        timer.start();
        counter = 0;
        while(timer.time() < duration && counter < NBQ) {
                outputs3.push_back(tree.closest_point_and_primitive(queries[counter]).second);
                ++counter;
        }
        timer.stop();
        speed = static_cast<double>(counter)/(counter == NBQ? timer.time(): duration);
        std::cout << "with KD-tree:      " << speed << " queries/s" << std::endl << std::endl;
        std::cout << "Consistency:" << std::endl;
        if((counter = check_outputs(outputs1, outputs2, Id())) == 0)
                std::cout << "without hint and spatial sort are consistent" << std::endl;
        else
                std::cout << "without hint and spatial sort have " << counter << " inconsistencies (closest point on vertex/edge?)" << std::endl;

        if((counter = check_outputs(outputs1, outputs3, Id())) == 0)
                std::cout << "without hint and with KD-tree are consistent (modulo hint case)" << std::endl;
        else
                std::cout << "without hint and with KD-tree have " << counter << " inconsistencies (closest point on vertex/edge? the hint case has been excluded)" << std::endl;

        std::cout << std::endl;
}

template<class K, class Tree, class Polyhedron, Primitive_type Type>
void test_impl(Tree& tree,
               Polyhedron& p,
               const double duration)
{
  test_hint_strategies<Tree,K>(tree, p, duration);
}

int main(void)
{
    const double duration = 0.1;
    std::cout << "AABB hint strategy tests" << std::endl;
    test_kernels<TRIANGLE>("data/cube.off",duration);
    test_kernels<TRIANGLE>("data/coverrear.off",duration);
    test_kernels<TRIANGLE>("data/finger.off",duration);
    return EXIT_SUCCESS;
}
