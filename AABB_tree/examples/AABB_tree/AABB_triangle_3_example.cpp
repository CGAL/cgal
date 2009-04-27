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
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;

typedef std::list<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

int main(void)
{
        Point a(1.0, 0.0, 0.0);
        Point b(0.0, 1.0, 0.0);
        Point c(0.0, 0.0, 1.0);
        Point d(0.0, 0.0, 0.0);

        std::list<Triangle> triangles;
        triangles.push_back(Triangle(a,b,c));
        triangles.push_back(Triangle(a,b,d));
        triangles.push_back(Triangle(a,d,c));

        Tree tree(triangles.begin(),triangles.end());

        // counts #intersections
        Ray ray(a,b);
        std::cout << tree.number_of_intersections(ray)
                << " intersections(s) with ray" << std::endl;

        // compute closest point
        Point hint(a);
        Point query(2.0, 2.0, 2.0);
        Point closest = tree.closest_point(query,hint);

        return 0;
}