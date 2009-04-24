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

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Line_3 Line;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;

template <class Iterator>
class Triangle_primitive
{
        // type
public:
        typedef K::Triangle_3 Object; // object type
        typedef Iterator Id; // Id type

        // member data
private:
        Id m_it; // iterator
        Object m_object; // 3D triangle

public:
        Triangle_primitive(Id it)
                : m_it(it)
        {
                m_object = *it; // copy triangle
        }
public:
        const Object& object() const { return m_object; }
        Object& object() { return m_object; }
        Id id() { return m_it; }
};

typedef std::list<typename Triangle>::iterator Iterator;
typedef Triangle_primitive<Iterator> Primitive;
typedef CGAL::AABB_traits<typename K, typename Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<typename AABB_triangle_traits> Polyhedron_tree;

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

        Polyhedron_tree tree(triangles.begin(),triangles.end());
        Line line(a,b);
        std::cout << tree.number_of_intersections(line) 
                << " intersections(s) with ray" << std::endl;
        return 0;
}