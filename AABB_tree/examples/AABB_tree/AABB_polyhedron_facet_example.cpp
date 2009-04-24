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
#include <Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_tree.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;

typedef K::FT FT;
typedef K::Ray_3 Ray;
typedef K::Point_3 Point;

class Polyhedron_triangle_primitive
{
        // type
public:
        typedef K::Triangle_3 Object; // object type
        typedef Polyhedron::Facet_handle Id; // Id type

        // member data
private:
        Id m_handle; // Facet_handle
        Object m_object; // 3D triangle

        Polyhedron_triangle_primitive(Id handle)
                : m_handle(handle)
        {
                const Point& a = handle->halfedge()->vertex()->point();
                const Point& b = handle->halfedge()->next()->vertex()->point();
                const Point& c = handle->halfedge()->prev()->vertex()->point();
                m_object = Object(a,b,c);
        }
public:
        Object object() { return m_object; }
        Id id() { return m_handle; }
};

typedef AABB_traits<K, Polyhedron_triangle_primitive> AABB_Polyhedron_traits;
typedef AABB_tree<AABB_Polyhedron_traits> Polyhedron_tree;

int main(void)
{
        Point p(1.0, 0.0, 0.0);
        Point q(0.0, 1.0, 0.0);
        Point r(0.0, 0.0, 1.0);
        Point s(0.0, 0.0, 0.0);
        Polyhedron polyhedron;
        polyhedron.make_tetrahedron( p, q, r, s);

        Polyhedron_tree tree(polyhedron.facets_begin(),polyhedron.facets_end());
        Point source(0.2, 0.2, 0.2);
        Ray ray(source, Vector(0.1, 0.2, 0.3));
        std::cout << tree.number_of_intersections(ray) 
                << " intersections(s) with ray" << std::endl;
        return 0;
}