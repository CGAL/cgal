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
// $URL:  $
// $Id:  $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Object_and_primitive_id Object_and_primitive_id;

int main()
{
        Point p(1.0, 0.0, 0.0);
        Point q(0.0, 1.0, 0.0);
        Point r(0.0, 0.0, 1.0);
        Point s(0.0, 0.0, 0.0);
        Polyhedron polyhedron;
        polyhedron.make_tetrahedron(p, q, r, s);

        // constructs AABB tree
        Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());

        // computes any intersection with segment query
        Point a(0.2, 0.2, 0.2);
        Point b(1.3, 1.3, 1.3);
        Segment query(a,b);
        boost::optional<Object_and_primitive_id> intersection = 
			tree.any_intersection(query);

		if(intersection)
		{
			// get intersection object
			Object_and_primitive_id op = *intersection;
			CGAL::Object object = op.first;
			Point point;
			if(CGAL::assign(point,object))
			{
				std::cout << "intersection object is a point" << std::endl;
			}
		}

        return EXIT_SUCCESS;
}
