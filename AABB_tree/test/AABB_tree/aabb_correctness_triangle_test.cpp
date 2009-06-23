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

template <class K>
void test()
{
	// load polyhedron
    typedef CGAL::Polyhedron_3<K> Polyhedron;
    Polyhedron polyhedron;
    std::ifstream ifs("./data/tetrahedron.off");
    ifs >> polyhedron;

	typedef K::Point_3 Point;
	typedef K::Segment_3 Segment;

	// construct tree
    typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
    typedef CGAL::AABB_traits<K,Primitive> Traits;
    typedef CGAL::AABB_tree<Traits> Tree;
	typedef Tree::Object_and_primitive_id Object_and_primitive_id;
	Tree tree(polyhedron.facets_begin(),polyhedron.facets_end());

	// segment query
	Point p(0.25, -0.25, 0.25);
	Point q(0.25,  0.25, 0.25);
	Segment pq(p,q);
	assert(tree.do_intersect(pq));
	assert(tree.number_of_intersected_primitives(pq) == 1);

	boost::optional<Object_and_primitive_id> any;
	any = tree.any_intersection(pq);
	assert(any);
	Object_and_primitive_id op = *any;
	CGAL::Object object = op.first;
	Point point;
	assert(CGAL::assign(object,point));
	std::cout << "Intersection point: (" << point.x() <<
		                             ";" << point.y() << 
									 ";" << point.z() << ")" << std::endl;
}

int main()
{
	test<CGAL::Simple_cartesian<float> >();
	test<CGAL::Simple_cartesian<double> >();
	test<CGAL::Exact_predicates_inexact_constructions_kernel>();

    return EXIT_SUCCESS;
}
