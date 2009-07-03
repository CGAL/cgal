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
// Author(s)     : Camille Wormser, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <list>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<double> K;


// My own types:
struct My_point {
  double x;
  double y;
  double z;
  
  My_point (double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};

struct My_triangle {
  My_point a;
  My_point b;
  My_point c;
  
  My_triangle (My_point _a, My_point _b, My_point _c) : a(_a), b(_b), c(_c) {}
};

// My triangles will be stored in a vector
typedef std::vector<My_triangle>::const_iterator Iterator;

// The following primitive provides the conversion facilities between
// my own triangle and point types and the CGAL ones
struct My_triangle_primitive {
public:
// this is the type of data that the queries will return. For this example, 
// we imagine that, for some reason, we do not want to store the iterators
// of the vector, but raw pointers. This is just to show that the Id type
// does not have to be the same as the one of the input parameter of the 
// constructor.
  typedef const My_triangle*    Id;
  
// the CGAL types that are returned
  typedef K::Point_3      Point;
  typedef K::Triangle_3   Datum;
private:
  Id m_pt;  // this is what the AABB tree will store internally

public:
  My_triangle_primitive() {} // default constructor is needed
  
  // the following constructor is the one that receives the iterators from the 
  // iterator range given as input to the AABB_tree
  My_triangle_primitive(Iterator a) : m_pt(&(*a)) {}
  
  const Id& id() const { return m_pt; }
  // on the fly conversion from the internal data to the CGAL types
  Datum datum() const { return Datum(Point(m_pt->a.x, m_pt->a.y, m_pt->a.z), 
                                     Point(m_pt->b.x, m_pt->b.y, m_pt->b.z),
                                     Point(m_pt->c.x, m_pt->c.y, m_pt->c.z)); }
  Point reference_point() const { return Point(m_pt->a.x, m_pt->a.y, m_pt->a.z); }
};



typedef CGAL::AABB_traits<K, My_triangle_primitive> My_AABB_traits;
typedef CGAL::AABB_tree<My_AABB_traits> Tree;

int main()
{
	My_point a(1.0, 0.0, 0.0);
	My_point b(0.0, 1.0, 0.0);
	My_point c(0.0, 0.0, 1.0);
	My_point d(0.0, 0.0, 0.0);

	std::vector<My_triangle> triangles;
	triangles.push_back(My_triangle(a,b,c));
	triangles.push_back(My_triangle(a,b,d));
	triangles.push_back(My_triangle(a,d,c));

	// constructs AABB tree
	Tree tree(triangles.begin(),triangles.end());

	// counts #intersections
	K::Ray_3 ray_query(K::Point_3(1.0, 0.0, 0.0), K::Point_3(0.0, 1.0, 0.0));
	std::cout << tree.number_of_intersected_primitives(ray_query)
		<< " intersections(s) with ray query" << std::endl;

	// compute closest point and squared distance
	K::Point_3 point_query(2.0, 2.0, 2.0);
	K::Point_3 closest_point = tree.closest_point(point_query);
	double sqd = tree.squared_distance(point_query);
	std::cout << "squared distance: " << sqd << std::endl;

	return EXIT_SUCCESS;
}
