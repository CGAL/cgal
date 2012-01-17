// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Naceur MESKINI.

#include <iostream>
#include <cassert>
#include <utility>

#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/barycenter.h>

#include <CGAL/natural_neighbor_coordinates_2.h>

template < class ForwardIterator >
bool test_norm(ForwardIterator first, ForwardIterator beyond,
	       typename std::iterator_traits<ForwardIterator>
	       ::value_type::second_type norm)
{
   typename
     std::iterator_traits<ForwardIterator>::value_type::second_type sum(0);
   for(; first !=beyond; first++)
     sum+= first->second;

  return norm==sum;
}

template < class ForwardIterator, class Point >
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
		     typename std::iterator_traits<ForwardIterator>
		     ::value_type::second_type /*norm*/, const Point& p)
{
  return p == CGAL::barycenter(first, beyond);
}


template <class Triangul>
void
_test_natural_neighbors_2( const Triangul & )
{

  typedef typename Triangul::Geom_traits          Gt;

  typedef typename Gt::Point_2                    Point_2;
  typedef typename Gt::FT                         Coord_type;

  typedef std::vector< std::pair<Point_2, Coord_type> > Point_coordinate_vector;

  //TESTING a GRID POINT SET
  std::cout << "NN2: Testing grid points." << std::endl;

  Triangul T;
  //Grid points:
  Point_2 p1_2(-2, -2);
  Point_2 p2_2(-2,2);
  Point_2 p3_2(2,-2);
  Point_2 p4_2(2,2);
  Point_2 p1(-1, -1);
  Point_2 p2(-1,1);
  Point_2 p3(1,-1);
  Point_2 p4(1,1);
  Point_2 p12(-1, 0);
  Point_2 p23(0,1);
  Point_2 p34(1,0);
  Point_2 p41(0,-1);

  T.insert(p1_2);
  T.insert(p2_2);
  T.insert(p3_2);
  T.insert(p4_2);
  T.insert(p1);
  T.insert(p2);
  T.insert(p3);
  T.insert(p4);
  T.insert(p12);
  T.insert(p23);
  T.insert(p34);
  T.insert(p41);
  T.insert(Point_2(0,0));


  Point_coordinate_vector coords;

  //point on a vertex;
  Point_2 p = Point_2(p34);
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>,Coord_type, bool> coordinate_result  = CGAL::natural_neighbor_coordinates_2(T,p,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  Coord_type norm = coordinate_result.second;
  assert(norm == Coord_type(1));
   typename std::vector< std::pair< Point_2, Coord_type >
    >::const_iterator ci= coords.begin();
  assert(ci->first == p);
  assert(ci->second == Coord_type(1));
  ci++;
  assert(ci==coords.end());
  coords.clear();

  //point on an edge:
  p = Point_2(0,0.5);
  coordinate_result = CGAL::natural_neighbor_coordinates_2
    (T,p,std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm, p));
  coords.clear();
  
  //outside convex hull:
  p= Point_2(3,0.5);
  coordinate_result = CGAL::natural_neighbor_coordinates_2
    (T,p,std::back_inserter(coords));
  assert(!coordinate_result.third);

  //on a convex hull edge:
  coords.clear();
  p= Point_2(2,1);
  coordinate_result = CGAL::natural_neighbor_coordinates_2
    (T,p,std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm, p));

}
