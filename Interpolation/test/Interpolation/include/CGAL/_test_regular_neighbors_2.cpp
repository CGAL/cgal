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
// Author(s)     : Julia Floetotto

#include <iostream>
#include <cassert>
#include <utility>

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/Origin.h>

#include <CGAL/regular_neighbor_coordinates_2.h>

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
		     ::value_type::second_type norm, const Point& p)
{
  typedef typename CGAL::Kernel_traits<Point>::Kernel::Point_2  Bare_point;
  Bare_point b = CGAL::ORIGIN;
  for(; first != beyond; ++first)
    b = b + (first->second/norm) * (first->first.point() - CGAL::ORIGIN);

  return p==b;
}


template <class Triangul>
void
_test_regular_neighbors_2( const Triangul & )
{
  Triangul T;

  int n=20, m=200;
  double r = 3;
  double max_weight =1;

  typedef typename Triangul::Geom_traits          Gt;
  typedef typename Gt::Weighted_point_2           Weighted_point;
  typedef typename Gt::Bare_point                 Bare_point;
  typedef typename Gt::Rep::FT                    Coord_type;

  typedef std::vector< std::pair< Weighted_point, Coord_type > >
                                                  Point_coordinate_vector;


  std::cout << "RN2: Testing random points." << std::endl;
  //test random points in a square of length r:
  std::vector<Bare_point> points;
  points.reserve(n+m);

  //put four bounding box points:
  points.push_back(Bare_point(-r,-r));
  points.push_back(Bare_point(r,-r));
  points.push_back(Bare_point(-r,r));
  points.push_back(Bare_point(r,r));

  // Create n+m-4 points within a disc of radius 2
  CGAL::Random_points_in_square_2<Bare_point> g(r);
  CGAL::cpp11::copy_n( g, n+m, std::back_inserter(points));


  CGAL::Random random;
  for(int i=0; i<n ; i++)
    T.insert(Weighted_point(points[i],random.get_double(-max_weight,
							max_weight)));


  Point_coordinate_vector coords;
  for(int i=n;i<n+m;i++){
    Weighted_point wp = Weighted_point(points[i],random.get_double
				       (2*max_weight,3*max_weight));
    CGAL::Triple<
      std::back_insert_iterator<Point_coordinate_vector>,
      Coord_type, bool> coordinate_result =
      CGAL::regular_neighbor_coordinates_2(T,wp,std::back_inserter(coords));
    assert(coordinate_result.third);
    Coord_type norm = coordinate_result.second;

    assert(norm>0);
    assert(test_barycenter( coords.begin(), coords.end(),norm, points[i]));
    coords.clear();
  }

  //TESTING a GRID POINT SET
  std::cout << "RN2: Testing grid points." << std::endl;

  Triangul T2;
  //Grid points:
  Bare_point p1_2(-2, -2);
  Bare_point p2_2(-2,2);
  Bare_point p3_2(2,-2);
  Bare_point p4_2(2,2);
  Bare_point p1(-1, -1);
  Bare_point p2(-1,1);
  Bare_point p3(1,-1);
  Bare_point p4(1,1);
  Bare_point p12(-1, 0);
  Bare_point p23(0,1);
  Bare_point p34(1,0);
  Bare_point p41(0,-1);

  T2.insert(Weighted_point(p1_2, 0));
  T2.insert(Weighted_point(p2_2,  0));
  T2.insert(Weighted_point(p3_2, 0));
  T2.insert(Weighted_point(p4_2, 0));
  T2.insert(Weighted_point(p1, 0));
  T2.insert(Weighted_point(p2,  0));
  T2.insert(Weighted_point(p3, 0));
  T2.insert(Weighted_point(p4, 0));
  T2.insert(Weighted_point(p12,  0));
  T2.insert(Weighted_point(p23,  0));
  T2.insert(Weighted_point(p34,   0));
  T2.insert(Weighted_point(p41,  0));

  //test with 0 weight:
  Weighted_point wp(Bare_point(0,0),0);
  CGAL::Triple<
    std::back_insert_iterator<Point_coordinate_vector>,
    Coord_type, bool> coordinate_result =
    CGAL::regular_neighbor_coordinates_2(T2,wp,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  Coord_type norm = coordinate_result.second;
  assert(norm == Coord_type(1));
  typename std::vector< std::pair< Weighted_point, Coord_type >
    >::const_iterator
    ci= coords.begin();
  for(; ci!= coords.end(); ci++)
    assert(ci->second == Coord_type(0.25));
  assert(test_barycenter( coords.begin(), coords.end(),norm, wp));
  coords.clear();

  //test with hidden_vertices:
  wp = Weighted_point(Bare_point(0,0),4);
  coordinate_result =
    CGAL::regular_neighbor_coordinates_2(T2,wp,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm,wp));
  coords.clear();

  //add the middle point of the grid
  T2.insert(Weighted_point(Bare_point(0,0),0));

  //point on a vertex;
  wp = Weighted_point(p34,0);
  coordinate_result  =
    CGAL::regular_neighbor_coordinates_2(T2,wp,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(norm == Coord_type(1));
  ci= coords.begin();
  assert(ci->first == wp);
  assert(ci->second == Coord_type(1));
  ci++;
  assert(ci==coords.end());
  coords.clear();

  //point on the vertex but creating a hole:
  wp = Weighted_point(p34,2);
  coordinate_result =
    CGAL::regular_neighbor_coordinates_2(T2,wp,
					 std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm, wp));
  coords.clear();

  //point on an edge:
  wp= Weighted_point(Bare_point(0,0.5), 3);
  coordinate_result = CGAL::regular_neighbor_coordinates_2
    (T2,wp,std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm, wp));
  coords.clear();

  //a vertex v in Reg(P\v->point()):
  typename Triangul::Vertex_iterator vit = T2.finite_vertices_end();
  coordinate_result = CGAL::regular_neighbor_coordinates_2
    (T2, --vit, std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  assert(test_barycenter( coords.begin(), coords.end(),norm,vit->point()));
  coords.clear();

  //outside convex hull:
  wp= Weighted_point(Bare_point(3,0.5), 3);
  coordinate_result = CGAL::regular_neighbor_coordinates_2
    (T2,wp,std::back_inserter(coords));
  assert(!coordinate_result.third);

  //on a convex hull edge:
  wp= Weighted_point(Bare_point(2,1), 3);
  coordinate_result = CGAL::regular_neighbor_coordinates_2
    (T2,wp,std::back_inserter(coords));
  assert(!coordinate_result.third);
}
