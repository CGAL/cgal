// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Naceur MESKINI.
//                 Mael Rouxel-Labbé

#include <CGAL/Interpolation/internal/helpers.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/barycenter.h>

#include <cassert>
#include <iostream>
#include <iterator>
#include <list>
#include <utility>
#include <vector>

template < class ForwardIterator >
bool test_norm(ForwardIterator first, ForwardIterator beyond,
               typename std::iterator_traits<ForwardIterator>::value_type::second_type norm)
{
  typename std::iterator_traits<ForwardIterator>::value_type::second_type sum(0);
  for(; first!=beyond; first++)
    sum += first->second;

  return norm == sum;
}

template < class ForwardIterator, class Dt >
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
                     const typename Dt::Point& p)
{
  return p == CGAL::barycenter(first, beyond);
}

template < class Dt >
void _test_natural_neighbors_2_without_outputfunctor(const Dt& T)
{
  std::cout << "Testing backward compatibility..." << std::endl;

  typedef typename Dt::Geom_traits                          Gt;
  typedef typename Gt::Point_2                              Point_2;
  typedef typename Gt::FT                                   Coord_type;

  typedef std::vector<std::pair<Point_2, Coord_type> >      Point_coordinate_vector;
  typedef typename Point_coordinate_vector::const_iterator  PCV_cit;

  typedef typename Dt::Vertex_handle                        Vertex_handle;
  typedef typename Dt::Face_handle                          Face_handle;
  typedef std::pair<Face_handle, int>                       Edge;

  Point_coordinate_vector coords;

  // point on a vertex
  Point_2 p = T.finite_vertices_begin()->point();
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> coordinate_result =
      CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords));
  assert(coordinate_result.third);

  Coord_type norm = coordinate_result.second;
  assert(norm == Coord_type(1));

  PCV_cit ci = coords.begin();
  assert(ci->first == p);
  assert(ci->second == Coord_type(1));
  ++ci;
  assert(ci == coords.end());
  coords.clear();

  // On a point, but as a vertex handle
  Vertex_handle vh = --(T.finite_vertices_end());
  coordinate_result = CGAL::natural_neighbor_coordinates_2(T, vh, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  bool is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);
  assert(norm == Coord_type(1));

  ci = coords.begin();
  while(ci != coords.end())
  {
    assert(ci->first != vh->point());
    assert(ci->second != Coord_type(1));
    ++ci;
  }
  coords.clear();

  // point on an edge
  p = Point_2(0,0.5);
  coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_barycenter<PCV_cit, Dt>(coords.begin(), coords.end(), p);
  assert(is_equal);
  coords.clear();

  // with the conflict hole
  std::list<Edge> hole_bd;
  T.get_boundary_of_conflicts(p, std::back_inserter(hole_bd));
  coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords),
                                                           hole_bd.begin(), hole_bd.end());
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<PCV_cit, Dt>(coords.begin(), coords.end(), p);
  assert(is_equal);
  coords.clear();

  // outside convex hull
  p = Point_2(3,0.5);
  coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords));
  assert(!coordinate_result.third);

  // same, but with the face hint API
  coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords),
                                                           T.finite_faces_begin());
  assert(!coordinate_result.third);

  // on a convex hull edge:
  coords.clear();
  p = Point_2(2, 1);
  coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<PCV_cit, Dt>(coords.begin(), coords.end(), p);
  assert(is_equal);
}

template < class Dt >
void _test_natural_neighbors_2_with_outputfunctor(const Dt& T)
{
  std::cout << "Testing with OutputFunctor..." << std::endl;

  typedef typename Dt::Geom_traits                               Gt;
  typedef typename Dt::Point                                     Point;
  typedef typename Gt::FT                                        Coord_type;

  typedef typename Dt::Vertex_handle                             Vertex_handle;
  typedef typename Dt::Face_handle                               Face_handle;
  typedef std::pair<Face_handle, int>                            Edge;

  typedef std::vector<std::pair<Vertex_handle, Coord_type> >     Vertex_coordinate_vector;
  typedef typename Vertex_coordinate_vector::const_iterator      VCV_cit;
  typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> >  Identity_output_functor;

  typedef std::vector<std::pair<Point, Coord_type> >                            Point_coordinate_vector;
  typedef typename Point_coordinate_vector::const_iterator                      PCV_cit;
  typedef CGAL::Interpolation::internal::Extract_point_in_pair<Dt, Coord_type>  Point_output_functor;

  Vertex_coordinate_vector vh_coords;
  Identity_output_functor vh_fct;

  Point_coordinate_vector pt_coords;
  Point_output_functor pt_fct;

  // point on a vertex
  Point p = T.finite_vertices_begin()->point();
  CGAL::Triple<std::back_insert_iterator<Vertex_coordinate_vector>, Coord_type, bool> vh_coordinate_result =
      CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(vh_coords), vh_fct);
  assert(vh_coordinate_result.third);

  Coord_type norm = vh_coordinate_result.second;
  assert(norm == Coord_type(1));

  VCV_cit ci = vh_coords.begin();
  assert(ci->first == Vertex_handle(T.finite_vertices_begin()));
  assert(ci->second == Coord_type(1));
  ++ci;
  assert(ci == vh_coords.end());
  vh_coords.clear();

  // On a point, but as a vertex handle
  Vertex_handle vh = --(T.finite_vertices_end());
  vh_coordinate_result = CGAL::natural_neighbor_coordinates_2(T, vh, std::back_inserter(vh_coords), vh_fct);
  assert(vh_coordinate_result.third);

  norm = vh_coordinate_result.second;
  assert(norm == Coord_type(1));

  ci = vh_coords.begin();
  while(ci != vh_coords.end())
  {
    assert(ci->first != vh);
    assert(T.tds().is_vertex(vh));
    assert(ci->second != Coord_type(1));
    ++ci;
  }
  vh_coords.clear();

  // point on an edge
  p = Point(0,0.5);
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> pt_coordinate_result =
    CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  bool is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<PCV_cit, Dt>(pt_coords.begin(), pt_coords.end(), p);
  assert(is_equal);
  pt_coords.clear();

  // with the conflict hole
  std::list<Edge> hole_bd;
  T.get_boundary_of_conflicts(p, std::back_inserter(hole_bd));
  pt_coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(pt_coords), pt_fct,
                                                              hole_bd.begin(), hole_bd.end());
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<PCV_cit, Dt>(pt_coords.begin(), pt_coords.end(), p);
  assert(is_equal);
  pt_coords.clear();

  // outside convex hull
  p = Point(3,0.5);
  pt_coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(pt_coords), pt_fct);
  assert(!pt_coordinate_result.third);
  pt_coords.clear();

  // same, but with the face hint API
  pt_coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(pt_coords), pt_fct,
                                                              T.finite_faces_begin());
  assert(!pt_coordinate_result.third);
  pt_coords.clear();

  // on a convex hull edge:
  p = Point(2, 1);
  pt_coordinate_result = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<PCV_cit, Dt>(pt_coords.begin(), pt_coords.end(), p);
  assert(is_equal);
}

template < class Dt >
void _test_natural_neighbors_2(const Dt&)
{
  typedef typename Dt::Geom_traits                Gt;
  typedef typename Gt::Point_2                    Point;

  Dt T;

  //Grid points:
  Point p1_2(-2, -2);
  Point p2_2(-2, 2);
  Point p3_2(2, -2);
  Point p4_2(2, 2);
  Point p1(-1, -1);
  Point p2(-1, 1);
  Point p3(1, -1);
  Point p4(1, 1);
  Point p12(-1, 0);
  Point p23(0, 1);
  Point p34(1, 0);
  Point p41(0, -1);

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
  T.insert(Point(0,0));

  _test_natural_neighbors_2_without_outputfunctor(T);
  _test_natural_neighbors_2_with_outputfunctor(T);

  std::cout << "done" << std::endl;
}
