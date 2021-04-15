// Copyright (c) 2003   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Julia Floetotto

#include <CGAL/regular_neighbor_coordinates_2.h>
#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
#include <CGAL/sibson_gradient_fitting.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/Origin.h>

#include <iostream>
#include <cassert>
#include <utility>

template < class ForwardIterator >
bool test_norm(ForwardIterator first, ForwardIterator beyond,
               typename std::iterator_traits<ForwardIterator>::value_type::second_type norm)
{
   typename std::iterator_traits<ForwardIterator>::value_type::second_type sum(0);
   for(; first !=beyond; first++)
     sum += first->second;

  return norm == sum;
}

template < class Rt, class ForwardIterator >
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
                     typename std::iterator_traits<ForwardIterator>::value_type::second_type norm,
                     const typename Rt::Weighted_point& p)
{
  typedef typename Rt::Geom_traits                          Gt;
  typedef typename Rt::Bare_point                           Bare_point;

  typename Gt::Construct_point_2 cp = Gt().construct_point_2_object();

  Bare_point b = CGAL::ORIGIN;
  for(; first != beyond; ++first)
    b = b + (first->second / norm) * (cp(first->first) - CGAL::ORIGIN);

  return (cp(p) == b);
}

template < class Rt >
void _test_natural_neighbors_2_without_outputfunctor(Rt T) // intentional copy because T is modified
{
  std::cout << "Testing backward compatibility..." << std::endl;

  typedef typename Rt::Bare_point                              Bare_point;
  typedef typename Rt::Weighted_point                          Weighted_point;
  typedef typename Rt::Geom_traits                             Gt;
  typedef typename Gt::FT                                      Coord_type;
  typedef std::vector<std::pair<Weighted_point, Coord_type> >  Point_coordinate_vector;

  typedef typename Rt::Vertex_handle                           Vertex_handle;
  typedef typename Rt::Face_handle                             Face_handle;
  typedef std::pair<Face_handle, int>                          Edge;

  Point_coordinate_vector coords;

  // test with 0 weight:
  Weighted_point wp(Bare_point(0,0), 0.);
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> coordinate_result =
    CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(coordinate_result.third);

  Coord_type norm = coordinate_result.second;
  assert(norm == Coord_type(1));

  typename std::vector< std::pair<Weighted_point, Coord_type> >::const_iterator ci = coords.begin();
  for(; ci!= coords.end(); ci++)
    assert(ci->second == Coord_type(0.25));

  bool is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, wp);
  assert(is_equal);
  coords.clear();

  // test with hidden_vertices:
  wp = Weighted_point(Bare_point(0,0), 4.);
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(coordinate_result.third);
  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, wp);
  assert(is_equal);
  coords.clear();

  // add the middle point of the grid
  T.insert(Weighted_point(Bare_point(0,0), 0.));

  // point on a vertex;
  wp = Weighted_point(Bare_point(1,0), 0.);
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  assert(norm == Coord_type(1));
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  ci = coords.begin();
  assert(ci->first == wp);
  assert(ci->second == Coord_type(1));
  ++ci;
  assert(ci == coords.end());
  coords.clear();

  // point on the vertex but creating a hole:
  wp = Weighted_point(Bare_point(1,0), 2.);
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, wp);
  assert(is_equal);
  coords.clear();

  // point on an edge:
  wp = Weighted_point(Bare_point(0,0.5), 3.);
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, wp);
  assert(is_equal);
  coords.clear();

  // a vertex v in Reg(P\v->point()):
  typename Rt::Vertex_iterator vit = T.finite_vertices_end();
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, --vit, std::back_inserter(coords));
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, vit->point());
  assert(is_equal);
  coords.clear();

  ci = coords.begin();
  while(ci != coords.end())
  {
    assert(ci->first != (--vit)->point());
    assert(ci->second != Coord_type(1));
    ++ci;
  }
  coords.clear();

  // with the conflict hole
  std::list<Edge> hole_bd;
  std::list<Vertex_handle> hidden_vertices;
  T.get_boundary_of_conflicts_and_hidden_vertices(wp, std::back_inserter(hole_bd), std::back_inserter(hidden_vertices));
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords),
                                                           hole_bd.begin(), hole_bd.end(),
                                                           hidden_vertices.begin(), hidden_vertices.end());
  assert(coordinate_result.third);

  norm = coordinate_result.second;
  is_equal = test_norm(coords.begin(), coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(coords.begin(), coords.end(), norm, wp);
  assert(is_equal);
  coords.clear();

  // outside convex hull:
  wp = Weighted_point(Bare_point(3,0.5), 3.);
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(!coordinate_result.third);

  // same but with the face hint API
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords),
                                                           T.finite_faces_begin());
  assert(!coordinate_result.third);

  // on a convex hull edge:
  wp = Weighted_point(Bare_point(2,1), 3.);
  coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(coords));
  assert(!coordinate_result.third);
}

template < class Rt >
void _test_natural_neighbors_2_with_outputfunctor(Rt T) // intentional copy because T is modified
{
  std::cout << "Testing with OutputFunctor..." << std::endl;

  typedef typename Rt::Bare_point                                Bare_point;
  typedef typename Rt::Weighted_point                            Weighted_point;
  typedef typename Rt::Geom_traits                               Gt;
  typedef typename Gt::FT                                        Coord_type;

  typedef typename Rt::Vertex_handle                             Vertex_handle;
  typedef typename Rt::Face_handle                               Face_handle;
  typedef std::pair<Face_handle, int>                            Edge;

  typedef std::vector<std::pair<Vertex_handle, Coord_type> >     Vertex_coordinate_vector;
  typedef typename Vertex_coordinate_vector::const_iterator      VCV_cit;
  typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> >  Identity_output_functor;

  typedef std::vector<std::pair<Weighted_point, Coord_type> >                   Point_coordinate_vector;
  typedef typename Point_coordinate_vector::const_iterator                      PCV_cit;
  typedef CGAL::Interpolation::internal::Extract_point_in_pair<Rt, Coord_type>  Point_output_functor;

  Vertex_coordinate_vector vh_coords;
  Identity_output_functor vh_fct;

  Point_coordinate_vector pt_coords;
  Point_output_functor pt_fct;

  // test with 0 weight:
  Weighted_point wp(Bare_point(0,0), 0.);
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, Coord_type, bool> pt_coordinate_result =
    CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);

  Coord_type norm = pt_coordinate_result.second;
  assert(norm == Coord_type(1));

  PCV_cit ci = pt_coords.begin();
  for(; ci!= pt_coords.end(); ci++)
    assert(ci->second == Coord_type(0.25));

  bool is_equal = test_barycenter<Rt>(pt_coords.begin(), pt_coords.end(), norm, wp);
  assert(is_equal);
  pt_coords.clear();

  // with an absurdly low weight to have a vertex that is hidden on insertion
  wp = Weighted_point(Bare_point(0,0), -1000.);
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coords.empty());
  assert(pt_coordinate_result.third);
  assert(pt_coordinate_result.second == 0);

  // test with hidden_vertices:
  wp = Weighted_point(Bare_point(0,0), 4.);
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);
  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(pt_coords.begin(), pt_coords.end(), norm, wp);
  assert(is_equal);
  pt_coords.clear();

  // add the middle point of the grid
  T.insert(Weighted_point(Bare_point(0,0), 0.));

  // point on a vertex;
  wp = Weighted_point(Bare_point(1,0), 0.);
  CGAL::Triple<std::back_insert_iterator<Vertex_coordinate_vector>, Coord_type, bool> vh_coordinate_result =
    CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(vh_coords), vh_fct);
  assert(vh_coordinate_result.third);

  norm = vh_coordinate_result.second;
  assert(norm == Coord_type(1));
  is_equal = test_norm(vh_coords.begin(), vh_coords.end(), norm);
  assert(is_equal);

  VCV_cit vh_ci = vh_coords.begin();
  assert(vh_ci->first->point() == wp);
  assert(vh_ci->second == Coord_type(1));
  ++vh_ci;
  assert(vh_ci == vh_coords.end());
  vh_coords.clear();

  // point on the vertex but creating a hole:
  wp = Weighted_point(Bare_point(1,0), 2.);
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(pt_coords.begin(), pt_coords.end(), norm, wp);
  assert(is_equal);
  pt_coords.clear();

  // point on an edge:
  wp = Weighted_point(Bare_point(0,0.5), 3.);
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(pt_coords.begin(), pt_coords.end(), norm, wp);
  assert(is_equal);
  pt_coords.clear();

  // a vertex v in Reg(P\v->point()):
  typename Rt::Vertex_iterator vit = T.finite_vertices_end();
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, --vit, std::back_inserter(pt_coords), pt_fct);
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(pt_coords.begin(), pt_coords.end(), norm, vit->point());
  assert(is_equal);
  pt_coords.clear();

  ci = pt_coords.begin();
  while(ci != pt_coords.end())
  {
    assert(ci->first != (--vit)->point());
    assert(ci->second != Coord_type(1));
    ++ci;
  }
  pt_coords.clear();

  // with the conflict hole
  std::list<Edge> hole_bd;
  std::list<Vertex_handle> hidden_vertices;
  T.get_boundary_of_conflicts_and_hidden_vertices(wp, std::back_inserter(hole_bd), std::back_inserter(hidden_vertices));
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct,
                                                              hole_bd.begin(), hole_bd.end(),
                                                              hidden_vertices.begin(), hidden_vertices.end());
  assert(pt_coordinate_result.third);

  norm = pt_coordinate_result.second;
  is_equal = test_norm(pt_coords.begin(), pt_coords.end(), norm);
  assert(is_equal);

  is_equal = test_barycenter<Rt>(pt_coords.begin(), pt_coords.end(), norm, wp);
  assert(is_equal);
  pt_coords.clear();

  // outside convex hull:
  wp = Weighted_point(Bare_point(3,0.5), 3.);
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(!pt_coordinate_result.third);

  // same but with the face hint API
  vh_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(vh_coords), vh_fct,
                                                              T.finite_faces_begin());
  assert(!vh_coordinate_result.third);

  // on a convex hull edge:
  wp = Weighted_point(Bare_point(2,1), 3.);
  pt_coordinate_result = CGAL::regular_neighbor_coordinates_2(T, wp, std::back_inserter(pt_coords), pt_fct);
  assert(!pt_coordinate_result.third);
}

template <class Rt>
void _test_regular_neighbors_2(const Rt &)
{
  typedef typename Rt::Bare_point                             Bare_point;
  typedef typename Rt::Weighted_point                         Weighted_point;

  Rt T;

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

  T.insert(Weighted_point(p1_2, 0));
  T.insert(Weighted_point(p2_2, 0));
  T.insert(Weighted_point(p3_2, 0));
  T.insert(Weighted_point(p4_2, 0));
  T.insert(Weighted_point(p1, 0));
  T.insert(Weighted_point(p2, 0));
  T.insert(Weighted_point(p3, 0));
  T.insert(Weighted_point(p4, 0));
  T.insert(Weighted_point(p12, 0));
  T.insert(Weighted_point(p23, 0));
  T.insert(Weighted_point(p34, 0));
  T.insert(Weighted_point(p41, 0));

  _test_natural_neighbors_2_without_outputfunctor(T);
  _test_natural_neighbors_2_with_outputfunctor(T);
}
