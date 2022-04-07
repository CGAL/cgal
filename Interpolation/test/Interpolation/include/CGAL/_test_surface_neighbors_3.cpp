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
// Author(s)     : Julia Floetotto

#include <CGAL/Voronoi_intersection_2_traits_3.h>
#include <CGAL/surface_neighbor_coordinates_3.h>
#include <CGAL/surface_neighbors_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/function_objects.h>
#include <CGAL/Origin.h>

#include <iostream>
#include <cassert>
#include <fstream>
#include <utility>


template < class ForwardIterator >
bool test_norm(ForwardIterator first, ForwardIterator beyond,
               typename std::iterator_traits<ForwardIterator>::value_type::second_type norm)
{
  typename
  std::iterator_traits<ForwardIterator>::value_type::second_type sum(0);
  for(; first !=beyond; first++)
    sum += first->second;

  return norm == sum;
}

template < class ForwardIterator >
bool test_barycenter(ForwardIterator first, ForwardIterator beyond,
                     typename std::iterator_traits<ForwardIterator>::value_type::second_type norm ,
                     const typename std::iterator_traits<ForwardIterator>::value_type::first_type& p,
                     const typename std::iterator_traits<ForwardIterator>::value_type::second_type& tolerance)
{
  typedef typename
  std::iterator_traits<ForwardIterator>::value_type::first_type  Point;
  Point b(0,0,0);
  ForwardIterator it=first;
  for(; it!=beyond; it++)
    b = b + (it->second/norm) * (it->first - CGAL::ORIGIN);
  return CGAL::squared_distance(p,b) <= tolerance;
}

template < class ForwardIteratorCoord, class ForwardIteratorPoint, class Kernel >
bool compare_neighbors(ForwardIteratorCoord first_coord,
                       ForwardIteratorCoord beyond_coord,
                       ForwardIteratorPoint first_point,
                       ForwardIteratorPoint beyond_point,
                       Kernel )
{
  typedef std::set<typename std::iterator_traits<ForwardIteratorPoint>::value_type,
      typename Kernel::Less_xyz_3 > Point_set;

  Point_set neighbors, coord_neighbors, diff;

  neighbors.insert(first_point, beyond_point);

  for(; first_coord !=beyond_coord; first_coord++)
    coord_neighbors.insert(first_coord->first);

  std::set_difference(neighbors.begin(),
                      neighbors.end(),
                      coord_neighbors.begin(),
                      coord_neighbors.end(),
                      std::inserter(diff, diff.begin()),
                      typename Kernel::Less_xyz_3());

  if(!diff.empty())
  {
    std::cout << "Compare neighbors -- Diff: " << std::endl;
    for(typename Point_set::const_iterator it = diff.begin(); it != diff.end(); ++it)
      std::cout << " point: " << *it;
    std::cout << std::endl;
  }

  return diff.empty();
}

template < class Tr, class OutputIterator>
OutputIterator
test_neighbors(const Tr& T, const typename Tr::Geom_traits::Point_3& p,
               const typename Tr::Geom_traits::Vector_3& n,
               const int& version,
               OutputIterator out)
{
  typedef CGAL::Voronoi_intersection_2_traits_3<typename Tr::Geom_traits>   I_traits;

  //the result type of the certified version:
  typedef std::pair<  OutputIterator, bool >  NeighborIt_bool_pair;

  typename Tr::Cell_handle start;
  //test different function calls
  switch(version)
  {
    case 0:{
      //certified call with Kernel:
      NeighborIt_bool_pair
          result_pair = CGAL::surface_neighbors_certified_3(T.points_begin(),
                                                            T.points_end(), p, n,
                                                            out,
                                                            T.geom_traits());
      assert(result_pair.second);
      out = result_pair.first; break;}
    case 1: {
      //certified call with instantiated traits::
      NeighborIt_bool_pair
          result_pair = CGAL::surface_neighbors_certified_3(T.points_begin(),
                                                            T.points_end(), p,
                                                            out, I_traits(p,n));
      assert(result_pair.second);
      out =result_pair.first; break;}
      //both versions with locate:
    case 2:{
      start = T.locate(p);
      //certified call with Kernel and locate:
      out =CGAL::surface_neighbors_3(T, p,n,out, start);
      break;}
    case 3: {
      start = T.locate(p);
      //with instantiated traits and locate:
      out =CGAL::surface_neighbors_3(T,p,out, I_traits(p,n),start);
      break;}
      //taking all points:
    case 4: {
      //with instantiated traits and locate:
      out =
          CGAL::surface_neighbors_3(T,p,out,I_traits(p,n));
      break;}
    case 5: {
      //certified call with Kernel and locate:
      out =
          CGAL::surface_neighbors_3(T, p,n,out);
      break;}
      //the last two with certification:
    case 6: {
      out =
          CGAL::surface_neighbors_3(T.points_begin(),
                                    T.points_end(),
                                    p, out,I_traits(p,n));
      break;
    }
    case 7: {
      out =
          CGAL::surface_neighbors_3(T.points_begin(),
                                    T.points_end(),
                                    p,n,out,T.geom_traits());
      break;
    }
    default:
      std::cout << "Switch function calls: Nothing is tested. " <<
                   std::endl;
  }
  return out;
}


template < class Tr, class OutputIterator>
std::pair< OutputIterator,  typename Tr::Geom_traits::FT>
test_coords(const Tr& T,
            const typename Tr::Geom_traits::Point_3& p,
            const typename Tr::Geom_traits::Vector_3& n,
            const int& version, OutputIterator out)
{
  typedef CGAL::Voronoi_intersection_2_traits_3<typename Tr::Geom_traits>     I_traits;

  //coordinate computation result types
  typedef CGAL::Triple< OutputIterator, typename Tr::Geom_traits::FT, bool >  Result_triple;
  //the result type of the certified version:
  typedef CGAL::Quadruple< OutputIterator, typename Tr::Geom_traits::FT, bool, bool > Result_quadruple;

  typename Tr::Cell_handle start;
  typename Tr::Geom_traits::FT  norm = 1; // 1 for that default doesn't trigger an assert
  //test different function calls
  switch(version){
    case 0:{
      Result_triple result
          = CGAL::surface_neighbor_coordinates_3(T, p,n,out);
      assert(result.third);
      norm =  result.second;
      break;}
    case 1: {
      Result_triple result  =
          CGAL::surface_neighbor_coordinates_3(T, p,out,I_traits(p,n));
      assert(result.third);
      norm =  result.second; break;}
      //both versions with locate:
    case 2:{
      start = T.locate(p);
      Result_triple result  = CGAL::surface_neighbor_coordinates_3(T, p, n, out, start);
      assert(result.third);
      norm =  result.second; break;}
    case 3: {
      start = T.locate(p);
      Result_triple result  =
          CGAL::surface_neighbor_coordinates_3(T, p, out, I_traits(p,n), start);
      assert(result.third);
      norm =  result.second; break;}
      //taking all points:
    case 4: {
      Result_triple result
          = CGAL::surface_neighbor_coordinates_3(T.points_begin(),
                                                 T.points_end(), p, n,
                                                 out,
                                                 T.geom_traits());
      assert(result.third);
      norm =  result.second; break;}
    case 5: {
      Result_triple result
          = CGAL::surface_neighbor_coordinates_3(T.points_begin(),
                                                 T.points_end(), p,
                                                 out ,I_traits(p,n));
      assert(result.third);
      norm =  result.second; break;}
      //the last two with certification:
    case 6: {
      Result_quadruple
          result = CGAL::surface_neighbor_coordinates_certified_3
                   (T.points_begin(), T.points_end(),p,n,
                    out, T.geom_traits());
      assert(result.third && result.fourth);
      norm =  result.second; break;
    }
    case 7: {
      Result_quadruple
          result = CGAL::surface_neighbor_coordinates_certified_3
                   (T.points_begin(), T.points_end(),p, out ,I_traits(p,n));
      assert(result.third && result.fourth);
      norm =  result.second;
      break;
    }
    default:
      std::cout << "Switch function calls: Nothing is tested. " <<
                   std::endl;
  }
  assert(norm > 0);

  return std::make_pair(out, norm);
}

//call for test functions:
template < class Tr>
void
test_coords_and_neighbors(const Tr& T, const  typename
                          Tr::Geom_traits::Point_3& p,
                          const typename Tr::Geom_traits::Vector_3& n,
                          const typename  Tr::Geom_traits::FT& tolerance, const int& version)
{
  CGAL::Set_ieee_double_precision pfr;

  typedef std::pair<typename Tr::Geom_traits::Point_3,
      typename Tr::Geom_traits::FT >    Point_coord_pair;

  std::vector<Point_coord_pair>  coords;
  typename Tr::Geom_traits::FT   norm;
  norm = test_coords(T, p, n,version, std::back_inserter(coords)).second;
  assert(test_norm( coords.begin(), coords.end(), norm));
  assert(test_barycenter(coords.begin(), coords.end(), norm, p, tolerance));
  //All function testing surface neighbors are
  // grouped together:
  std::vector< typename Tr::Geom_traits::Point_3 >  neighbors;
  test_neighbors(T, p, n,version, std::back_inserter(neighbors));
  assert(compare_neighbors(coords.begin(),
                           coords.end(),neighbors.begin(),
                           neighbors.end(), T.geom_traits()));
  //done
}

template <class Tr>
void
_test_surface_neighbors_3_sphere( const Tr & )
{
  Tr T;

  int n=200, m=20;
  double r = 3;

  typedef typename Tr::Geom_traits                 Gt;
  typedef typename Gt::Point_3                     Point;

  std::vector<Point> points;
  points.reserve(n+m);

  // Create n+m-4 points on a sphere of radius 2
  CGAL::Random_points_on_sphere_3<Point> g(r);
  std::copy_n( g, n+m, std::back_inserter(points));

  for(int i=0; i<n ; i++)
    T.insert(points[i]);

  //test with different calls:
  int k=0;
  for(int i=n;i<n+m;i++)
  {
    test_coords_and_neighbors(T, points[i],
                              typename Gt::Vector_3(points[i]-CGAL::ORIGIN),
                              typename Gt::FT(0.5),++k % 8);
  }
}

////cube case in the interior of a face
template <class Tr, class Transformation>
void
_test_surface_neighbors_3_cube(const Tr &, const Transformation&
                               transform, const int n = 75,
                               typename Tr::Geom_traits::FT tolerance = typename Tr::Geom_traits::FT(1e-29),
                               bool grid=true)
{
  Tr T;

  int m=10;
  double r = 3;

  typedef typename Tr::Geom_traits                Gt;
  typedef typename Gt::FT                         Coord_type;
  typedef typename Gt::Point_3                    Point;
  typedef typename Gt::Point_2                    Point_2;
  typedef typename Gt::Vector_3                   Vector;

  //data points: generate random points in a square of length r
  std::vector<Point_2> points_2_data;
  points_2_data.reserve(n);

  if(grid)
  {
    CGAL::points_on_square_grid_2(r, n, std::back_inserter(points_2_data),
                                  CGAL::Creator_uniform_2<Coord_type,Point_2>());
  }
  else
  {
    CGAL::Random_points_in_square_2<Point_2> g(r);
    std::copy_n(g, n, std::back_inserter(points_2_data));
  }
  for(int i=0; i < n; i++)
  {
    T.insert(transform(Point(points_2_data[i].x(),points_2_data[i].y(), -r)));
    T.insert(transform(Point(points_2_data[i].x(),points_2_data[i].y(), r)));
    T.insert(transform(Point(-r, points_2_data[i].x(),points_2_data[i].y())));
    T.insert(transform(Point(r, points_2_data[i].x(), points_2_data[i].y())));
    T.insert(transform(Point(points_2_data[i].x(), -r, points_2_data[i].y())));
    T.insert(transform(Point(points_2_data[i].x(), r, points_2_data[i].y())));
  }

  //test_points: generate random points in a square of length r
  std::vector<Point_2> points_2_test;
  points_2_test.reserve(m);
  CGAL::Random_points_in_square_2<Point_2> g2(r-1.0);
  std::copy_n(g2, m, std::back_inserter(points_2_test));

  int k=0;
  for(int i=0;i<m;i++)
  {
    //test point on z=r plane:
    test_coords_and_neighbors(T,transform(Point(points_2_test[i].x(),
                                                points_2_test[i].y(), r)),
                              transform(Vector(0,0,1)),tolerance, ++k % 8);
    //test point on x=-r plane:
    test_coords_and_neighbors(T,transform(Point(-r, points_2_test[i].x(),
                                                points_2_test[i].y())),
                              transform(Vector(-1,0,0)),tolerance, ++k % 8 );
    //test point on x=r plane:
    test_coords_and_neighbors(T,transform(Point(r, points_2_test[i].x(),
                                                points_2_test[i].y())),
                              transform(Vector(1,0,0)),tolerance,++k % 8 );
    //test point on y=-r plane:
    test_coords_and_neighbors(T,transform(Point(points_2_test[i].x(),
                                                -r,points_2_test[i].y())),
                              transform(Vector(0,-1,0)),tolerance,++k % 8 );
    //test point on y=r plane:
    test_coords_and_neighbors(T,transform(Point(points_2_test[i].x(),
                                                r,points_2_test[i].y())),
                              transform(Vector(0,1,0)),tolerance,++k % 8);
  }

  //test a sample point:
  //with Delaunay triangulation filering:
  test_coords_and_neighbors(T,transform(Point(points_2_data[n/2].x(),
                                        points_2_data[n/2].y(), r)),
                            transform(Vector(0,0,1)),Coord_type(0),0);

  //considering all points:
  test_coords_and_neighbors(T,transform(Point(points_2_data[n/2].x(),
                                        points_2_data[n/2].y(), r)),
                            transform(Vector(0,0,1)),Coord_type(0),4);
}
