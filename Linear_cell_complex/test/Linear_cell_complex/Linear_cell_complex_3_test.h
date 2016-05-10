// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LCC_3_TEST_H
#define CGAL_LCC_3_TEST_H

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include "Linear_cell_complex_2_test.h"
#include <fstream>

template<typename LCC>
bool check_number_of_cells_3(LCC& lcc, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbvol,
                             unsigned int nbcc)
{
  if ( !lcc.is_valid() )
  {
    std::cout<<"ERROR: the lcc is not valid."<<std::endl;
    assert(false);
    return false;
  }

  std::vector<unsigned int> nbc;
  nbc=lcc.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbvol!=nbc[3] ||
      nbcc!=nbc[4])
  {
    std::cout<<"ERROR: the number of cells is not correct. We must have "
             <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbvol<<", "<<nbcc
             <<") and we have"<<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "
             <<nbc[3]<<", "<<nbc[4]<<")."
             <<std::endl;
    assert(false);
    return false;
  }

  if ( nbv!=lcc.number_of_vertex_attributes() )
  {
    std::cout<<"ERROR: the number of vertices ("<<nbv<<") is different than "
             <<"the number of vertex attributes ("
             <<lcc.number_of_vertex_attributes()<<")"<<std::endl;

    assert(false);
    return false;
  }

  trace_test_end();

  return true;
}

template<typename LCC>
typename LCC::Dart_handle make_loop(LCC& lcc, const typename LCC::Point& p1)
{
  typename LCC::Dart_handle dh1 = lcc.create_dart(p1);
  lcc.template sew<1>(dh1, dh1);
  return dh1;
}

template<typename LCC>
typename LCC::Dart_handle make_face_two_edges(LCC& lcc,
                                              const typename LCC::Point& p1,
                                              const typename LCC::Point& p2)
{
  typename LCC::Dart_handle dh1 = lcc.create_dart(p1);
  lcc.template sew<1>(dh1, lcc.create_dart(p2));
  lcc.template sew<0>(dh1, lcc.beta(dh1, 1));
  return dh1;
}

template<typename LCC>
bool test_LCC_3()
{
  LCC lcc;

  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;

  // Construction operations
  trace_test_begin();
  Dart_handle dh1=lcc.make_segment(Point(0,0,0),Point(1,0,0));
  Dart_handle dh2=lcc.make_segment(Point(2,0,0),Point(2,1,0));
  Dart_handle dh3=lcc.make_segment(Point(2,2,0),Point(3,1,0));
  if ( !check_number_of_cells_3(lcc, 6, 3, 6, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<0>(dh2,dh1);
  lcc.template sew<1>(dh2,dh3);
  if ( !check_number_of_cells_3(lcc, 4, 3, 4, 1, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh5=lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  Dart_handle dh6=lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  if ( !check_number_of_cells_3(lcc, 10, 9, 6, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_3(lcc, 8, 8, 6, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 9, 9, 6, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_3(lcc, 10, 12, 8, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0,3));
  if ( !check_number_of_cells_3(lcc, 11, 13, 8, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5,3));
  if ( !check_number_of_cells_3(lcc, 12, 16, 10, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2,3));
  if ( !check_number_of_cells_3(lcc, 13, 17, 10, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh12 = lcc.make_tetrahedron(Point(-1, 0, 0),Point(0, 2, 0),
                                          Point(1, 0, 0),Point(1, 1, 2));
  Dart_handle dh13 = lcc.make_tetrahedron(Point(0, 2, -1),Point(-1, 0, -1),
                                          Point(1, 0, -1),Point(1, 1, -3));
  if ( !check_number_of_cells_3(lcc, 21, 29, 18, 4, 4) )
    return false;

  trace_test_begin();
  lcc.template sew<3>(dh12, dh13);
  if ( !check_number_of_cells_3(lcc, 18, 26, 17, 4, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh14=lcc.template insert_barycenter_in_cell<2>(dh12);
  if ( !check_number_of_cells_3(lcc, 19, 29, 19, 4, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh15=lcc.template insert_barycenter_in_cell<1>(dh14);
  if ( !check_number_of_cells_3(lcc, 20, 30, 19, 4, 3) )
    return false;

  // Removal operations
  trace_test_begin();
  lcc.template remove_cell<0>(dh15);
  if ( !check_number_of_cells_3(lcc, 19, 29, 19, 4, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(lcc.beta(dh14, 2, 1));
  lcc.template remove_cell<1>(lcc.beta(dh14, 0));
  lcc.template remove_cell<1>(dh14);
  if ( !check_number_of_cells_3(lcc, 18, 26, 17, 4, 3) )
    return false;

  trace_test_begin();
  lcc.template unsew<3>(dh12);
  if ( !check_number_of_cells_3(lcc, 21, 29, 18, 4, 4) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<3>(dh13);
  lcc.template remove_cell<3>(dh12);
  if ( !check_number_of_cells_3(lcc, 13, 17, 10, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(dh11);
  if ( !check_number_of_cells_3(lcc, 12, 16, 10, 2, 2) )
    return false;

  trace_test_begin();
  std::vector<Dart_handle> toremove;
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh10).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh10).end();
        it!=itend; ++it )
    toremove.push_back( it );

  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    lcc.template remove_cell<1>(*it);
  toremove.clear();
  if ( !check_number_of_cells_3(lcc, 11, 13, 8, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh9);
  if ( !check_number_of_cells_3(lcc, 10, 12, 8, 2, 2) )
    return false;

  trace_test_begin();
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh8).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh8).end();
        it!=itend; ++it )
    toremove.push_back( it );

  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    lcc.template remove_cell<1>(*it);
  toremove.clear();
  if ( !check_number_of_cells_3(lcc, 9, 9, 6, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh7);
  if ( !check_number_of_cells_3(lcc, 8, 8, 6, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<2>(dh5);
  if ( !check_number_of_cells_3(lcc, 10, 9, 6, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<2>(dh6);
  lcc.template remove_cell<2>(dh5);
  if ( !check_number_of_cells_3(lcc, 4, 3, 4, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template unsew<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 5, 3, 5, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<0>(dh2);
  if ( !check_number_of_cells_3(lcc, 6, 3, 6, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(dh1);
  lcc.template remove_cell<1>(dh2);
  lcc.template remove_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  // Edge contraction
  trace_test_begin();
  dh1 = lcc.create_dart(Point(0,0,0));
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = make_loop(lcc, Point(0,0,0));
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  lcc.template sew<1>(lcc.beta(dh1, 2), lcc.beta(dh1, 2));
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.beta(dh1,0); dh3 = lcc.beta(dh1,1);
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 1, 1) ||
       !lcc.is_face_combinatorial_polygon(dh2, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) ||
       !lcc.is_face_combinatorial_polygon(dh3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  lcc.template sew<3>(dh1, dh2); dh2 = lcc.beta(dh1, 0); dh3 = lcc.beta(dh1, 1);

  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 2, 1) ||
       !lcc.is_face_combinatorial_polygon(dh2, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 2, 1) ||
       !lcc.is_face_combinatorial_polygon(dh3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  lcc.template sew<2>(dh1, dh2);

  dh2 = lcc.beta(dh2, 1);
  dh3 = lcc.beta(dh1, 1);

  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 4, 4, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh2, 1));
  if ( !check_number_of_cells_3(lcc, 3, 3, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh3, 1));
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  lcc.template sew<2>(dh1, dh2);

  dh2 = lcc.beta(dh2, 1);
  dh3 = lcc.beta(dh1, 1);

  lcc.template contract_cell<1>(lcc.beta(dh2, 1));
  if ( !check_number_of_cells_3(lcc, 3, 4, 2, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh3, 1));
  if ( !check_number_of_cells_3(lcc, 2, 3, 2, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 1, 2, 2, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 1, 1, 2, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  lcc.template sew<2>(dh1, dh2);

  dh3 = lcc.make_triangle(Point(5,5,4),Point(7,5,4),Point(6,6,4));
  lcc.template sew<3>(dh1, dh3);

  dh3 = lcc.make_triangle(Point(5,4,4),Point(7,4,4),Point(6,3,4));
  lcc.template sew<3>(dh2, dh3);

  lcc.template sew<2>(lcc.beta(dh1, 3), dh3);

  dh2 = lcc.beta(dh2,1);
  dh3 = lcc.beta(dh1,1);

  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 4, 4, 2, 4, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh2,1));
  if ( !check_number_of_cells_3(lcc, 3, 3, 2, 4, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh3,1));
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  lcc.template sew<2>(dh1, dh2);

  dh3 = lcc.make_triangle(Point(5,5,4),Point(7,5,4),Point(6,6,4));
  lcc.template sew<3>(dh1, dh3);

  dh3 = lcc.make_triangle(Point(5,4,4),Point(7,4,4),Point(6,3,4));
  lcc.template sew<3>(dh2, dh3);

  dh2 = lcc.beta(dh2,1);
  dh3 = lcc.beta(dh1,1);

  lcc.template contract_cell<1>(lcc.beta(dh2,1));
  if ( !check_number_of_cells_3(lcc, 3, 4, 2, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh3,1));
  if ( !check_number_of_cells_3(lcc, 2, 3, 2, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 1, 2, 2, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 1, 1, 2, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.
      make_hexahedron(Point(0,0,0),Point(1,0,0),Point(1,1,0),Point(0,1,0),
                      Point(0,1,1),Point(0,0,1),Point(1,0,1),Point(1,1,1));
  dh2 = lcc.
      make_hexahedron(Point(0,3,0),Point(1,3,0),Point(1,4,0),Point(0,4,0),
                      Point(0,4,1),Point(0,3,1),Point(1,3,1),Point(1,4,1));
  dh2 = lcc.beta(dh2, 2,1,1,2);
  lcc.template sew<3>(dh1,dh2);

  lcc.template contract_cell<1>(lcc.beta(dh1,0));
  if ( !check_number_of_cells_3(lcc, 11, 19, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh1,1,1));
  if ( !check_number_of_cells_3(lcc, 10, 18, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.beta(dh1,1));
  if ( !check_number_of_cells_3(lcc, 9, 17, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 10, 16, 10, 2, 2 ) )
    return false;
  lcc.clear();

  // Face contraction
  trace_test_begin();
  dh1 = lcc.create_dart(Point(0,0,0));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = make_loop<LCC>(lcc, Point(0,0,0));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  lcc.template sew<1>(lcc.beta(dh1,2), lcc.beta(dh1,2));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = make_face_two_edges(lcc, Point(0,0,0), Point(1,0,0));
  dh2 = make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1));
  lcc.template sew<2>(dh1, dh2);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(dh2);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  dh3 = lcc.make_triangle(Point(5,3,3),Point(7,3,3),Point(6,0,3));
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<2>(lcc.beta(dh2,1), dh3);

  lcc.template contract_cell<1>(lcc.beta(dh2,0));
  if ( !check_number_of_cells_3(lcc, 4, 6, 3, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(dh2);
  if ( !check_number_of_cells_3(lcc, 4, 5, 2, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.create_dart(Point(0,0,0));
  lcc.template sew<3>(dh1, lcc.create_dart(Point(1,0,0)));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = make_loop(lcc, Point(0,0,0));
  dh2 = make_loop(lcc, Point(0,0,1));
  lcc.template sew<3>(dh1, dh2);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<3>(dh1, lcc.make_segment(Point(0,0,1),Point(1,0,1)));
  lcc.template sew<3>(lcc.beta(dh1,2),lcc.beta(dh1,3,2));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 1, 1, 2, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  dh2 = lcc.make_segment(Point(0,0,1),Point(1,0,1));
  lcc.template sew<1>(dh2, dh2);
  lcc.template sew<3>(dh1, dh2);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 2, 2, 2) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  lcc.template sew<1>(lcc.beta(dh1,2), lcc.beta(dh1,2));
  dh2 = lcc.make_segment(Point(0,0,1),Point(1,0,1));
  lcc.template sew<1>(dh2, dh2);
  lcc.template sew<1>(lcc.beta(dh2,2), lcc.beta(dh2,2));
  lcc.template sew<3>(dh1, dh2);
  lcc.template sew<3>(lcc.beta(dh1,2), lcc.beta(dh2,2));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 2, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = make_face_two_edges(lcc, Point(0,0,0), Point(1,0,0));
  dh2 = make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1));
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<3>(dh1,
                      make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1)));
  lcc.template sew<3>(dh2,
                      make_face_two_edges(lcc, Point(1,0,1), Point(1,0,2)));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(dh2);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = make_face_two_edges(lcc, Point(0,0,0), Point(1,0,0));
  dh2 = make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1));
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<3>(dh1,
                      make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1)));
  lcc.template sew<3>(dh2,
                      make_face_two_edges(lcc, Point(1,0,1), Point(1,0,2)));
  lcc.template sew<2>(lcc.beta(dh1,3), lcc.beta(dh2,3));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(dh2);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  dh3 = lcc.make_triangle(Point(5,3,3),Point(7,3,3),Point(6,0,3));
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<2>(lcc.beta(dh2,1), dh3);
  lcc.template sew<3>(dh1, lcc.make_triangle(Point(5,5,4),Point(7,5,4),
                                             Point(6,6,4)));
  lcc.template sew<3>(dh2, lcc.make_triangle(Point(5,4,4),Point(7,4,4),
                                             Point(6,3,4)));
  lcc.template sew<3>(dh3, lcc.make_triangle(Point(5,3,4),Point(7,3,4),
                                             Point(6,0,4)));
  lcc.template sew<2>(lcc.beta(dh1,3), lcc.beta(dh2,3));
  lcc.template sew<2>(lcc.beta(dh2,1,3), lcc.beta(dh3,3));
  lcc.template contract_cell<1>(lcc.beta(dh2,0));
  if ( !check_number_of_cells_3(lcc, 4, 6, 3, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(dh2);
  if ( !check_number_of_cells_3(lcc, 4, 5, 2, 2, 1) )
    return false;
  lcc.clear();

  // Volume contraction
  trace_test_begin();
  dh1 = lcc.create_dart(Point(0,0,0));
  lcc.template contract_cell<3>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = make_loop<LCC>(lcc, Point(0,0,0));
  lcc.template contract_cell<3>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template contract_cell<3>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  lcc.template sew<1>(dh1, dh1);
  lcc.template sew<1>(lcc.beta(dh1,2), lcc.beta(dh1,2));
  lcc.template contract_cell<3>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.
    make_hexahedron(Point(0,0,0),Point(1,0,0),Point(1,1,0),Point(0,1,0),
                    Point(0,1,1),Point(0,0,1),Point(1,0,1),Point(1,1,1));
  dh2 = lcc.
    make_hexahedron(Point(0,3,0),Point(1,3,0),Point(1,4,0),Point(0,4,0),
                    Point(0,4,1),Point(0,3,1),Point(1,3,1),Point(1,4,1));
  dh2 = lcc.beta(dh2, 2,1,1,2);
  lcc.template sew<3>(dh1,dh2);

  lcc.template contract_cell<1>(lcc.beta(dh2,2,1));
  lcc.template contract_cell<1>(lcc.beta(dh2,2,0));
  lcc.template contract_cell<1>(lcc.beta(dh2,1,1,2,0));
  lcc.template contract_cell<1>(lcc.beta(dh2,1,1,2,1));

  if ( !check_number_of_cells_3(lcc, 8, 16, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,0,2));
  if ( !check_number_of_cells_3(lcc, 8, 15, 10, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,1,1,2));
  if ( !check_number_of_cells_3(lcc, 8, 14, 9, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,1,2));
  if ( !check_number_of_cells_3(lcc, 8, 13, 8, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,2));
  if ( !check_number_of_cells_3(lcc, 8, 12, 7, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<3>(dh2);
  if ( !check_number_of_cells_3(lcc, 8, 12, 6, 1, 1 ) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.
    make_hexahedron(Point(0,0,0),Point(1,0,0),Point(1,1,0),Point(0,1,0),
                    Point(0,1,1),Point(0,0,1),Point(1,0,1),Point(1,1,1));
  dh2 = lcc.
    make_hexahedron(Point(0,3,0),Point(1,3,0),Point(1,4,0),Point(0,4,0),
                    Point(0,4,1),Point(0,3,1),Point(1,3,1),Point(1,4,1));
  dh3 = lcc.
    make_hexahedron(Point(0,6,0),Point(1,6,0),Point(1,7,0),Point(0,7,0),
                    Point(0,7,1),Point(0,6,1),Point(1,6,1),Point(1,7,1));
  dh3 = lcc.beta(dh3, 2,1,1,2);
  lcc.template sew<3>(dh2,dh3);
  dh2 = lcc.beta(dh2, 2,1,1,2);
  lcc.template sew<3>(dh1,dh2);

  lcc.template contract_cell<1>(lcc.beta(dh2,2,1));
  lcc.template contract_cell<1>(lcc.beta(dh2,2,0));
  lcc.template contract_cell<1>(lcc.beta(dh2,1,1,2,0));
  lcc.template contract_cell<1>(lcc.beta(dh2,1,1,2,1));

  if ( !check_number_of_cells_3(lcc, 12, 24, 16, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,0,2));
  if ( !check_number_of_cells_3(lcc, 12, 23, 15, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,1,1,2));
  if ( !check_number_of_cells_3(lcc, 12, 22, 14, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,1,2));
  if ( !check_number_of_cells_3(lcc, 12, 21, 13, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.beta(dh2,2));
  if ( !check_number_of_cells_3(lcc, 12, 20, 12, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<3>(dh2);
  if ( !check_number_of_cells_3(lcc, 12, 20, 11, 2, 1 ) )
    return false;
  lcc.clear();

  // Construction from Polyhedron_3
  {
    trace_test_begin();
    CGAL::Polyhedron_3<typename LCC::Traits> P;
    std::ifstream in("data/armadillo.off");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/armadillo.off'"<<std::endl;
      return false;
    }
    in >> P;
    CGAL::import_from_polyhedron_3<LCC>(lcc,P);
    if ( !check_number_of_cells_3(lcc, 26002, 78000, 52000, 1, 1) )
      return false;
    lcc.clear();
  }

  // Construction from Triangulation_3
  {
    trace_test_begin();
    CGAL::Triangulation_3<typename LCC::Traits> T;
    std::ifstream in("data/points.txt");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/points.txt'"<<std::endl;
      return false;
    }
    T.insert ( std::istream_iterator < Point >(in),
               std::istream_iterator < Point >() );
    CGAL::import_from_triangulation_3<LCC>(lcc,T);
    // Pb: the triangulation_3 is not the same on different machines ?
    // if ( !check_number_of_cells_3(lcc, 795, 4156, 6722, 3361, 1) )
    if ( !lcc.is_valid() )
      return false;
    lcc.clear();
    trace_test_end();
  }

  return true;
}

#endif // CGAL_LCC_3_TEST_H
