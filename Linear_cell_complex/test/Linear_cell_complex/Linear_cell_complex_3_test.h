// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_LCC_3_TEST_H
#define CGAL_LCC_3_TEST_H

#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3_to_lcc.h>
#include <CGAL/Triangulation_3_to_lcc.h>
#include "Linear_cell_complex_2_test.h"
#include <fstream>
#include <typeinfo>
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

template<typename Map>
void create_attributes_3(Map& map)
{
  create_attributes_2(map);
  CreateAttributes<Map, 3>::run(map);
}

template<typename LCC>
typename LCC::Dart_descriptor make_loop(LCC& lcc, const typename LCC::Point& p1)
{
  typename LCC::Dart_descriptor dh1 = lcc.make_half_edge();
  lcc.set_vertex_attribute(dh1, lcc.create_vertex_attribute(p1));
  lcc.template sew<1>(dh1, lcc.other_orientation(dh1));
  return dh1;
}

template<typename LCC>
typename LCC::Dart_descriptor make_face_two_edges(LCC& lcc,
                                              const typename LCC::Point& p1,
                                              const typename LCC::Point& p2)
{
  typename LCC::Dart_descriptor dh1 = lcc.make_combinatorial_polygon(2);
  lcc.set_vertex_attribute(dh1, lcc.create_vertex_attribute(p1));
  lcc.set_vertex_attribute(lcc.next(dh1), lcc.create_vertex_attribute(p2));
  return dh1;
}

template<typename LCC>
bool test_LCC_3()
{
  LCC lcc;

  typedef typename LCC::Dart_descriptor Dart_descriptor;
  typedef typename LCC::Point Point;

  // Construction operations
  trace_test_begin();
  Dart_descriptor dh1=lcc.make_segment(Point(0,0,0),Point(1,0,0), true);
  Dart_descriptor dh2=lcc.make_segment(Point(2,0,0),Point(2,1,0), true);
  Dart_descriptor dh3=lcc.make_segment(Point(2,2,0),Point(3,1,0), true);
  create_attributes_3(lcc);
  if ( !check_number_of_cells_3(lcc, 6, 3, 6, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<1>(dh1,dh2);
  lcc.template sew<1>(lcc.other_orientation(dh2),dh3);
  if ( !check_number_of_cells_3(lcc, 4, 3, 4, 1, 1) )
    return false;

  trace_test_begin();
  Dart_descriptor dh5=lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  Dart_descriptor dh6=lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  if ( !check_number_of_cells_3(lcc, 10, 9, 6, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_3(lcc, 8, 8, 6, 2, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 9, 9, 6, 2, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_3(lcc, 10, 12, 8, 2, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0,3));
  if ( !check_number_of_cells_3(lcc, 11, 13, 8, 2, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5,3));
  if ( !check_number_of_cells_3(lcc, 12, 16, 10, 2, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2,3));
  if ( !check_number_of_cells_3(lcc, 13, 17, 10, 2, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh12 = lcc.make_tetrahedron(Point(-1, 0, 0),Point(0, 2, 0),
                                          Point(1, 0, 0),Point(1, 1, 2));
  Dart_descriptor dh13 = lcc.make_tetrahedron(Point(0, 2, -1),Point(-1, 0, -1),
                                          Point(1, 0, -1),Point(1, 1, -3));
  if ( !check_number_of_cells_3(lcc, 21, 29, 18, 4, 4) )
    return false;

  trace_test_begin();
  lcc.template sew<3>(dh12, dh13);
  if ( !check_number_of_cells_3(lcc, 18, 26, 17, 4, 3) )
    return false;

  trace_test_begin();
  Dart_descriptor dh14=lcc.template insert_barycenter_in_cell<2>(dh12);
  if ( !check_number_of_cells_3(lcc, 19, 29, 19, 4, 3) )
    return false;

  trace_test_begin();
  Dart_descriptor dh15=lcc.template insert_barycenter_in_cell<1>(dh14);
  if ( !check_number_of_cells_3(lcc, 20, 30, 19, 4, 3) )
    return false;

  // Removal operations
  trace_test_begin();
  lcc.template remove_cell<0>(dh15);
  if ( !check_number_of_cells_3(lcc, 19, 29, 19, 4, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(lcc.next(lcc.template opposite<2>(dh14)));
  lcc.template remove_cell<1>(lcc.previous(dh14));
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
  std::vector<Dart_descriptor> toremove;
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh10).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh10).end();
        it!=itend; ++it )
    toremove.push_back( it );

  for ( typename std::vector<Dart_descriptor>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    if (lcc.is_dart_used(*it)) // For GMap because we have 2 dart per edge incident to the vertex
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

  for ( typename std::vector<Dart_descriptor>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    if (lcc.is_dart_used(*it)) // For GMap because we have 2 dart per edge incident to the vertex
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
  lcc.template unsew<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 5, 3, 5, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<1>(lcc.other_orientation(dh2));
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
  lcc.template sew<1>(dh1, lcc.other_orientation(dh1));
  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.previous(dh1); dh3 = lcc.next(dh1);
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
  lcc.template sew<3>(dh1, dh2);
  create_attributes_3(lcc);
  dh2 = lcc.previous(dh1); dh3 = lcc.next(dh1);

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
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);

  dh2 = lcc.next(dh2);
  dh3 = lcc.next(dh1);

  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 4, 4, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh2));
  if ( !check_number_of_cells_3(lcc, 3, 3, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh3));
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);

  dh2 = lcc.next(dh2);
  dh3 = lcc.next(dh1);

  lcc.template contract_cell<1>(lcc.next(dh2));
  if ( !check_number_of_cells_3(lcc, 3, 4, 2, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh3));
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
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);

  dh3 = lcc.make_triangle(Point(5,5,4),Point(7,5,4),Point(6,6,4));
  lcc.template sew<3>(dh1, dh3);

  dh3 = lcc.make_triangle(Point(5,4,4),Point(7,4,4),Point(6,3,4));
  lcc.template sew<3>(dh2, dh3);

  lcc.template sew<2>(lcc.template opposite<3>(dh1), dh3);

  dh2 = lcc.next(dh2);
  dh3 = lcc.next(dh1);

  lcc.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 4, 4, 2, 4, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh2));
  if ( !check_number_of_cells_3(lcc, 3, 3, 2, 4, 2) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 2, 2, 1, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh3));
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  dh2 = lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);

  dh3 = lcc.make_triangle(Point(5,5,4),Point(7,5,4),Point(6,6,4));
  lcc.template sew<3>(dh1, dh3);

  dh3 = lcc.make_triangle(Point(5,4,4),Point(7,4,4),Point(6,3,4));
  lcc.template sew<3>(dh2, dh3);

  dh2 = lcc.next(dh2);
  dh3 = lcc.next(dh1);

  lcc.template contract_cell<1>(lcc.next(dh2));
  if ( !check_number_of_cells_3(lcc, 3, 4, 2, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh3));
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
  dh1 = lcc.make_tetrahedron(Point(9, 9, 0),Point(9, 0, 9),
                             Point(0, 9, 9),Point(0, 0, 0));
  create_attributes_3(lcc);
  typename LCC::Vector v=CGAL::compute_normal_of_cell_0(lcc, dh1);
  if (v!=typename LCC::Vector(-9,-9,9))
  {
    assert(false);
    return false;
  }
  trace_test_end();
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.
      make_hexahedron(Point(0,0,0),Point(1,0,0),Point(1,2,0),Point(0,2,0),
                      Point(0,3,4),Point(0,0,4),Point(6,0,4),Point(6,3,4));
  create_attributes_3(lcc);
  v=CGAL::compute_normal_of_cell_2(lcc, lcc.template
                                   opposite<2>(lcc.previous(dh1)));
  if (v!=typename LCC::Vector(0,0,1))
  {
    assert(false);
    return false;
  }

  if (lcc.template barycenter<1>(dh1)!=typename LCC::Point(0, 0, 2))
  {
    assert(false);
    return false;
  }

  if (lcc.template barycenter<2>(dh1)!=typename LCC::Point(1.75, 0, 2))
  {
    assert(false);
    return false;
  }

  if (lcc.template barycenter<3>(dh1)!=typename LCC::Point(1.75, 1.25, 2))
  {
    assert(false);
    return false;
  }
  trace_test_end();

  trace_test_begin();
  dh2 = lcc.
      make_hexahedron(Point(0,3,0),Point(1,3,0),Point(1,4,0),Point(0,4,0),
                      Point(0,4,1),Point(0,3,1),Point(1,3,1),Point(1,4,1));
  dh2 = lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh2))));
  create_attributes_3(lcc);
  lcc.template sew<3>(dh1,dh2);

  lcc.template contract_cell<1>(lcc.previous(dh1));
  if ( !check_number_of_cells_3(lcc, 11, 19, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(lcc.next(dh1)));
  if ( !check_number_of_cells_3(lcc, 10, 18, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<1>(lcc.next(dh1));
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
  create_attributes_3(lcc);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0), true);
  create_attributes_3(lcc);
  lcc.template sew<1>(lcc.template opposite<2>(dh1),
                      lcc.other_orientation(lcc.template opposite<2>(dh1)));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0),true);
  create_attributes_3(lcc);
  lcc.template sew<1>(dh1, lcc.other_orientation(dh1));
  lcc.template sew<1>(lcc.template opposite<2>(dh1),
                      lcc.other_orientation(lcc.template opposite<2>(dh1)));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = make_face_two_edges(lcc, Point(0,0,0), Point(1,0,0));
  dh2 = make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1));
  create_attributes_3(lcc);
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
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<2>(lcc.next(dh2), dh3);

  lcc.template contract_cell<1>(lcc.previous(dh2));
  if ( !check_number_of_cells_3(lcc, 4, 6, 3, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(dh2);
  if ( !check_number_of_cells_3(lcc, 4, 5, 2, 1, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.create_dart(Point(0,0,0));
  create_attributes_3(lcc);
  lcc.template sew<3>(dh1, lcc.create_dart(Point(1,0,0)));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = make_loop(lcc, Point(0,0,0));
  dh2 = make_loop(lcc, Point(0,0,1));
  create_attributes_3(lcc);
  lcc.template sew<3>(dh1, dh2);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0), true);
  create_attributes_3(lcc);
  lcc.template sew<3>(dh1, lcc.make_segment(Point(0,0,1),Point(1,0,1), true));
  lcc.template sew<3>(lcc.template opposite<2>(dh1),
                      lcc.template opposite<2>(lcc.template opposite<3>(dh1)));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 1, 1, 2, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0), true);
  create_attributes_3(lcc);
  lcc.template sew<1>(dh1, lcc.other_orientation(dh1));
  dh2 = lcc.make_segment(Point(0,0,1),Point(1,0,1), true);
  create_attributes_3(lcc);
  lcc.template sew<1>(dh2, lcc.other_orientation(dh2));
  lcc.template sew<3>(dh1, dh2);
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 2, 2, 2, 2, 2) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0), true);
  create_attributes_3(lcc);
  lcc.template sew<1>(dh1, lcc.other_orientation(dh1));
  lcc.template sew<1>(lcc.template opposite<2>(dh1),
                      lcc.other_orientation(lcc.template opposite<2>(dh1)));
  dh2 = lcc.make_segment(Point(0,0,1),Point(1,0,1), true);
  create_attributes_3(lcc);
  lcc.template sew<1>(dh2, lcc.other_orientation(dh2));
  lcc.template sew<1>(lcc.template opposite<2>(dh2),
                      lcc.other_orientation(lcc.template opposite<2>(dh2)));
  lcc.template sew<3>(dh1, dh2);
  lcc.template sew<3>(lcc.template opposite<2>(dh1),
                      lcc.template opposite<2>(dh2));
  lcc.template contract_cell<2>(dh1);
  if ( !check_number_of_cells_3(lcc, 1, 1, 1, 2, 1) )
    return false;
  lcc.clear();

  trace_test_begin();
  dh1 = make_face_two_edges(lcc, Point(0,0,0), Point(1,0,0));
  dh2 = make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1));
  create_attributes_3(lcc);
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
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<3>(dh1,
                      make_face_two_edges(lcc, Point(0,0,1), Point(1,0,1)));
  lcc.template sew<3>(dh2,
                      make_face_two_edges(lcc, Point(1,0,1), Point(1,0,2)));
  lcc.template sew<2>(lcc.template opposite<3>(dh1),
                      lcc.template opposite<3>(dh2));
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
  create_attributes_3(lcc);
  lcc.template sew<2>(dh1, dh2);
  lcc.template sew<2>(lcc.next(dh2), dh3);
  lcc.template sew<3>(dh1, lcc.make_triangle(Point(5,5,4),Point(7,5,4),
                                             Point(6,6,4)));
  lcc.template sew<3>(dh2, lcc.make_triangle(Point(5,4,4),Point(7,4,4),
                                             Point(6,3,4)));
  lcc.template sew<3>(dh3, lcc.make_triangle(Point(5,3,4),Point(7,3,4),
                                             Point(6,0,4)));
  lcc.template sew<2>(lcc.template opposite<3>(dh1),
                      lcc.template opposite<3>(dh2));
  lcc.template sew<2>(lcc.template opposite<3>(lcc.next(dh2)),
                      lcc.template opposite<3>(dh3));
  lcc.template contract_cell<1>(lcc.previous(dh2));
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
  create_attributes_3(lcc);
  lcc.template contract_cell<3>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0));
  create_attributes_3(lcc);
  lcc.template contract_cell<3>(dh1);
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_segment(Point(0,0,0),Point(1,0,0), true);
  create_attributes_3(lcc);
  lcc.template sew<1>(dh1, lcc.other_orientation(dh1));
  lcc.template sew<1>(lcc.template opposite<2>(dh1),
                      lcc.other_orientation(lcc.template opposite<2>(dh1)));
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
  create_attributes_3(lcc);
  dh2 = lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh2))));
  lcc.template sew<3>(dh1,dh2);

  lcc.template contract_cell<1>(lcc.next(lcc.template opposite<2>(dh2)));
  lcc.template contract_cell<1>(lcc.previous(lcc.template opposite<2>(dh2)));
  lcc.template contract_cell<1>(lcc.previous(lcc.template opposite<2>
                                             (lcc.next(lcc.next(dh2)))));
  lcc.template contract_cell<1>(lcc.next(lcc.template opposite<2>
                                         (lcc.next(lcc.next(dh2)))));

  if ( !check_number_of_cells_3(lcc, 8, 16, 11, 2, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(lcc.previous(dh2)));
  if ( !check_number_of_cells_3(lcc, 8, 15, 10, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(lcc.next(lcc.next(dh2))));
  if ( !check_number_of_cells_3(lcc, 8, 14, 9, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(lcc.next(dh2)));
  if ( !check_number_of_cells_3(lcc, 8, 13, 8, 2, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(dh2));
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
  create_attributes_3(lcc);
  dh3 = lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh3))));
  lcc.template sew<3>(dh2,dh3);
  dh2 = lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh2))));
  lcc.template sew<3>(dh1,dh2);

  lcc.template contract_cell<1>(lcc.next(lcc.template opposite<2>(dh2)));
  lcc.template contract_cell<1>(lcc.previous(lcc.template opposite<2>(dh2)));
  lcc.template contract_cell<1>(lcc.previous(lcc.template opposite<2>
                                             (lcc.next(lcc.next(dh2)))));
  lcc.template contract_cell<1>(lcc.next(lcc.template opposite<2>
                                         (lcc.next(lcc.next(dh2)))));

  if ( !check_number_of_cells_3(lcc, 12, 24, 16, 3, 1) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(lcc.previous(dh2)));
  if ( !check_number_of_cells_3(lcc, 12, 23, 15, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(lcc.next(lcc.next(dh2))));
  if ( !check_number_of_cells_3(lcc, 12, 22, 14, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(lcc.next(dh2)));
  if ( !check_number_of_cells_3(lcc, 12, 21, 13, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<2>(lcc.template opposite<2>(dh2));
  if ( !check_number_of_cells_3(lcc, 12, 20, 12, 3, 1 ) )
    return false;

  trace_test_begin();
  lcc.template contract_cell<3>(dh2);
  if ( !check_number_of_cells_3(lcc, 12, 20, 11, 2, 1 ) )
    return false;
  lcc.clear();

  trace_test_begin();
  lcc.clear();
  dh1 = lcc.
      make_hexahedron(Point(0,0,0),Point(1,0,0),
                      Point(1,2,0),Point(0,2,0),
                      Point(0,3,4),Point(0,0,4),
                      Point(6,0,4),Point(6,3,4));

  dh2=lcc.insert_cell_1_in_cell_2(lcc.next(dh1), lcc.previous(dh1));
  if ( !check_number_of_cells_3(lcc, 8, 13, 7, 1, 1) )
    return false;
  lcc.template remove_cell<1>(dh2);

  dh2 = lcc.
      make_hexahedron(Point(0,0,4),Point(1,0,4),
                      Point(1,2,4),Point(0,2,4),
                      Point(0,3,8),Point(0,0,8),
                      Point(6,0,8),Point(6,3,8));
  dh3 = lcc.
      make_hexahedron(Point(5,0,4),Point(5,0,4),
                      Point(6,2,4),Point(5,2,4),
                      Point(5,3,8),Point(5,0,8),
                      Point(11,0,8),Point(11,3,8));
  create_attributes_3(lcc);
  lcc.template sew<3>(dh1,lcc.template opposite<2>
                      (lcc.next(lcc.next(lcc.template opposite<2>(dh2)))));
  lcc.template sew<3>(lcc.template opposite<2>(lcc.next(dh1)),
                      lcc.template opposite<2>(lcc.previous(dh3)));

  lcc.template close<3>();
  if ( !check_number_of_cells_3(lcc, 16, 28, 16, 4, 1) )
    return false;

  lcc.insert_cell_1_in_cell_2(lcc.next(dh1), lcc.previous(dh1));
  if ( !check_number_of_cells_3(lcc, 16, 29, 17, 4, 1) )
    return false;

  dh2=lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh1))));
  lcc.insert_cell_1_in_cell_2(dh2, lcc.next(lcc.next(dh2)));
  if ( !check_number_of_cells_3(lcc, 16, 30, 18, 4, 1) )
    return false;

  std::vector<Dart_descriptor> path;
  path.push_back(lcc.next(dh1));
  path.push_back(lcc.next(lcc.template opposite<2>(lcc.previous(dh1))));
  path.push_back(lcc.previous(dh2));
  path.push_back(lcc.next(lcc.template opposite<2>(dh2)));
  lcc.insert_cell_2_in_cell_3(path.begin(),path.end());
  if ( !check_number_of_cells_3(lcc, 16, 30, 19, 5, 1) )
    return false;

  // Test insertion between two different 2-cells
  trace_test_begin();
  lcc.clear();
  dh1 = lcc.
      make_hexahedron(Point(0,0,0),Point(1,0,0),
                      Point(1,2,0),Point(0,2,0),
                      Point(0,3,4),Point(0,0,4),
                      Point(6,0,4),Point(6,3,4));
  dh2 = lcc.
      make_hexahedron(Point(0,0,4),Point(1,0,4),
                      Point(1,2,4),Point(0,2,4),
                      Point(0,3,8),Point(0,0,8),
                      Point(6,0,8),Point(6,3,8));
  create_attributes_3(lcc);
  lcc.template sew<3>(dh1,lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh2)))));

  lcc.insert_cell_1_between_two_cells_2(lcc.template opposite<2>(dh1),
                                        lcc.template opposite<2>(lcc.next(dh1)));
  if ( !check_number_of_cells_3(lcc, 12, 21, 10, 2, 1) )
    return false;

  trace_test_begin();
  lcc.clear();
  dh1=lcc.make_hexahedron(Point(0,0,0), Point(5,0,0),
                          Point(5,5,0), Point(0,5,0),
                          Point(0,5,4), Point(0,0,4),
                          Point(5,0,4), Point(5,5,4));
  dh2=lcc.make_hexahedron(Point(5,0,0), Point(10,0,0),
                          Point(10,5,0), Point(5,5,0),
                          Point(5,5,4), Point(5,0,4),
                          Point(10,0,4), Point(10,5,4));
  dh3=lcc.make_quadrangle(Point(5,2,2), Point(5,1,2),
                          Point(5,1,1), Point(5,2,1));
  lcc.template sew<3>(lcc.template opposite<2>(lcc.next(lcc.next(dh1))),
                      lcc.other_orientation(lcc.template opposite<2>(dh2)));
  lcc.template sew<3>(dh3, lcc.make_combinatorial_polygon(4));
  create_attributes_3(lcc);

  // Create an hole in the face between the two cubes
  lcc.insert_cell_1_between_two_cells_2(lcc.template opposite<2>(lcc.next(lcc.next(dh1))),
                                        lcc.next(lcc.next(dh3)));

  if (!check_number_of_cells_3(lcc, 16, 25, 11, 2, 1) )
    return false;

  // Construction from Polyhedron_3
  {
    trace_test_begin();
    lcc.clear();
    CGAL::Polyhedron_3<typename LCC::Traits> P;
    std::ifstream in("data/head.off");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/head.off'"<<std::endl;
      return false;
    }
    in >> P;

    CGAL::import_from_polyhedron_3<LCC>(lcc,P);
    if ( !check_number_of_cells_3(lcc, 1539, 4434, 2894, 2, 2) )
      return false;

    CGAL::write_off(lcc, "copy-head.off");

    LCC lcc2; CGAL::load_off(lcc2, "copy-head.off");
    if ( !check_number_of_cells_3(lcc2, 1539, 4434, 2894, 2, 2) )
      return false;

    if (!lcc.is_isomorphic_to(lcc2, false, false, true)) // dartinfo, attrib, point
    {
      assert(false);
      return false;
    }

    lcc.clear();
    trace_test_end();
  }

  // Construction from Triangulation_3
  {
    trace_test_begin();
    CGAL::Delaunay_triangulation_3<typename LCC::Traits> T;
    std::ifstream in("data/points3D.txt");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/points.txt'"<<std::endl;
      return false;
    }
    T.insert ( std::istream_iterator < Point >(in),
               std::istream_iterator < Point >() );
    CGAL::import_from_triangulation_3<LCC>(lcc,T);
    if ( !lcc.is_valid() )
      return false;

    // Pb: the triangulation_3 is not the same on different machines ?
    if ( !check_number_of_cells_3(lcc, 286, 2386, 4200, 2100, 1) )
      return false;

    std::ofstream os("save.map");
    os<<lcc;
    os.close();

    LCC lcc2;
    std::ifstream is("save.map");
    assert(is.is_open());
    try
    {
      is>>lcc2;
    }
    catch(...)
    { // Problem during the load (boost assertion)
      std::cout<<"Problem to load combinatorial map save.map"<<std::endl;
      lcc2=lcc;
    }

    if ( !check_number_of_cells_3(lcc2, 286, 2386, 4200, 2100, 1) )
      return false;

    if (!lcc.is_isomorphic_to(lcc2, false, false, true))
    {
      std::cout<<"Different geometries after load for "
               <<typeid(LCC).name()<<std::endl;
    }

    if (!lcc.is_isomorphic_to(lcc2, false, false, false))
    {
      assert(false);
      return false;
    }

    // dual o dual is isomorphic to the initial map
    lcc.dual_points_at_barycenter(lcc2);
    LCC lcc3;
    lcc2.dual_points_at_barycenter(lcc3);

    if ( !check_number_of_cells_3(lcc3, 286, 2386, 4200, 2100, 1) )
      return false;

    if (!lcc3.is_isomorphic_to(lcc, true, true, false))
    {
      assert(false);
      return false;
    }

    lcc.clear();
    trace_test_end();
  }

  return true;
}

#endif // CGAL_LCC_3_TEST_H
