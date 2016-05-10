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
#ifndef CGAL_LCC_2_TEST_H
#define CGAL_LCC_2_TEST_H

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <fstream>

// #define LCC_TRACE_TEST_BEGIN 1

void trace_test_begin()
{
#ifdef LCC_TRACE_TEST_BEGIN
  static unsigned int nbtest = 0;

  std::cout<<"Test "<<nbtest++<<" ..."<<std::flush;
#endif
}

void trace_test_end()
{
#ifdef LCC_TRACE_TEST_BEGIN
  std::cout<<"Ok."<<std::endl;
#endif
}

void trace_display_msg(const char*
#ifdef LCC_TRACE_TEST_BEGIN
                       msg
#endif
                       )
{
#ifdef LCC_TRACE_TEST_BEGIN
  std::cout<<"***************** "<<msg<<"***************** "<<std::endl;
#endif
}

template<typename LCC>
bool check_number_of_cells_2(LCC& lcc, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbcc)
{
  if ( !lcc.is_valid() )
    {
      std::cout<<"ERROR: the lcc is not valid."<<std::endl;
      assert(false);
      return false;
    }

  std::vector<unsigned int> nbc;
  nbc=lcc.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbcc!=nbc[3])
    {
      std::cout<<"ERROR: the number of cells is not correct. We must have "
               <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbcc<<") and we have"
               <<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "<<nbc[3]<<")."
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
void display_lcc(LCC& lcc)
{
  unsigned int nb = 0;
  for ( typename LCC::Dart_range::const_iterator it=lcc.darts().begin();
        it!=lcc.darts().end(); ++it)
  {
    std::cout << " dart " << &(*it) << "; beta[i]=";
    for ( unsigned int i=0; i<=LCC::dimension; ++i)
    {
      std::cout << &(*it->beta(i)) << ",\t";
      if (it->is_free(i)) std::cout << "\t";
    }
    std::cout<<it->template attribute<0>()->point();
    std::cout << std::endl;
    ++nb;
  }
  std::cout << "Number of darts: " << nb <<"(sizeofdarts="
     <<lcc.number_of_darts()<<")" << std::endl;
}

template<typename LCC>
bool test_LCC_2()
{
  LCC lcc;

  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;

  // Construction operations
  trace_test_begin();
  Dart_handle dh1=lcc.make_segment(Point(0,0),Point(1,0));
  Dart_handle dh2=lcc.make_segment(Point(2,0),Point(2,1));
  Dart_handle dh3=lcc.make_segment(Point(2,2),Point(3,1));
  if ( !check_number_of_cells_2(lcc, 6, 3, 6, 3) )
    return false;

  typename LCC::Vertex_attribute_handle vh=lcc.template attribute<0>(dh1);
  if (!lcc.template is_attribute_used<0>(vh)) return false;

  trace_test_begin();
  lcc.template sew<0>(dh2,dh1);
  lcc.template sew<1>(dh2,dh3);
  if ( !check_number_of_cells_2(lcc, 4, 3, 4, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh5=lcc.make_triangle(Point(5,5),Point(7,5),Point(6,6));
  Dart_handle dh6=lcc.make_triangle(Point(5,4),Point(7,4),Point(6,3));
  if ( !check_number_of_cells_2(lcc, 10, 9, 6, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_2(lcc, 8, 8, 6, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_2(lcc, 9, 9, 6, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_2(lcc, 10, 12, 8, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0));
  if ( !check_number_of_cells_2(lcc, 11, 13, 8, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5));
  if ( !check_number_of_cells_2(lcc, 12, 16, 10, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2));
  if ( !check_number_of_cells_2(lcc, 13, 17, 10, 2) )
    return false;

  // Removal operations
  trace_test_begin();
  lcc.template remove_cell<1>(dh11);
  if ( !check_number_of_cells_2(lcc, 12, 16, 10, 2) )
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
  if ( !check_number_of_cells_2(lcc, 11, 13, 8, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh9);
  if ( !check_number_of_cells_2(lcc, 10, 12, 8, 2) )
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
  if ( !check_number_of_cells_2(lcc, 9, 9, 6, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh7);
  if ( !check_number_of_cells_2(lcc, 8, 8, 6, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<2>(dh5);
  if ( !check_number_of_cells_2(lcc, 10, 9, 6, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<2>(dh6);
  lcc.template remove_cell<2>(dh5);
  if ( !check_number_of_cells_2(lcc, 4, 3, 4, 1) )
    return false;

  trace_test_begin();
  lcc.template unsew<1>(dh2);
  if ( !check_number_of_cells_2(lcc, 5, 3, 5, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<0>(dh2);
  if ( !check_number_of_cells_2(lcc, 6, 3, 6, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(dh1);
  lcc.template remove_cell<1>(dh2);
  lcc.template remove_cell<1>(dh3);
  if ( !check_number_of_cells_2(lcc, 0, 0, 0, 0) )
    return false;

  if (lcc.template is_attribute_used<0>(vh)) return false;

  {
    trace_test_begin();
    std::ifstream in("data/graph.txt");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/graph.txt'"<<std::endl;
      return false;
    }
    CGAL:: import_from_plane_graph<LCC>(lcc,in);
    if ( !check_number_of_cells_2(lcc, 61, 160, 101, 1) )
      return false;
    lcc.clear();
  }

  trace_test_begin();
  lcc.clear();
  dh1=lcc.make_triangle(Point(5,5),Point(7,5),Point(6,6));
  dh2=lcc.make_triangle(Point(5,4),Point(7,4),Point(6,3));
  lcc.template sew<2>(dh1,dh2);

  LCC lcc2(lcc);
  if ( !lcc.is_valid() ) { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  lcc.reverse_orientation();
  if ( !lcc.is_valid() ) { assert(false); return false; }
  if ( lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc, false) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  lcc.reverse_orientation();
  if ( !lcc.is_valid() ) { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc, false) )
  { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  lcc.reverse_orientation_connected_component(dh1);
  if ( !lcc.is_valid() ) { assert(false); return false; }
  if ( lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc, false) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  lcc.reverse_orientation_connected_component(dh1);
  if ( !lcc.is_valid() ) { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  trace_test_end();

  return true;
}

#endif // CGAL_LCC_2_TEST_H
