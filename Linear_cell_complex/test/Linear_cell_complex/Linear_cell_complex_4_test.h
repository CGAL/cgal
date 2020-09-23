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
#ifndef CGAL_LCC_4_TEST_H
#define CGAL_LCC_4_TEST_H

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Combinatorial_map_operations.h>
#include "Linear_cell_complex_2_test.h"

// Test orientation specialized below only for CMap. For GMap return true.
template<typename LCC, typename Map=typename LCC::Combinatorial_data_structure>
struct Test_change_orientation_LCC_4
{
  static bool run()
  { return true; }
};

template<typename LCC>
bool check_number_of_cells_4(LCC& lcc, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbvol,
                             unsigned int nbhvol, unsigned int nbcc)
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
      nbhvol!=nbc[4] || nbcc!=nbc[5])
    {
      std::cout<<"ERROR: the number of cells is not correct. We must have "
               <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbvol<<", "<<nbhvol
               <<", "<<nbcc<<") and we have"<<" ("<<nbc[0]<<", "<<nbc[1]<<", "
               <<nbc[2]<<", "<<nbc[3]<<", "<<nbc[4]<<", "<<nbc[5]<<")."
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
typename LCC::Point apoint(typename LCC::FT x, typename LCC::FT y,
                           typename LCC::FT z, typename LCC::FT t)
{
  std::vector<typename LCC::FT> tab;
  tab.push_back(x); tab.push_back(y);
  tab.push_back(z); tab.push_back(t);
  typename LCC::Point p(4,tab.begin(),tab.end());
  return p;
}

template<typename LCC>
bool test_LCC_4()
{
  LCC lcc;

  typedef typename LCC::Dart_handle Dart_handle;

  // Construction operations
  trace_test_begin();
  Dart_handle dh1=lcc.make_segment(apoint<LCC>(0,0,0,0),apoint<LCC>(1,0,0,0), true);
  Dart_handle dh2=lcc.make_segment(apoint<LCC>(2,0,0,0),apoint<LCC>(2,1,0,0), true);
  Dart_handle dh3=lcc.make_segment(apoint<LCC>(2,2,0,0),apoint<LCC>(3,1,0,0), true);
  if ( !check_number_of_cells_4(lcc, 6, 3, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<1>(dh1,dh2);
  lcc.template sew<1>(lcc.other_orientation(dh2),dh3);
  if ( !check_number_of_cells_4(lcc, 4, 3, 4, 1, 1, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh5=lcc.make_triangle(apoint<LCC>(5,5,3,0),apoint<LCC>(7,5,3,0),apoint<LCC>(6,6,3,0));
  Dart_handle dh6=lcc.make_triangle(apoint<LCC>(5,4,3,0),apoint<LCC>(7,4,3,0),apoint<LCC>(6,3,3,0));
  if ( !check_number_of_cells_4(lcc, 10, 9, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_4(lcc, 8, 8, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_4(lcc, 9, 9, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_4(lcc, 10, 12, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh9=lcc.template insert_point_in_cell<1>(dh2,apoint<LCC>(1,0,3,0));
  if ( !check_number_of_cells_4(lcc, 11, 13, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh10=lcc.template insert_point_in_cell<2>(dh6,apoint<LCC>(6,5,3,0));
  if ( !check_number_of_cells_4(lcc, 12, 16, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,apoint<LCC>(6,5.2,3,0));
  if ( !check_number_of_cells_4(lcc, 13, 17, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh12 = lcc.make_tetrahedron(apoint<LCC>(-1, 0, 0,0),apoint<LCC>(0, 2, 0,0),
                                          apoint<LCC>(1, 0, 0,0),apoint<LCC>(1, 1, 2,0));
  Dart_handle dh13 = lcc.make_tetrahedron(apoint<LCC>(0, 2, -1,0),apoint<LCC>(-1, 0, -1,0),
                                          apoint<LCC>(1, 0, -1,0),apoint<LCC>(1, 1, -3,0));
  if ( !check_number_of_cells_4(lcc, 21, 29, 18, 4, 4, 4) )
    return false;

  trace_test_begin();
  lcc.template sew<3>(dh12, dh13);
  if ( !check_number_of_cells_4(lcc, 18, 26, 17, 4, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh14=lcc.template insert_barycenter_in_cell<2>(dh12);
  if ( !check_number_of_cells_4(lcc, 19, 29, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh15=lcc.template insert_barycenter_in_cell<1>(dh14);
  if ( !check_number_of_cells_4(lcc, 20, 30, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<4>(dh12, dh13);
  if ( !check_number_of_cells_4(lcc, 19, 27, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh16=lcc.template insert_barycenter_in_cell<1>(dh15);
  if ( !check_number_of_cells_4(lcc, 20, 28, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh17=lcc.template insert_barycenter_in_cell<2>(dh16);
  if ( !check_number_of_cells_4(lcc, 21, 33, 20, 3, 3, 3) )
    return false;

  // Removal operations
  trace_test_begin();
  std::stack<Dart_handle> toremove;
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh17).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh17).end();
          it!=itend; ++it )
    toremove.push( it );
  while ( !toremove.empty() )
  {
    if (lcc.is_dart_used(toremove.top()))
      lcc.template remove_cell<1>(toremove.top());
    toremove.pop();
  }
  if ( !check_number_of_cells_4(lcc, 20, 28, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh16);
  if ( !check_number_of_cells_4(lcc, 19, 27, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template unsew<4>(dh12);
  if ( !check_number_of_cells_4(lcc, 20, 30, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh15);
  if ( !check_number_of_cells_4(lcc, 19, 29, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(lcc.next(lcc.template opposite<2>(dh14)));
  lcc.template remove_cell<1>(lcc.previous(dh14));
  lcc.template remove_cell<1>(dh14);
  if ( !check_number_of_cells_4(lcc, 18, 26, 17, 4, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template unsew<3>(dh12);
  if ( !check_number_of_cells_4(lcc, 21, 29, 18, 4, 4, 4) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<3>(dh13);
  lcc.template remove_cell<3>(dh12);
  if ( !check_number_of_cells_4(lcc, 13, 17, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(dh11);
  if ( !check_number_of_cells_4(lcc, 12, 16, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh10).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh10).end();
          it!=itend; ++it )
    toremove.push( it );
  while ( !toremove.empty() )
  {
    if (lcc.is_dart_used(toremove.top()))
      lcc.template remove_cell<1>(toremove.top());
    toremove.pop();
  }
  if ( !check_number_of_cells_4(lcc, 11, 13, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh9);
  if ( !check_number_of_cells_4(lcc, 10, 12, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh8).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh8).end();
        it!=itend; ++it )
    toremove.push( it );
  while ( !toremove.empty() )
  {
    if (lcc.is_dart_used(toremove.top()))
      lcc.template remove_cell<1>(toremove.top());
    toremove.pop();
  }
  if ( !check_number_of_cells_4(lcc, 9, 9, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<0>(dh7);
  if ( !check_number_of_cells_4(lcc, 8, 8, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<2>(dh5);
  if ( !check_number_of_cells_4(lcc, 10, 9, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<2>(dh6);
  lcc.template remove_cell<2>(dh5);
  if ( !check_number_of_cells_4(lcc, 4, 3, 4, 1, 1, 1) )
    return false;

  trace_test_begin();
  lcc.template unsew<1>(dh1);
  if ( !check_number_of_cells_4(lcc, 5, 3, 5, 2, 2, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<1>(lcc.other_orientation(dh2));
  if ( !check_number_of_cells_4(lcc, 6, 3, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  lcc.template remove_cell<1>(dh1);
  lcc.template remove_cell<1>(dh2);
  lcc.template remove_cell<1>(dh3);
  if ( !check_number_of_cells_4(lcc, 0, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = lcc.make_tetrahedron(apoint<LCC>(-1, 0, 0,0),apoint<LCC>(0, 2, 0,0),
                             apoint<LCC>(1, 0, 0,0),apoint<LCC>(1, 1, 2,0));
  dh2 = lcc.make_tetrahedron(apoint<LCC>(0, 2, -1,0),apoint<LCC>(-1, 0, -1,0),
                             apoint<LCC>(1, 0, -1,0),apoint<LCC>(1, 1, -3,0));

  if ( !lcc.template is_sewable<4>(dh1, dh2) )
  {
    std::cout<<"ERROR: the two 3-cells are not sewable."<<std::endl;
    assert(false);
    return false;
  }
  trace_test_end();

  trace_test_begin();
  dh3 = lcc.template opposite<2>(dh1);
  dh5 = lcc.template opposite<2>(lcc.next(dh1));

  lcc.template unsew<2>(dh3);
  lcc.template unsew<2>(dh5);
  lcc.template sew<2>(dh1, dh5);
  lcc.template sew<2>(lcc.next(dh1), dh3);

  if ( lcc.template is_sewable<4>(dh1, dh2) )
  {
    std::cout<<"ERROR: the two 3-cells are sewable."<<std::endl;
    assert(false);
    return false;
  }
  trace_test_end();

  trace_test_begin();
  lcc.clear();
  dh1 = lcc.make_tetrahedron(apoint<LCC>(-1, 0, 0,0),apoint<LCC>(0, 2, 0,0),
                             apoint<LCC>(1, 0, 0,0),apoint<LCC>(1, 1, 2,0));
  dh2 = lcc.make_tetrahedron(apoint<LCC>(0, 2, -1,0),apoint<LCC>(-1, 0, -1,0),
                             apoint<LCC>(1, 0, -1,0),apoint<LCC>(1, 1, -3,0));
  dh3 = lcc.make_tetrahedron(apoint<LCC>(0, 2, -4,0),apoint<LCC>(-1, 0, -4,0),
                             apoint<LCC>(1, 0, -4,0),apoint<LCC>(1, 1, -5,0));
  lcc.template sew<3>(dh1, dh2);
  lcc.template sew<4>(dh1, dh3);

  LCC lcc2(lcc);
  if ( !lcc2.is_valid() ) { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  lcc2=lcc;
  if ( !lcc2.is_valid() ) { assert(false); return false; }
  if ( !lcc2.is_isomorphic_to(lcc) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  lcc.clear();
  dh1 = lcc.
      make_hexahedron(apoint<LCC>(0,0,0,0),apoint<LCC>(1,0,0,0),
                      apoint<LCC>(1,2,0,0),apoint<LCC>(0,2,0,0),
                      apoint<LCC>(0,3,4,0),apoint<LCC>(0,0,4,0),
                      apoint<LCC>(6,0,4,0),apoint<LCC>(6,3,4,0));
  dh2 = lcc.
      make_hexahedron(apoint<LCC>(0,0,4,0),apoint<LCC>(1,0,4,0),
                      apoint<LCC>(1,2,4,0),apoint<LCC>(0,2,4,0),
                      apoint<LCC>(0,3,8,0),apoint<LCC>(0,0,8,0),
                      apoint<LCC>(6,0,8,0),apoint<LCC>(6,3,8,0));
  dh3 = lcc.
      make_hexahedron(apoint<LCC>(5,0,4,0),apoint<LCC>(5,0,4,0),
                      apoint<LCC>(6,2,4,0),apoint<LCC>(5,2,4,0),
                      apoint<LCC>(5,3,8,0),apoint<LCC>(5,0,8,0),
                      apoint<LCC>(11,0,8,0),apoint<LCC>(11,3,8,0));
  lcc.template sew<3>(dh1,lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh2)))));
  lcc.template sew<3>(lcc.template opposite<2>(lcc.next(dh1)), lcc.template opposite<2>(lcc.previous(dh3)));

  lcc.template close<3>();
  lcc.template close<4>();

  if ( !check_number_of_cells_4(lcc, 16, 28, 16, 4, 2, 1) )
    return false;

  lcc.insert_cell_1_in_cell_2(lcc.next(dh1), Alpha1<LCC>::run(lcc, lcc.previous(dh1)));
  dh2=lcc.template opposite<2>(lcc.next(lcc.next(lcc.template opposite<2>(dh1))));
  lcc.insert_cell_1_in_cell_2(dh2, Alpha1<LCC>::run(lcc, lcc.next(lcc.next(dh2))));

  std::vector<Dart_handle> path;
  path.push_back(lcc.next(dh1));
  path.push_back(lcc.next(lcc.template opposite<2>(lcc.previous(dh1))));
  path.push_back(lcc.previous(dh2));
  path.push_back(lcc.next(lcc.template opposite<2>(dh2)));
  lcc.insert_cell_2_in_cell_3(path.begin(),path.end());
  if ( !check_number_of_cells_4(lcc, 16, 30, 19, 5, 2, 1) )
    return false;

  if ( !Test_change_orientation_LCC_4<LCC>::run() )
    return false;

  return true;
}

template<typename LCC>
struct Test_change_orientation_LCC_4<LCC, CGAL::Combinatorial_map_tag>
{
  static bool run()
  {
    LCC lcc;

    trace_test_begin();
    typename LCC::Dart_handle dh1 = lcc.make_tetrahedron(apoint<LCC>(-1, 0, 0,0),
                                                         apoint<LCC>(0, 2, 0,0),
                                                         apoint<LCC>(1, 0, 0,0),
                                                         apoint<LCC>(1, 1, 2,0));
    typename LCC::Dart_handle dh2 = lcc.make_tetrahedron(apoint<LCC>(0, 2, -1,0),
                                                         apoint<LCC>(-1, 0, -1,0),
                                                         apoint<LCC>(1, 0, -1,0),
                                                         apoint<LCC>(1, 1, -3,0));
    typename LCC::Dart_handle dh3 = lcc.make_tetrahedron(apoint<LCC>(0, 2, -4,0),
                                                         apoint<LCC>(-1, 0, -4,0),
                                                         apoint<LCC>(1, 0, -4,0),
                                                         apoint<LCC>(1, 1, -5,0));
    lcc.make_tetrahedron(apoint<LCC>(0, 2, -10,0),
                         apoint<LCC>(-1, 0, -10,0),
                         apoint<LCC>(1, 0, -10,0),
                         apoint<LCC>(1, 1, -12,0));
    lcc.template sew<3>(dh1, dh2);
    lcc.template sew<4>(dh1, dh3);

    LCC lcc2(lcc);

    trace_test_begin();
    lcc.reverse_orientation();
    if ( !lcc.is_valid() ) { assert(false); return false; }
    if ( lcc2.is_isomorphic_to(lcc) )
    { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc, false, false, false) )
    { assert(false); return false; }
    trace_test_end();

    trace_test_begin();
    lcc.reverse_orientation();
    if ( !lcc.is_valid() ) { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc) )
    { assert(false); return false; }
    trace_test_end();

    trace_test_begin();
    lcc.reverse_orientation_connected_component(dh1);
    if ( !lcc.is_valid() ) { assert(false); return false; }
    if ( lcc2.is_isomorphic_to(lcc) )
    { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc, false, false, false) )
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
};

#endif // CGAL_LCC_4_TEST_H
