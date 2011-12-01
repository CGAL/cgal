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
  
  return true;
}

template<typename LCC>
bool test_LCC_2()
{
  LCC lcc;

  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  typedef typename LCC::Vector Vector;

  // Construction operations
  Dart_handle dh1=lcc.make_segment(Point(0,0),Point(1,0));
  Dart_handle dh2=lcc.make_segment(Point(2,0),Point(2,1));
  Dart_handle dh3=lcc.make_segment(Point(2,2),Point(3,1));
  if ( !check_number_of_cells_2(lcc, 6, 3, 6, 3) )
    return false;
  
  lcc.template sew<0>(dh2,dh1);
  lcc.template sew<1>(dh2,dh3);
  if ( !check_number_of_cells_2(lcc, 4, 3, 4, 1) )
    return false;

  Dart_handle dh5=lcc.make_triangle(Point(5,5),Point(7,5),Point(6,6));
  Dart_handle dh6=lcc.make_triangle(Point(5,4),Point(7,4),Point(6,3));    
  if ( !check_number_of_cells_2(lcc, 10, 9, 6, 3) )
    return false;

  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_2(lcc, 8, 8, 6, 2) )
    return false;

  Dart_handle dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_2(lcc, 9, 9, 6, 2) )
    return false;

  Dart_handle dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_2(lcc, 10, 12, 8, 2) )
    return false;

  Dart_handle dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0));
  if ( !check_number_of_cells_2(lcc, 11, 13, 8, 2) )
    return false;

  Dart_handle dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5));
  if ( !check_number_of_cells_2(lcc, 12, 16, 10, 2) )
    return false;

  Dart_handle dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2));
  if ( !check_number_of_cells_2(lcc, 13, 17, 10, 2) )
    return false;

  // Removal operations
  CGAL::remove_cell<LCC,1>(lcc, dh11);
  if ( !check_number_of_cells_2(lcc, 12, 16, 10, 2) )
    return false;

  std::vector<Dart_handle> toremove;
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh10).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh10).end();
          it!=itend; ++it )
    toremove.push_back( it );
  
  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    CGAL::remove_cell<LCC,1>(lcc, *it);
  toremove.clear();
  if ( !check_number_of_cells_2(lcc, 11, 13, 8, 2) )
    return false;

  CGAL::remove_cell<LCC,0>(lcc, dh9);
  if ( !check_number_of_cells_2(lcc, 10, 12, 8, 2) )
    return false;

  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh8).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh8).end();
        it!=itend; ++it )
    toremove.push_back( it );
  
  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    CGAL::remove_cell<LCC,1>(lcc, *it);
  toremove.clear();
  if ( !check_number_of_cells_2(lcc, 9, 9, 6, 2) )
    return false;
  
  CGAL::remove_cell<LCC,0>(lcc, dh7);
  if ( !check_number_of_cells_2(lcc, 8, 8, 6, 2) )
    return false;

  lcc.template unsew<2>(dh5);
  if ( !check_number_of_cells_2(lcc, 10, 9, 6, 3) )
    return false;

  CGAL::remove_cell<LCC,2>(lcc, dh6);
  CGAL::remove_cell<LCC,2>(lcc, dh5);
  if ( !check_number_of_cells_2(lcc, 4, 3, 4, 1) )
    return false;

  lcc.template unsew<1>(dh2);
  if ( !check_number_of_cells_2(lcc, 5, 3, 5, 2) )
    return false;

  lcc.template unsew<0>(dh2);
  if ( !check_number_of_cells_2(lcc, 6, 3, 6, 3) )
    return false;

  CGAL::remove_cell<LCC,1>(lcc, dh1);
  CGAL::remove_cell<LCC,1>(lcc, dh2);
  CGAL::remove_cell<LCC,1>(lcc, dh3);
  if ( !check_number_of_cells_2(lcc, 0, 0, 0, 0) )
    return false;

  {
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
    
  return true;
}

#endif // CGAL_LCC_2_TEST_H
