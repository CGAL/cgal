// Copyright (c) 2014 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GMAP_4_TEST_H
#define CGAL_GMAP_4_TEST_H

#include <CGAL/Generalized_map_operations.h>
#include "Generalized_map_2_test.h"

template<typename GMAP>
bool check_number_of_cells_4(GMAP& gmap, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbvol,
                             unsigned int nbhvol, unsigned int nbcc)
{
  if ( !gmap.is_valid() )
    {
      std::cout<<"ERROR: the gmap is not valid."<<std::endl;
      assert(false);
      return false;
    }

  std::vector<unsigned int> nbc;
  nbc=gmap.count_all_cells();

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

  trace_test_end();

  return true;
}

template<typename GMAP>
bool test_GMAP_4()
{
  GMAP gmap;

  typedef typename GMAP::Dart_handle Dart_handle;

  // Construction operations
  trace_test_begin();
  Dart_handle dh1=CGAL::make_edge(gmap);
  Dart_handle dh2=CGAL::make_edge(gmap);
  Dart_handle dh3=CGAL::make_edge(gmap);
  if ( !check_number_of_cells_4(gmap, 6, 3, 3, 3, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template sew<1>(dh1,dh2);
  gmap.template sew<1>(gmap.alpha(dh2,0),dh3);
  if ( !check_number_of_cells_4(gmap, 4, 3, 1, 1, 1, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh5=CGAL::make_combinatorial_polygon(gmap, 3);
  Dart_handle dh6=CGAL::make_combinatorial_polygon(gmap, 3);
  if ( !check_number_of_cells_4(gmap, 10, 9, 3, 3, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_4(gmap, 8, 8, 3, 2, 2, 2) )
    return false;

  trace_test_begin();
  gmap.clear();
  Dart_handle dh7=CGAL::make_combinatorial_hexahedron(gmap); // f1
  Dart_handle dh8=gmap.template alpha<2,1,0,1,2>(dh7); // f2 opposite to f1
  Dart_handle dh9=gmap.template alpha<2>(dh7); // face incident to f1 and d2

  CGAL::remove_cell<GMAP,2>(gmap, dh7);
  if ( !check_number_of_cells_4(gmap, 8, 12, 5, 1, 1, 1) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,2>(gmap, dh8);
  if ( !check_number_of_cells_4(gmap, 8, 12, 4, 1, 1, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh10=CGAL::make_combinatorial_hexahedron(gmap);
  gmap.template sew<3>(dh9,dh10);
  if ( !check_number_of_cells_4(gmap, 12, 20, 9, 2, 1, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh11=CGAL::make_combinatorial_hexahedron(gmap);
  gmap.template sew<4>(dh10,dh11);
  if ( !check_number_of_cells_4(gmap, 12, 20, 9, 2, 2, 1) )
    return false;

  trace_test_begin();
  GMAP gmap2(gmap);
  if ( !check_number_of_cells_4(gmap2, 12, 20, 9, 2, 2, 1) )
    return false;
  if ( !gmap.is_isomorphic_to(gmap2) )
  {
    std::cout<<"Error: gmap and gmap2 are not isomorphic (after copy).\n";
    assert(false);
    return false;
  }

  trace_test_begin();
  gmap.template close<2>();
  if ( !check_number_of_cells_4(gmap, 12, 20, 11, 2, 2, 1) )
    return false;
  if ( !CGAL::is_volume_combinatorial_hexahedron(gmap, dh9) )
  {
    std::cout<<"Error: the closed volume is not a combinatorial hexahedron.\n";
    assert(false);
    return false;
  }

  trace_test_begin();
  gmap.template close<3>();
  if ( !check_number_of_cells_4(gmap, 12, 20, 11, 4, 2, 1) )
    return false;

  trace_test_begin();
  gmap.template close<4>();
  if ( !check_number_of_cells_4(gmap, 12, 20, 11, 4, 3, 1) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,4>(gmap, gmap.alpha(dh9, 2, 3, 4));
  if ( !check_number_of_cells_4(gmap, 12, 20, 11, 4, 2, 1) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,3>(gmap, gmap.alpha(dh9, 2, 3));
  CGAL::remove_cell<GMAP,3>(gmap, gmap.alpha(dh10, 2, 4, 3));
  if ( !check_number_of_cells_4(gmap, 12, 20, 11, 2, 2, 1) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,2>(gmap, gmap.alpha(dh9, 2));
  CGAL::remove_cell<GMAP,2>(gmap, gmap.alpha(dh9, 1, 0, 1, 2));
  if ( !check_number_of_cells_4(gmap, 12, 20, 9, 2, 2, 1) )
    return false;

  if ( !gmap.is_isomorphic_to(gmap2) )
  {
    std::cout<<"Error: gmap and gmap2 are not isomorphic (after close and removals).\n";
    assert(false);
    return false;
  }


/*  trace_test_begin();
  Dart_handle dh7=gmap.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_4(gmap, 9, 9, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh8=gmap.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_4(gmap, 10, 12, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh9=gmap.template insert_point_in_cell<1>(dh2,apoint<GMAP>(1,0,3,0));
  if ( !check_number_of_cells_4(gmap, 11, 13, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh10=gmap.template insert_point_in_cell<2>(dh6,apoint<GMAP>(6,5,3,0));
  if ( !check_number_of_cells_4(gmap, 12, 16, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh11=gmap.insert_dangling_cell_1_in_cell_2(dh8,apoint<GMAP>(6,5.2,3,0));
  if ( !check_number_of_cells_4(gmap, 13, 17, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  Dart_handle dh12 = gmap.make_tetrahedron(apoint<GMAP>(-1, 0, 0,0),apoint<GMAP>(0, 2, 0,0),
                                          apoint<GMAP>(1, 0, 0,0),apoint<GMAP>(1, 1, 2,0));
  Dart_handle dh13 = gmap.make_tetrahedron(apoint<GMAP>(0, 2, -1,0),apoint<GMAP>(-1, 0, -1,0),
                                          apoint<GMAP>(1, 0, -1,0),apoint<GMAP>(1, 1, -3,0));
  if ( !check_number_of_cells_4(gmap, 21, 29, 18, 4, 4, 4) )
    return false;

  trace_test_begin();
  gmap.template sew<3>(dh12, dh13);
  if ( !check_number_of_cells_4(gmap, 18, 26, 17, 4, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh14=gmap.template insert_barycenter_in_cell<2>(dh12);
  if ( !check_number_of_cells_4(gmap, 19, 29, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh15=gmap.template insert_barycenter_in_cell<1>(dh14);
  if ( !check_number_of_cells_4(gmap, 20, 30, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template sew<4>(dh12, dh13);
  if ( !check_number_of_cells_4(gmap, 19, 27, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh16=gmap.template insert_barycenter_in_cell<1>(dh15);
  if ( !check_number_of_cells_4(gmap, 20, 28, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  Dart_handle dh17=gmap.template insert_barycenter_in_cell<2>(dh16);
  if ( !check_number_of_cells_4(gmap, 21, 33, 20, 3, 3, 3) )
    return false;

  // Removal operations
  trace_test_begin();
  std::stack<Dart_handle> toremove;
  for ( typename GMAP::template Dart_of_cell_range<0,2>::iterator
          it=gmap.template darts_of_cell<0,2>(dh17).begin(),
          itend=gmap.template darts_of_cell<0,2>(dh17).end();
          it!=itend; ++it )
    toremove.push( it );
  while ( !toremove.empty() )
  {
    CGAL::remove_cell<GMAP,1>(gmap, toremove.top());
    toremove.pop();
  }
  if ( !check_number_of_cells_4(gmap, 20, 28, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,0>(gmap, dh16);
  if ( !check_number_of_cells_4(gmap, 19, 27, 16, 3, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template unsew<4>(dh12);
  if ( !check_number_of_cells_4(gmap, 20, 30, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,0>(gmap, dh15);
  if ( !check_number_of_cells_4(gmap, 19, 29, 19, 4, 3, 3) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,1>(gmap, gmap.beta(dh14,2,1));
  CGAL::remove_cell<GMAP,1>(gmap, gmap.beta(dh14,0));
  CGAL::remove_cell<GMAP,1>(gmap, dh14);
  if ( !check_number_of_cells_4(gmap, 18, 26, 17, 4, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template unsew<3>(dh12);
  if ( !check_number_of_cells_4(gmap, 21, 29, 18, 4, 4, 4) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,3>(gmap, dh13);
  CGAL::remove_cell<GMAP,3>(gmap, dh12);
  if ( !check_number_of_cells_4(gmap, 13, 17, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,1>(gmap, dh11);
  if ( !check_number_of_cells_4(gmap, 12, 16, 10, 2, 2, 2) )
    return false;

  trace_test_begin();
  for ( typename GMAP::template Dart_of_cell_range<0,2>::iterator
          it=gmap.template darts_of_cell<0,2>(dh10).begin(),
          itend=gmap.template darts_of_cell<0,2>(dh10).end();
          it!=itend; ++it )
    toremove.push( it );
  while ( !toremove.empty() )
  {
    CGAL::remove_cell<GMAP,1>(gmap, toremove.top());
    toremove.pop();
  }
  if ( !check_number_of_cells_4(gmap, 11, 13, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,0>(gmap, dh9);
  if ( !check_number_of_cells_4(gmap, 10, 12, 8, 2, 2, 2) )
    return false;

  trace_test_begin();
  for ( typename GMAP::template Dart_of_cell_range<0,2>::iterator
          it=gmap.template darts_of_cell<0,2>(dh8).begin(),
          itend=gmap.template darts_of_cell<0,2>(dh8).end();
        it!=itend; ++it )
    toremove.push( it );
  while ( !toremove.empty() )
  {
    CGAL::remove_cell<GMAP,1>(gmap, toremove.top());
    toremove.pop();
  }
  if ( !check_number_of_cells_4(gmap, 9, 9, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,0>(gmap, dh7);
  if ( !check_number_of_cells_4(gmap, 8, 8, 6, 2, 2, 2) )
    return false;

  trace_test_begin();
  gmap.template unsew<2>(dh5);
  if ( !check_number_of_cells_4(gmap, 10, 9, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,2>(gmap, dh6);
  CGAL::remove_cell<GMAP,2>(gmap, dh5);
  if ( !check_number_of_cells_4(gmap, 4, 3, 4, 1, 1, 1) )
    return false;

  trace_test_begin();
  gmap.template unsew<1>(dh2);
  if ( !check_number_of_cells_4(gmap, 5, 3, 5, 2, 2, 2) )
    return false;

  trace_test_begin();
  gmap.template unsew<0>(dh2);
  if ( !check_number_of_cells_4(gmap, 6, 3, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  CGAL::remove_cell<GMAP,1>(gmap, dh1);
  CGAL::remove_cell<GMAP,1>(gmap, dh2);
  CGAL::remove_cell<GMAP,1>(gmap, dh3);
  if ( !check_number_of_cells_4(gmap, 0, 0, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  dh1 = gmap.make_tetrahedron(apoint<GMAP>(-1, 0, 0,0),apoint<GMAP>(0, 2, 0,0),
                             apoint<GMAP>(1, 0, 0,0),apoint<GMAP>(1, 1, 2,0));
  dh2 = gmap.make_tetrahedron(apoint<GMAP>(0, 2, -1,0),apoint<GMAP>(-1, 0, -1,0),
                             apoint<GMAP>(1, 0, -1,0),apoint<GMAP>(1, 1, -3,0));

  if ( !gmap.template is_sewable<4>(dh1, dh2) )
  {
    std::cout<<"ERROR: the two 3-cells are not sewable."<<std::endl;
    assert(false);
    return false;
  }
  trace_test_end();

  trace_test_begin();
  dh3 = gmap.beta(dh1,2);
  dh5 = gmap.template beta<1, 2>(dh1);

  gmap.template unsew<2>(dh3);
  gmap.template unsew<2>(dh5);
  gmap.template sew<2>(dh1, dh5);
  gmap.template sew<2>(gmap.beta(dh1,1), dh3);

  if ( gmap.template is_sewable<4>(dh1, dh2) )
  {
    std::cout<<"ERROR: the two 3-cells are sewable."<<std::endl;
    assert(false);
    return false;
  }
  trace_test_end();

  trace_test_begin();
  gmap.clear();
  dh1 = gmap.make_tetrahedron(apoint<GMAP>(-1, 0, 0,0),apoint<GMAP>(0, 2, 0,0),
                             apoint<GMAP>(1, 0, 0,0),apoint<GMAP>(1, 1, 2,0));
  dh2 = gmap.make_tetrahedron(apoint<GMAP>(0, 2, -1,0),apoint<GMAP>(-1, 0, -1,0),
                             apoint<GMAP>(1, 0, -1,0),apoint<GMAP>(1, 1, -3,0));
  dh3 = gmap.make_tetrahedron(apoint<GMAP>(0, 2, -4,0),apoint<GMAP>(-1, 0, -4,0),
                             apoint<GMAP>(1, 0, -4,0),apoint<GMAP>(1, 1, -5,0));
  gmap.template sew<3>(dh1, dh2);
  gmap.template sew<4>(dh1, dh3);

  GMAP gmap2(gmap);
  if ( !gmap.is_valid() ) { assert(false); return false; }
  if ( !gmap2.is_isomorphic_to(gmap) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  gmap.reverse_orientation();
  if ( !gmap.is_valid() ) { assert(false); return false; }
  if ( gmap2.is_isomorphic_to(gmap) )
  { assert(false); return false; }
  if ( !gmap2.is_isomorphic_to(gmap, false) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  gmap.reverse_orientation();
  if ( !gmap.is_valid() ) { assert(false); return false; }
  if ( !gmap2.is_isomorphic_to(gmap, false) )
  { assert(false); return false; }
  if ( !gmap2.is_isomorphic_to(gmap) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  gmap.reverse_orientation_connected_component(dh1);
  if ( !gmap.is_valid() ) { assert(false); return false; }
  if ( gmap2.is_isomorphic_to(gmap) )
  { assert(false); return false; }
  if ( !gmap2.is_isomorphic_to(gmap, false) )
  { assert(false); return false; }
  trace_test_end();

  trace_test_begin();
  gmap.reverse_orientation_connected_component(dh1);
  if ( !gmap.is_valid() ) { assert(false); return false; }
  if ( !gmap2.is_isomorphic_to(gmap) )
  { assert(false); return false; }
  trace_test_end();

  //    import_from_polyhedron<GMAP>(gmap,ap);

        gmap.clear();

  //      import_from_plane_graph<GMAP>(gmap,ais);
  */
        gmap.clear();

  return true;
}

#endif // CGAL_GMAP_4_TEST_H
