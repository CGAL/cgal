// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_2_TEST
#define CGAL_GENERALIZED_MAP_2_TEST 1

#include <CGAL/Generalized_map_operations.h>
#include <CGAL/Random.h>

#include <iostream>

using namespace std;

// #define GMAP_TRACE_TEST_BEGIN 1

void trace_test_begin()
{
#ifdef GMAP_TRACE_TEST_BEGIN
  static unsigned int nbtest = 0;

  std::cout<<"Test "<<nbtest++<<" ..."<<std::flush;
#endif
}

void trace_test_end()
{
#ifdef GMAP_TRACE_TEST_BEGIN
  std::cout<<"Ok."<<std::endl;
#endif
}

void trace_display_msg(const char*
#ifdef GMAP_TRACE_TEST_BEGIN
                       msg
#endif
                       )
{
#ifdef GMAP_TRACE_TEST_BEGIN
  std::cout<<"***************** "<<msg<<"***************** "<<std::endl;
#endif
}

template<typename GMap, typename Info=typename GMap::Dart_info>
struct InitDartInfo
{
  static void run(GMap& gmap)
  {
    long long int nb=0;
    for(typename GMap::Dart_range::iterator it=gmap.darts().begin(),
        itend=gmap.darts().end(); it!=itend; ++it)
    {
      nb=CGAL::get_default_random().get_int(0,20000);
      gmap.info(it)=Info(nb);
    }
  }
};

template<typename GMap>
struct InitDartInfo<GMap, CGAL::Void>
{
  static void run(GMap&)
  {}
};

template<typename GMAP>
bool check_number_of_cells_2(GMAP& gmap, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbcc)
{
  if ( !gmap.is_valid() )
    {
      std::cout<<"ERROR: the gmap is not valid."<<std::endl;
      assert(false);
      return false;
    }

  std::vector<unsigned int> nbc;
  nbc=gmap.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbcc!=nbc[3])
    {
      std::cout<<"ERROR: the number of cells is not correct. We must have "
               <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbcc<<") and we have"
               <<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "<<nbc[3]<<")."
               <<std::endl;
      assert(false);
      return false;
    }

  trace_test_end();

  return true;
}

template<class Gmap>
bool test_GMAP_2()
{
  Gmap gmap;

  typedef typename Gmap::Dart_handle Dart_handle;

  // Construction operations
  trace_test_begin();
  Dart_handle dh1=gmap.make_edge();
  Dart_handle dh2=gmap.make_edge();
  Dart_handle dh3=gmap.make_edge();
  if ( !check_number_of_cells_2(gmap, 6, 3, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template sew<1>(dh1,dh2);
  gmap.template sew<1>(gmap.alpha(dh2, 0),dh3);
  if ( !check_number_of_cells_2(gmap, 4, 3, 1, 1) )
    return false;

  trace_test_begin();
  Dart_handle dh5=gmap.make_combinatorial_polygon(3);
  Dart_handle dh6=gmap.make_combinatorial_polygon(3);
  if ( !check_number_of_cells_2(gmap, 10, 9, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_2(gmap, 8, 8, 3, 2) )
    return false;

  trace_test_begin();
  dh5=gmap.template alpha<1>(dh5);
  dh6=gmap.template alpha<1>(dh6);
  gmap.template contract_cell<1>(gmap.template alpha<1>(dh5));
  if ( !check_number_of_cells_2(gmap, 8, 7, 3, 3) )
    return false;

  trace_test_begin();
  gmap.template contract_cell<2>(dh6);
  if ( !check_number_of_cells_2(gmap, 6, 5, 2, 2) )
    return false;

  trace_test_begin();
  gmap.template contract_cell<1>(gmap.template alpha<1>(dh5));
  if ( !check_number_of_cells_2(gmap, 5, 4, 2, 2) )
    return false;

  trace_test_begin();
  gmap.template contract_cell<1>(dh5);
  if ( !check_number_of_cells_2(gmap, 4, 3, 1, 1) )
    return false;

  trace_test_begin();
  gmap.template contract_cell<1>(dh2);
  if ( !check_number_of_cells_2(gmap, 3, 2, 1, 1) )
    return false;

  trace_test_begin();
  gmap.template contract_cell<1>(dh1);
  if ( !check_number_of_cells_2(gmap, 2, 1, 1, 1) )
    return false;

  trace_test_begin();
  gmap.template contract_cell<1>(dh3);
  if ( !check_number_of_cells_2(gmap, 0, 0, 0, 0) )
    return false;

  trace_test_begin();
  Dart_handle dh7=gmap.make_combinatorial_hexahedron(); // f1
  Dart_handle dh8=gmap.template alpha<2,1,0,1,2>(dh7); // f2 opposite to f1
  Dart_handle dh9=gmap.template alpha<2>(dh7); // face incident to f1 and d2

  gmap.template remove_cell<2>(dh7);
  if ( !check_number_of_cells_2(gmap, 8, 12, 5, 1) )
    return false;

  trace_test_begin();
  gmap.template remove_cell<2>(dh8);
  if ( !check_number_of_cells_2(gmap, 8, 12, 4, 1) )
    return false;

  trace_test_begin();
  gmap.template close<2>();
  if ( !check_number_of_cells_2(gmap, 8, 12, 6, 1) )
    return false;
  if ( !gmap.is_volume_combinatorial_hexahedron(dh9) )
  {
    std::cout<<"Error: the closed volume is not a combinatorial hexahedron.\n";
    assert(false);
    return false;
  }

  return true;
}

#endif // CGAL_GENERALIZED_MAP_2_TEST
