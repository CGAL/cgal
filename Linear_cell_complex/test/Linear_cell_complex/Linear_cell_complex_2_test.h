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
#ifndef CGAL_LCC_2_TEST_H
#define CGAL_LCC_2_TEST_H

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Triangulation_2_to_lcc.h>
#include <CGAL/Delaunay_triangulation_2.h>
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

void trace_display_msg(const char* msg)
{
  CGAL_USE(msg);
#ifdef LCC_TRACE_TEST_BEGIN
  std::cout<<"***************** "<<msg<<"***************** "<<std::endl;
#endif
}

template<typename Map, int i, typename Info=
         typename Map::template Attribute_type<i>::type::Info>
struct SetInfoIfNonVoid
{
  static void run(Map& map,
                  typename Map::template Attribute_descriptor<i>::type attr,
                  long long int nb)
  {
    map.template info_of_attribute<i>(attr)=
      typename Map::template Attribute_type<i>::type::Info(nb);
  }
};
template<typename Map, int i>
struct SetInfoIfNonVoid<Map, i, void>
{
  static void run(Map&, typename Map::template Attribute_descriptor<i>::type,
                  long long int)
  {}
};

template<typename Map, unsigned int i, typename Attr=typename Map::
         template Attribute_type<i>::type>
struct CreateAttributes
{
  static void run(Map& map)
  {
    long long int nb=0;
    for(typename Map::Dart_range::iterator it=map.darts().begin(),
        itend=map.darts().end(); it!=itend; ++it)
    {
      if ( map.template attribute<i>(it)==map.null_descriptor )
      {
        map.template set_attribute<i>(it, map.template create_attribute<i>());
        SetInfoIfNonVoid<Map, i>::run(map, map.template attribute<i>(it), ++nb);
      }
    }
  }
};

template<typename Map, typename Attr>
struct CreateAttributes<Map, 0, Attr>
{
  static void run(Map& amap)
  {
    long long int nb=0;
    for ( typename Map::template Attribute_range<0>::type::iterator
          it=amap.template attributes<0>().begin(),
          itend=amap.template attributes<0>().end(); it!=itend; ++it )
      SetInfoIfNonVoid<Map, 0>::run(amap, it, ++nb);
  }
};

template<typename Map, unsigned int i>
struct CreateAttributes<Map, i, CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map>
struct CreateAttributes<Map, 0, CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map, typename Info=typename Map::Dart_info>
struct InitDartInfo
{
  static void run(Map& map)
  {
    long long int nb=0;
    for(typename Map::Dart_range::iterator it=map.darts().begin(),
        itend=map.darts().end(); it!=itend; ++it)
    {
      nb=CGAL::get_default_random().get_int(0,20000);
      map.info(it)=Info(nb);
    }
  }
};

template<typename Map>
struct InitDartInfo<Map, CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map>
void create_attributes_2(Map& map)
{
  CreateAttributes<Map, 0>::run(map);
  CreateAttributes<Map, 1>::run(map);
  CreateAttributes<Map, 2>::run(map);
  InitDartInfo<Map>::run(map);
}

// Test orientation specialized below only for CMap. For GMap return true.
template<typename LCC, typename Map=typename LCC::Combinatorial_data_structure>
struct Test_change_orientation_LCC_2
{
  static bool run()
  { return true; }
};

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
    std::cout << " dart " << lcc.darts().index(it) << "; beta[i]=";
    /* for ( unsigned int i=0; i<=LCC::dimension; ++i)
    {
      std::cout << &(*it->beta(i)) << ",\t";
      if (lcc.is_free(it, i)) std::cout << "\t";
      } */
    std::cout<<lcc.point(it);
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

  typedef typename LCC::Dart_descriptor Dart_descriptor;
  typedef typename LCC::Point Point;

  // Construction operations
  trace_test_begin();
  Dart_descriptor dh1=lcc.make_segment(Point(0,0),Point(1,0), true);
  Dart_descriptor dh2=lcc.make_segment(Point(2,0),Point(2,1), true);
  Dart_descriptor dh3=lcc.make_segment(Point(2,2),Point(3,1), true);
  if ( !check_number_of_cells_2(lcc, 6, 3, 6, 3) )
    return false;

  { // Test swap operator
    LCC lcc2;
    lcc2.swap(lcc);
    if ( !check_number_of_cells_2(lcc, 0, 0, 0, 0) )
      return false;
    if ( !check_number_of_cells_2(lcc2, 6, 3, 6, 3) )
      return false;

    lcc.swap(lcc2);
    if ( !check_number_of_cells_2(lcc2, 0, 0, 0, 0) )
      return false;
    if ( !check_number_of_cells_2(lcc, 6, 3, 6, 3) )
      return false;

    // And test operator=
    LCC lcc3;
    lcc3=lcc;
    if ( !check_number_of_cells_2(lcc3, 6, 3, 6, 3) )
      return false;
    if (!lcc.is_isomorphic_to(lcc3))
      return false;
  }

  typename LCC::Vertex_attribute_descriptor vh=lcc.template attribute<0>(dh1);
  if (!lcc.template is_attribute_used<0>(vh)) return false;

  trace_test_begin();
  lcc.template sew<1>(dh1, dh2);
  lcc.template sew<1>(lcc.other_orientation(dh2), dh3);
  if ( !check_number_of_cells_2(lcc, 4, 3, 4, 1) )
    return false;

  trace_test_begin();
  Dart_descriptor dh5=lcc.make_triangle(Point(5,5),Point(7,5),Point(6,6));
  Dart_descriptor dh6=lcc.make_triangle(Point(5,4),Point(7,4),Point(6,3));
  if ( !check_number_of_cells_2(lcc, 10, 9, 6, 3) )
    return false;

  trace_test_begin();
  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_2(lcc, 8, 8, 6, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_2(lcc, 9, 9, 6, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_2(lcc, 10, 12, 8, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0));
  if ( !check_number_of_cells_2(lcc, 11, 13, 8, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5));
  if ( !check_number_of_cells_2(lcc, 12, 16, 10, 2) )
    return false;

  trace_test_begin();
  Dart_descriptor dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2));
  if ( !check_number_of_cells_2(lcc, 13, 17, 10, 2) )
    return false;

  // Removal operations
  trace_test_begin();
  lcc.template remove_cell<1>(dh11);
  if ( !check_number_of_cells_2(lcc, 12, 16, 10, 2) )
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

  for ( typename std::vector<Dart_descriptor>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    if (lcc.is_dart_used(*it)) // For GMap because we have 2 dart per edge incident to the vertex
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
  lcc.template unsew<1>(dh1);
  if ( !check_number_of_cells_2(lcc, 5, 3, 5, 2) )
    return false;

  trace_test_begin();
  lcc.template unsew<1>(lcc.other_orientation(dh2));
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
    CGAL::import_from_plane_graph<LCC>(lcc,in);
    if ( !check_number_of_cells_2(lcc, 66, 166, 104, 2) )
      return false;
    lcc.clear();
  }

  // Construction from Triangulation_2
  {
    trace_test_begin();
    CGAL::Triangulation_2<typename LCC::Traits> T;
    std::ifstream in("data/points2D.txt");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/points2D.txt'"<<std::endl;
      return false;
    }
    T.insert ( std::istream_iterator < Point >(in),
               std::istream_iterator < Point >() );
    CGAL::import_from_triangulation_2<LCC>(lcc,T);
    if ( !lcc.is_valid() )
      return false;

    // Pb: the triangulation_2 is not the same on different machines ?
    if ( !check_number_of_cells_2(lcc, 501, 1497, 998, 1) )
      return false;

    lcc.clear();
    trace_test_end();
  }

  if ( !Test_change_orientation_LCC_2<LCC>::run() )
    return false;

  return true;
}

template<typename LCC>
struct Test_change_orientation_LCC_2<LCC, CGAL::Combinatorial_map_tag>
{
  static bool run()
  {
    LCC lcc;

    std::ifstream in("data/graph.txt");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/graph.txt'"<<std::endl;
      return false;
    }
    CGAL::import_from_plane_graph<LCC>(lcc,in);

    trace_test_begin();

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
    if ( !lcc2.is_isomorphic_to(lcc, false, false, false) )
    { assert(false); return false; }
    trace_test_end();

    trace_test_begin();
    lcc.reverse_orientation();
    if ( !lcc.is_valid() ) { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc, false, false, false) )
    { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc) )
    { assert(false); return false; }
    trace_test_end();

    trace_test_begin();
    lcc.reverse_orientation_connected_component(lcc.darts().begin());
    if ( !lcc.is_valid() ) { assert(false); return false; }
    if ( lcc2.is_isomorphic_to(lcc) )
    { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc, false, false, false) )
    { assert(false); return false; }
    trace_test_end();

    trace_test_begin();
    lcc.reverse_orientation_connected_component(lcc.darts().begin());
    if ( !lcc.is_valid() ) { assert(false); return false; }
    if ( !lcc2.is_isomorphic_to(lcc) )
    { assert(false); return false; }
    trace_test_end();

    return true;
  }
};

#endif // CGAL_LCC_2_TEST_H
