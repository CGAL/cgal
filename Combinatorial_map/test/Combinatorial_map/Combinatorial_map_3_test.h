// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_COMBINATORIAL_MAP_3_TEST
#define CGAL_COMBINATORIAL_MAP_3_TEST 1

#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include "Combinatorial_map_test_iterators.h"

#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

template<typename CMap>
bool check_number_of_cells_3(CMap& cmap, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbvol,
                             unsigned int nbcc)
{
  if ( !cmap.is_valid() )
    {
      std::cout<<"ERROR: the cmap is not valid."<<std::endl;
      assert(false);
      return false;
    }

  std::vector<unsigned int> nbc;
  nbc=cmap.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbvol!=nbc[3] ||
      nbcc!=nbc[4])
    {
      std::cout<<"ERROR: the number of cells is not correct. We must have "
               <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbvol
               <<", "<<nbcc<<") and we have"
               <<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "<<nbc[3]<<", "
               <<nbc[4]<<")."
               <<std::endl;
      assert(false);
      return false;
    }

  return true;
}

template<typename CMap>
struct Count_nb_darts
{
  Count_nb_darts(std::size_t &nb) : m_nb(nb)
  {}

  void operator() (const typename CMap::Dart&)
  { ++m_nb; }
protected:
  std::size_t& m_nb;
};

template<class Map, class Functor>
void createAllBasicCases1()
{
  Map map;
  Functor functor;

  typename Map::Dart_handle dh, dh2, dh3;

  // 1) isolated dart
  dh = map.create_dart();
  functor(map,dh);
  map.clear();

  // 2) loop
  dh = map.create_dart();
  map.template sew<1>(dh, dh);
  functor(map,dh);
  map.clear();

  // 3) 2 darts v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh, dh2);
  functor(map,dh);
  map.clear();

  // 4) 2 darts v2
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh, dh2);
  functor(map,dh2);
  map.clear();

  // 5) 2 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh);
  functor(map,dh);
  map.clear();

  // 6) 2 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh);
  functor(map,dh2);
  map.clear();

  // 7) 3 darts v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  functor(map,dh);
  map.clear();

  // 8) 3 darts v2
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  functor(map,dh2);
  map.clear();

  // 9) 3 darts v3
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  functor(map,dh3);
  map.clear();

  // 10) 3 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  map.template sew<1>(dh3, dh);
  functor(map,dh);
  map.clear();

  // 11) 3 darts cycle v2
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  map.template sew<1>(dh3, dh);
  functor(map,dh2);
  map.clear();

  // 12) 3 darts cycle v3
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  map.template sew<1>(dh3, dh);
  functor(map,dh3);
  map.clear();

 }


template<class Map, int close>
struct Myclose
{
  static void run(Map& map)
  { map.template close<close>(); }
};

template<class Map>
struct Myclose<Map,-1>
{
  static void run(Map&)
  {}
};


template<class Map, class Functor, int close1, int close2>
void createAllBasicCases2()
{
  Map map;
  Functor functor;

  typename Map::Dart_handle dh, dh2, dh3, dh4;

  // 1) isolated dart
  dh  = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<2>(dh,dh2);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh);
  map.clear();

  // 2) loop
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh,dh);
  map.template sew<2>(dh,dh2);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh);
  map.clear();

  // 3) 2 darts v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<2>(dh2,dh3);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh);
  map.clear();

  // 4) 2 darts v2
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<2>(dh2,dh3);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh2);
  map.clear();

  // 5) 2 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh);
  map.template sew<2>(dh2,dh3);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh);
  map.clear();

  // 6) 2 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh);
  dh3 = map.create_dart();
  map.template sew<2>(dh2,dh3);
  dh3 = map.create_dart();
  map.template sew<2>(dh,dh3);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh2);
  map.clear();

  // 7) 3 darts v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  dh4 = map.create_dart();
  map.template sew<2>(dh2,dh4);
  dh4 = map.create_dart();
  map.template sew<2>(dh3,dh4);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh);
  map.clear();

  // 8) 3 darts v2
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  dh4 = map.create_dart();
  map.template sew<2>(dh2,dh4);
  dh4 = map.create_dart();
  map.template sew<2>(dh3,dh4);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh2);
  map.clear();

  // 9) 3 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh3);
  map.template sew<1>(dh3, dh);
  dh4 = map.create_dart();
  map.template sew<2>(dh,dh4);
  Myclose<Map, close1>::run(map);
  Myclose<Map, close2>::run(map);
  functor(map,dh);
  map.clear();
 }

template<class Map>
struct InsertVertex
{
  InsertVertex() : nb(1)
  {};
  void operator() (Map& map, typename Map::Dart_handle d)
  {
    std::cout<<"Before: ";
    map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
    std::cout<<"After insert_vertex "<<nb++<<" : "<<std::flush;
    map.insert_cell_0_in_cell_1(d);
    map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  }
private:
  unsigned int nb;
};

template<class Map>
bool test3D()
{
  typedef typename Map::Dart_handle Dart_handle;
  Dart_handle d,dh,dh2,d1,d2,d3,d4;
  typename Map::size_type mark;
  unsigned int nbc,nb2;
  Map map;
  cout << "***************************** TEST BASIC CREATION 3D:"
       << endl;

  dh = map.create_dart();
  dh2 = map.create_dart();

  map.template sew<1>(dh, dh2);

  if (map.is_valid()) cout << "Map valid." << endl;
  cout << "Nombre de brins : " << map.number_of_darts() << endl;
  map.clear();

  cout << "***************************** TEST BASIC CREATION 3D DONE."
       << endl;

  cout << "***************************** TEST CREATION WITH EMBEDDING 3D:"
       << endl;

  dh = map.create_dart();
  dh2 = map.create_dart();

  cout << "Parcours all : "; ;

  map.template sew<1>(dh, dh2);
  cout << "Parcours all : "; ;
  cout << endl << endl;

  map.make_edge();

  dh = map.make_combinatorial_tetrahedron();
  cout << "Nombre de brins : " << map.number_of_darts() << endl;

  cout << "Parcours de CC : ";
  for (CGAL::CMap_dart_iterator_of_orbit<Map,1,2,3> it(map, dh);it.cont();++it)
    {
      if (it.prev_operation() != CGAL::OP_BETAI)  cout << "\nNew facet : ";
      //    cout << &*it << ", ";
    }
  cout << endl;

  cout << "Parcours all : "; ;

  if (map.is_valid()) cout << "Map valid." << endl;
  map.clear();

  cout << "***************************** TEST CREATION WITH EMBEDDING 3D DONE."
       << endl;

  cout << "***************************** TEST SEW 3D:" << endl;

  map.make_edge();

  d1 = map.make_combinatorial_tetrahedron();

  d2 = map.make_combinatorial_tetrahedron();

  cout << "Parcours all : "; ;

  map.template sew<3>(d1, d2);

  cout << "Parcours all : "; ;

  if (map.is_valid()) cout << "Map valid." << endl;
  map.clear();

  cout << "***************************** TEST SEW 3D DONE." << endl;

  cout << "***************************** TEST TRIANGULATION_2 3D:" << endl;

  d1 = map.make_combinatorial_tetrahedron();

  cout << "Before insert_cell_0_in_cell_2: "; ;

  map.insert_cell_0_in_cell_2(d1);

  cout << "After insert_cell_0_in_cell_2: "; ;

  if (map.is_valid()) cout << "Map valid." << endl;

  map.clear();

  cout << "***************************** TEST TRIANGULATION_2 3D DONE."
       << endl;

  cout << "***************************** TEST TRIANGULATION_3 3D:" << endl;

  // Create 2 tetrahedra
  d1 = map.make_combinatorial_tetrahedron();

  d2 = map.make_combinatorial_tetrahedron();

  // Sew the 2 tetrahedra along one facet
  map.template sew<3>(d1, d2);

  cout << "map valid: " << map.is_valid() << endl;
  // Display all the vertices
  cout << "Before insert_cell_0_in_cell_2: "; ;

  // Triangulate the facet between the two tetrahedra
  map.insert_cell_0_in_cell_2(d1);

  cout << "map valid: " << map.is_valid() << endl;
  cout << "After insert_cell_0_in_cell_2: "; ;

  if (map.is_valid()) cout << "Map valid." << endl;

  map.template unsew<3>(d1);
  cout << "After unsew: "; ;

  if (map.is_valid()) cout << "Map valid." << endl;
  map.clear();

  cout << "***************************** TEST TRIANGULATION_3 3D DONE."
       << endl;


  cout << "***************************** TEST ITERATORS 3D:" << endl;

  for (int i = 0; i < 1000; ++i)
    {
       d1 = map.make_combinatorial_tetrahedron();
       d2 = map.make_combinatorial_tetrahedron();
       d3 = map.make_combinatorial_tetrahedron();
      // Sew the 2 tetrahedra along one facet
      map.template sew<3>(d1, d2);
      map.template sew<3>(map.beta(d2,2), d3);
      map.template sew<3>(map.beta(d1,2), map.beta(d3,2));
    }

  if ( !test_iterators_3(map) )
  {  assert(false); return false; }

  // Two nested iterators
  cout << "Nombre de brins : " << map.number_of_darts() << ", "
       << "Nombre de CC : " << flush;
  mark = map.get_new_mark();
  nbc = 0;
  for (typename Map::Dart_range::const_iterator it1(map.darts().begin());
       it1!=map.darts().end(); ++it1)
    {
      if (!map.is_marked(it1, mark))
        {
          ++nbc;
          for (typename Map::template Dart_of_orbit_range<1,2,3>::const_iterator it2(map, it1);
               it2.cont(); ++it2)
            { map.mark(it2, mark); }
        }
    }
  cout << nbc << endl;
  cout << "All the darts marked ? " << map.is_whole_map_marked(mark)
       << endl;
  map.unmark_all(mark);

  // Iterator stl like
  {
    nbc = 0, nb2 = 0;
    unsigned int nbtest=0;
    cout << "Iterator stl like: #cellules=" << flush;
   for (typename Map::Dart_range::const_iterator it1(map.darts().begin());
        it1!=map.darts().end(); ++it1)
      {
        ++nb2;
        if (!map.is_marked(it1, mark))
          {
            ++nbc;
            for (typename Map::template Dart_of_orbit_range<2>::const_iterator
                   it2(map.template darts_of_orbit<2>(it1).begin());
                 it2 != map.template darts_of_orbit<2>(it1).end(); ++it2)
              { map.mark(it2, mark); }
            for (typename Map::template Dart_of_cell_range<2>::const_iterator
                   it2(map.template darts_of_cell<2>(it1).begin());
                 it2 != map.template darts_of_cell<2>(it1).end(); ++it2)
              { ++nbtest; }
            for (typename Map::template Dart_of_involution_range<2>::const_iterator
                   it2(map.template darts_of_involution<2>(it1).begin());
                 it2 != map.template darts_of_involution<2>(it1).end(); ++it2)
              { ++nbtest; }
            for (typename Map::template One_dart_per_incident_cell_range<2,0>::const_iterator
                   it2(map.template one_dart_per_incident_cell<2,0>(it1).begin());
                 it2 != map.template one_dart_per_incident_cell<2,0>(it1).end(); ++it2)
              { ++nbtest; }
          }
      }
    cout << nbc << "." << ", #brins=" << nb2 << "." << endl
         << "All the darts marked ? " << map.is_whole_map_marked(mark) << endl;
    {
      for (typename Map::template One_dart_per_cell_range<0>::const_iterator
             it2(map.template one_dart_per_cell<0>().begin());
           it2 != map.template one_dart_per_cell<0>().end(); ++it2)
      { ++nbtest; }
    }
    cout<<"Different const_iterators: "<<nbtest<<std::endl;
    map.unmark_all(mark);
  }
  map.free_mark(mark);

  {
    cout<<"Test operator= between iterators."<<std::endl;
    typename Map::Dart_range::const_iterator it1(map.darts().begin());
    it1=map.darts().begin();

    typename Map::template Dart_of_orbit_range<2>::const_iterator
      it2(map.template darts_of_orbit<2>(it1).begin());
    it2 = map.template darts_of_orbit<2>(it1).begin();

    typename Map::template Dart_of_cell_range<2>::const_iterator
      it3(map.template darts_of_cell<2>(it1).begin());
    it3 = map.template darts_of_cell<2>(it1).begin();

    typename Map::template Dart_of_involution_range<2>::const_iterator
      it4(map.template darts_of_involution<2>(it1).begin());
    it4 = map.template darts_of_involution<2>(it1).begin();

    typename Map::template One_dart_per_incident_cell_range<2,0>::const_iterator
      it5(map.template one_dart_per_incident_cell<2,0>(it1).begin());
    it5 = map.template one_dart_per_incident_cell<2,0>(it1).begin();

    typename Map::template One_dart_per_cell_range<0>::const_iterator
      it6(map.template one_dart_per_cell<0>().begin());
    it6 = map.template one_dart_per_cell<0>().begin();
  }

  {
    std::size_t nbtest=0;
    cout<<"std::for_each iterators : ";
    std::for_each(map.darts().begin(), map.darts().end(),
                  Count_nb_darts<Map>(nbtest));
    typename Map::Dart_range::const_iterator it1(map.darts().begin());
    std::for_each(map.template darts_of_orbit<1,2>(it1).begin(),
                  map.template darts_of_orbit<1,2>(it1).end(),
                  Count_nb_darts<Map>(nbtest));
    std::for_each(map.template darts_of_cell<3>(it1).begin(),
                  map.template darts_of_cell<3>(it1).end(),
                  Count_nb_darts<Map>(nbtest));
    std::for_each(map.template darts_of_involution<2>(it1).begin(),
                  map.template darts_of_involution<2>(it1).end(),
                  Count_nb_darts<Map>(nbtest));
    std::for_each(map.template one_dart_per_incident_cell<0,2>(it1).begin(),
                  map.template one_dart_per_incident_cell<0,2>(it1).end(),
                  Count_nb_darts<Map>(nbtest));
    std::for_each(map.template one_dart_per_cell<1>().begin(),
                  map.template one_dart_per_cell<1>().end(),
                  Count_nb_darts<Map>(nbtest));
    std::cout<<nbtest<<", comparison using size(): ";

    std::size_t nbtest2 = map.darts().size()+
      map.template darts_of_orbit<1,2>(it1).size()+
      map.template darts_of_cell<3>(it1).size()+
      map.template darts_of_involution<2>(it1).size()+
      map.template one_dart_per_incident_cell<0,2>(it1).size()+
      map.template one_dart_per_cell<1>().size();
    std::cout<<nbtest2<<(nbtest==nbtest2?". OK":". PROBLEM")<<std::endl;
  }

  test_iterators_3(map);

  map.clear();

  cout << "***************************** TEST ITERATORS 3D DONE." << endl;

  cout << "***************************** TEST INCIDENCE ITERATORS 3D:"
       << endl;

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.make_combinatorial_hexahedron();
  map.template sew<3>(d1, d2);
  cout << "Map valid : " << map.is_valid() << endl;

  mark = map.get_new_mark();

  map.free_mark(mark);
  map.clear();

  cout << "***************************** TEST INCIDENCE ITERATORS 3D DONE."
       << endl;

  cout << "***************************** TEST VERTEX REMOVAL 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex1: " << flush; map.template remove_cell<0>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex2: " << flush; map.template remove_cell<0>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_edge();
  d2 = map.beta(d1,2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex3: " << flush; map.template remove_cell<0>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex4: " << flush; map.template remove_cell<0>(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_edge();
  map.template sew<1>(d1, d1);
  map.template sew<1>(map.beta(d1,2), map.beta(d1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex5: " << flush; map.template remove_cell<0>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.beta(d1,0); d3 = map.beta(d1,1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex6: " << flush; map.template remove_cell<0>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex7: " << flush; map.template remove_cell<0>(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex8: " << flush; map.template remove_cell<0>(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<3>(d1, d2); d2 = map.beta(d1,0); d3 = map.beta(d1,1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex9: " << flush; map.template remove_cell<0>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex10: " << flush; map.template remove_cell<0>(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex11: " << flush; map.template remove_cell<0>(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  cout << "***************************** TEST VERTEX REMOVAL 3D DONE."
       << endl;

  cout << "***************************** TEST EDGE REMOVAL 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge1: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge2: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_edge();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge3: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_edge();
  map.template sew<1>(d1, d1);
  map.template sew<1>(map.beta(d1,2), map.beta(d1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge4: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_edge();
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge5: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.beta(d1,0); d3 = map.beta(d1,1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge6: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge7: " << flush; map.template remove_cell<1>(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge8: " << flush; map.template remove_cell<1>(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<3>(d1, d2); d2 = map.beta(d1,0); d3 = map.beta(d1,1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge9: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge10: " << flush; map.template remove_cell<1>(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge11: " << flush; map.template remove_cell<1>(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge12: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge13: " << flush;
  map.template remove_cell<1>(map.beta(d1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge14: " << flush; map.template remove_cell<1>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.make_combinatorial_hexahedron();
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d3 = map.insert_cell_0_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  {
    std::vector<Dart_handle> V;
    {
      for ( typename Map::template Dart_of_cell_range<0,2>::iterator
              it =  map.template darts_of_cell<0,2>(d3).begin();
            it!=map.template darts_of_cell<0,2>(d3).end(); ++it)
        V.push_back(it);
    }

    typedef typename std::vector<Dart_handle>::iterator vector_it;
    for ( vector_it it=V.begin(); it!=V.end(); ++it)
      {
        cout << "remove edge15: " << flush; map.template remove_cell<1>(*it);
        map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
      }
  }
  map.clear();

  cout << "***************************** TEST EDGE REMOVAL 3D DONE."
       << endl;

  cout << "***************************** TEST FACET REMOVAL 3D:"
       << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet1: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet2: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet3: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_hexahedron();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet4: " << flush;
  map.template remove_cell<2>(map.beta(d1,2,1,1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet5: " << flush;
  map.template remove_cell<2>(map.beta(d1,1,1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet6: " << flush;
  map.template remove_cell<2>(map.beta(d1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet7: " << flush;
  map.template remove_cell<2>(map.beta(d1,1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet8: " << flush;
  map.template remove_cell<2>(map.beta(d1,0,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet9: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.make_combinatorial_hexahedron();
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  d3 = map.beta(d1,1,2);

  cout << "remove facet10: " << flush;
  map.template remove_cell<2>(map.beta(d1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet11: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.beta(d3,0,2);
  cout << "remove edge12: " << flush; map.template remove_cell<1>(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.make_combinatorial_hexahedron();
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d3 = map.insert_cell_0_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  {
    std::vector<Dart_handle> V;
    {
      for ( typename Map::template Dart_of_cell_range<0,2>::iterator it =
              map.template darts_of_cell<0,2>(d3).begin();
            it!=map.template darts_of_cell<0,2>(d3).end(); ++it)
        V.push_back(it);
    }

    for (typename std::vector<Dart_handle>::iterator it=V.begin(); it!=V.end(); ++it)
      {
        cout << "remove facet13: " << flush; map.template remove_cell<2>(*it);
        map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
      }
  }
  std::cout<<"Size of map="<<map.bytes()<<std::endl;
  map.clear();

  cout << "***************************** TEST FACET REMOVAL 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT VERTEX 3D:"
       << endl;


  std::cout<<"************************* Test createAllBasicCases1 *************************"<<std::endl;
  createAllBasicCases1<Map,InsertVertex<Map> >();

  std::cout<<"************************* Test createAllBasicCases2 *************************"<<std::endl;
  createAllBasicCases2<Map,InsertVertex<Map>,-1,-1>();

  std::cout<<"************************* Test createAllBasicCases2<3> *************************"<<std::endl;
  createAllBasicCases2<Map,InsertVertex<Map>,-1,3>();

  std::cout<<"************************* Test different cases *************************"<<std::endl;
  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex1: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex2: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_edge();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex3: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_edge();
  map.template sew<1>(d1, d1);
  map.template sew<1>(map.beta(d1,2), map.beta(d1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex4: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_edge();
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex5: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.beta(d1,0); d3 = map.beta(d1,1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex6: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex7: " << flush; map.insert_cell_0_in_cell_1(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex8: " << flush; map.insert_cell_0_in_cell_1(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<3>(d1, d2); d2 = map.beta(d1,0); d3 = map.beta(d1,1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex9: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex10: " << flush; map.insert_cell_0_in_cell_1(d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex11: " << flush; map.insert_cell_0_in_cell_1(d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex12: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex13: " << flush;
  map.insert_cell_0_in_cell_1(map.beta(d1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex14: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.make_combinatorial_hexahedron();
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "insert vertex15: " << flush;
  map.insert_cell_0_in_cell_1(map.beta(d1,1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex16: " << flush;
  map.insert_cell_0_in_cell_1(map.beta(d1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex17: " << flush;
  map.insert_cell_0_in_cell_1(map.beta(d1,0));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex18: " << flush; map.insert_cell_0_in_cell_1(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  cout << "***************************** TEST INSERT VERTEX 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT EDGE 3D:"
       << endl;

  d1 = map.make_edge(); map.template sew<1>(d1, map.make_edge());
  d2 = map.make_edge(); map.template sew<0>(d2, map.make_edge());
  map.template sew<0>(map.template beta<0>(d2), d2);
  map.template sew<1>(map.template beta<1>(d1), d1);
  map.template sew<3>(d1, d2);
  map.template sew<3>(map.beta(d1,2), map.beta(d2,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge3: " << flush;
  map.insert_cell_1_in_cell_2(d1, map.template beta<1>(d1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(4 );
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge4: " << flush;
  map.insert_cell_1_in_cell_2(d1, map.beta(d1,1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(4);
  d2 = map.make_combinatorial_polygon(4);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge5: " << flush;
  map.insert_cell_1_in_cell_2(d1, map.beta(d1,1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(4);
  d2 = map.make_combinatorial_polygon(4);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge6: " << flush;
  map.insert_cell_1_in_cell_2(d1, map.beta(d1,1,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
  map.clear();

  cout << "***************************** TEST INSERT EDGE 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT DANGLING EDGE 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge1: " << flush;
  map.insert_dangling_cell_1_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_edge();
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge2: " << flush;
  map.insert_dangling_cell_1_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_edge();
  d2 = map.make_edge();
  map.template sew<3>(d1, d2);
  map.template sew<3>(map.beta(d1,2), map.beta(d2,2));
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge3: " << flush;
  map.insert_dangling_cell_1_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(4);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge4: " << flush;
  map.insert_dangling_cell_1_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.make_combinatorial_polygon(4);
  d2 = map.make_combinatorial_polygon(4);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge5: " << flush;
  map.insert_dangling_cell_1_in_cell_2(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  cout << "***************************** TEST INSERT DANGLING EDGE 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT FACET 3D:"
       << endl;

  std::vector<Dart_handle> v;

  d1 = map.make_combinatorial_polygon(4);
  v.push_back(d1); v.push_back(map.beta(v[0],1));
  v.push_back(map.beta(v[1],1)); v.push_back(map.beta(v[2],1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert facet1: " << flush;
  map.insert_cell_2_in_cell_3(v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();

  d1 = map.make_combinatorial_polygon(3);
  d2 = map.make_combinatorial_polygon(3);
  map.template sew<2>(d1, d2);
  v.push_back(d1); v.push_back(map.beta(v[0],1));
  v.push_back(map.beta(v[1],1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert facet2: " << flush;
  map.insert_cell_2_in_cell_3(v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.beta(d1,2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet3: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  v.push_back(d2); v.push_back(map.beta(v[0],1,2,1));
  v.push_back(map.beta(v[1],1,2,1)); v.push_back(map.beta(v[2],1,2,1));
  cout << "insert facet3: " << flush;
  map.insert_cell_2_in_cell_3(v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();

  d1 = map.make_combinatorial_hexahedron();
  d2 = map.make_combinatorial_hexahedron();
  map.template sew<3>(d1,d2);
  d3 = map.beta(d1, 2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet4: " << flush; map.template remove_cell<2>(d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  v.push_back(d3); v.push_back(map.beta(v[0],1,2,1));
  v.push_back(map.beta(v[1],1,2,1)); v.push_back(map.beta(v[2],1,2,1));
  cout << "insert facet4: " << flush; map.insert_cell_2_in_cell_3(v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  Map map2;
  d1 = map2.make_combinatorial_hexahedron();
  d2 = map2.make_combinatorial_hexahedron();
  map2.template sew<3>(d1,d2);
  if ( !map.is_isomorphic_to(map2, false) )
  {
    std::cout<<"Error: map and map2 are not isomorphic (after insertion/removal).\n";
    assert(false);
    return false;
  }

  if (CGAL::degree<Map, 0>(map2, d1)!=4)
  {
    std::cout<<"Error: 0-degree is wrong: "<<CGAL::degree<Map, 0>(map2, d1)<<" instead of 4."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::degree<Map, 1>(map2, d1)!=3)
  {
    std::cout<<"Error: 1-degree is wrong: "<<CGAL::degree<Map, 1>(map2, d1)<<" instead of 3."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::degree<Map, 2>(map2, d1)!=2)
  {
    std::cout<<"Error: 2-degree is wrong: "<<CGAL::degree<Map, 2>(map2, d1)<<" instead of 2."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::codegree<Map, 1>(map2, d1)!=2)
  {
    std::cout<<"Error: 1-codegree is wrong: "<<CGAL::codegree<Map, 1>(map2, d1)<<" instead of 2."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::codegree<Map, 2>(map2, d1)!=4)
  {
    std::cout<<"Error: 2-codegree is wrong: "<<CGAL::codegree<Map, 2>(map2, d1)<<" instead of 4."<<std::endl;
    assert(false);
    return false;
  }

  if (CGAL::codegree<Map, 3>(map2, d1)!=6)
  {
    std::cout<<"Error: 3-codegree is wrong: "<<CGAL::codegree<Map, 3>(map2, d1)<<" instead of 6."<<std::endl;
    assert(false);
    return false;
  }

  map.clear(); v.clear();

  cout << "***************************** TEST INSERT FACET 3D DONE."
       << endl;

  return true;
}

#endif // CGAL_COMBINATORIAL_MAP_3_TEST
