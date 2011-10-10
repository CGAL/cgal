// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#ifndef CGAL_COMBINATORIAL_MAP_3_TEST
#define CGAL_COMBINATORIAL_MAP_3_TEST 1

#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>

#include <iostream>
#include <fstream>

using namespace std;

template<class Map>
void drawCell3(Map& amap, typename Map::Dart_handle adart, int aorbit, int mark)
{
  cout << "Orbite " << Map::ORBIT_NAME[aorbit] << " (";
  /*  CGAL::CMap_dart_iterator_basic_of_orbit_3<Map> it2(amap, adart, aorbit, mark);
  for (;it2.cont(); ++it2)
    {
      cout << (*it2)->vertex()->point() << ", ";
      }*/
  cout << ")" << flush;
}

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

template<class Map, class Functor>
void createAllBasicCases2(int close1, int close2)
{
  Map map; 
  Functor functor;

  typename Map::Dart_handle dh, dh2, dh3, dh4;

  // 1) isolated dart
  dh  = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<2>(dh,dh2);
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
  functor(map,dh);
  map.clear();

  // 2) loop
  dh = map.create_dart();
  dh2 = map.create_dart();
  map.template sew<1>(dh,dh);
  map.template sew<2>(dh,dh2);
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
  functor(map,dh);
  map.clear();

  // 3) 2 darts v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<2>(dh2,dh3);
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
  functor(map,dh);
  map.clear();

  // 4) 2 darts v2
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<2>(dh2,dh3);
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
  functor(map,dh2);
  map.clear();

  // 5) 2 darts cycle v1
  dh = map.create_dart();
  dh2 = map.create_dart();
  dh3 = map.create_dart();
  map.template sew<1>(dh, dh2);
  map.template sew<1>(dh2, dh);
  map.template sew<2>(dh2,dh3);
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
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
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
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
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
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
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
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
  if (close1!=-1) map.close(close1);
  if (close2!=-1) map.close(close2);
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
    insert_cell_0_in_cell_1(map,d);
    map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  }
private:
  unsigned int nb;
};

template<class Map>
  void test3D()
{
  typedef typename Map::Dart_handle Dart_handle;
  Dart_handle d,dh,dh2,d1,d2,d3,d4;
  int mark;
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

  make_edge(map );

  dh = make_combinatorial_tetrahedron(map);
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

  make_edge(map);

  d1 = make_combinatorial_tetrahedron(map);

  d2 = make_combinatorial_tetrahedron(map);

  cout << "Parcours all : "; ;

  map.template sew<3>(d1, d2);

  cout << "Parcours all : "; ;

  if (map.is_valid()) cout << "Map valid." << endl;
  map.clear();
    
  cout << "***************************** TEST SEW 3D DONE." << endl;

  cout << "***************************** TEST TRIANGULATION_2 3D:" << endl;

  d1 = make_combinatorial_tetrahedron(map);

  cout << "Before insert_cell_0_in_cell_2: "; ;

  insert_cell_0_in_cell_2(map,d1);

  cout << "After insert_cell_0_in_cell_2: "; ;

  if (map.is_valid()) cout << "Map valid." << endl;
  map.clear();
    
  cout << "***************************** TEST TRIANGULATION_2 3D DONE."
       << endl;

  cout << "***************************** TEST TRIANGULATION_3 3D:" << endl;

  // Create 2 tetrahedra
  d1 = make_combinatorial_tetrahedron(map);

  d2 = make_combinatorial_tetrahedron(map);

  // Sew the 2 tetrahedra along one facet
  map.template sew<3>(d1, d2);

  cout << "map valid: " << map.is_valid() << endl;
  // Display all the vertices
  cout << "Before insert_cell_0_in_cell_2: "; ;

  // Triangulate the facet between the two tetrahedra
  insert_cell_0_in_cell_2(map,d1);

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
       d1 = make_combinatorial_tetrahedron(map);
       d2 = make_combinatorial_tetrahedron(map);
       d3 = make_combinatorial_tetrahedron(map);
      // Sew the 2 tetrahedra along one facet
      map.template sew<3>(d1, d2);
      map.template sew<3>(d2->beta(2), d3);
      map.template sew<3>(d1->beta(2), d3->beta(2));
    }

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

  // Tout les parcours possibles :
  /*  for (int i = CGAL::SELF_ORBIT; i <= CGAL::ALL_DARTS_ORBIT; ++i)
    {
      cout << "Parcours orbite " << Map::ORBIT_NAME[i] << " : #cellules=" << flush;
      unsigned int nbc = 0, nb2 = 0;
      for (CGAL::CMap_dart_iterator_of_all<Map> it1(map);it1.cont(); ++it1)
        {
	  ++nb2;
	  if (!map.is_marked(*it1, mark))
            {
	      ++nbc;
	      CGAL::CMap_dart_iterator_basic_of_orbit_3<Map> it2(map, *it1, i, mark);
	      for (;it2.cont(); ++it2)
		{}
            }
        }
      cout << nbc << "." << ", #brins=" << nb2 << "." << endl
	   << "All the darts marked ? " << map.is_whole_map_marked(mark)
	   << endl;
	   map.unmark_all(mark);
    }*/

  // Iterator stl like
  {
    nbc = 0, nb2 = 0;
	unsigned nbtest=0;
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
  map.clear();
    
  cout << "***************************** TEST ITERATORS 3D DONE." << endl;

  cout << "***************************** TEST INCIDENCE ITERATORS 3D:"
       << endl;

  d1 = make_combinatorial_hexahedron(map);
  d2 = make_combinatorial_hexahedron(map);
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
  cout << "remove vertex1: " << flush; CGAL::remove_cell<Map,0>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex2: " << flush; CGAL::remove_cell<Map,0>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map);
  d2 = d1->beta(2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex3: " << flush; CGAL::remove_cell<Map,0>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex4: " << flush; CGAL::remove_cell<Map,0>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map);
  map.template sew<1>(d1, d1); map.template sew<1>(d1->beta(2), d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove vertex5: " << flush; CGAL::remove_cell<Map,0>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex6: " << flush; CGAL::remove_cell<Map,0>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex7: " << flush; CGAL::remove_cell<Map,0>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex8: " << flush; CGAL::remove_cell<Map,0>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<3>(d1, d2); d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex9: " << flush; CGAL::remove_cell<Map,0>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex10: " << flush; CGAL::remove_cell<Map,0>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove vertex11: " << flush; CGAL::remove_cell<Map,0>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
    
  cout << "***************************** TEST VERTEX REMOVAL 3D DONE."
       << endl;

  cout << "***************************** TEST EDGE REMOVAL 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge1: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge2: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge3: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map);
  map.template sew<1>(d1, d1); map.template sew<1>(d1->beta(2), d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge4: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map);
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge5: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge6: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge7: " << flush; CGAL::remove_cell<Map,1>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge8: " << flush; CGAL::remove_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<3>(d1, d2); d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge9: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge10: " << flush; CGAL::remove_cell<Map,1>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove edge11: " << flush; CGAL::remove_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge12: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge13: " << flush; CGAL::remove_cell<Map,1>(map,d1->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove edge14: " << flush; CGAL::remove_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_hexahedron(map);
  d2 = make_combinatorial_hexahedron(map);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d3 = insert_cell_0_in_cell_2(map,d1);
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
	cout << "remove edge15: " << flush; CGAL::remove_cell<Map,1>(map,*it);
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
  cout << "remove facet1: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet2: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet3: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_hexahedron(map);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet4: " << flush; CGAL::remove_cell<Map,2>(map,d1->beta(2)->beta(1)->beta(1)->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet5: " << flush; CGAL::remove_cell<Map,2>(map,d1->beta(1)->beta(1)->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet6: " << flush; CGAL::remove_cell<Map,2>(map,d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet7: " << flush; CGAL::remove_cell<Map,2>(map,d1->beta(1)->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet8: " << flush; CGAL::remove_cell<Map,2>(map,d1->beta(0)->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet9: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_hexahedron(map);
  d2 = make_combinatorial_hexahedron(map);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  d3 = d1->beta(1)->beta(2);

  cout << "remove facet10: " << flush; CGAL::remove_cell<Map,2>(map,d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "remove facet11: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = d3->beta(0)->beta(2);
  cout << "remove edge12: " << flush; CGAL::remove_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_hexahedron(map);
  d2 = make_combinatorial_hexahedron(map);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d3 = insert_cell_0_in_cell_2(map,d1);
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
	cout << "remove facet13: " << flush; CGAL::remove_cell<Map,2>(map,*it);
        map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
      }
  }
  map.clear();
    
  cout << "***************************** TEST FACET REMOVAL 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT VERTEX 3D:"
       << endl;


  std::cout<<"************************* Test createAllBasicCases1 *************************"<<std::endl;
  createAllBasicCases1<Map,InsertVertex<Map> >();

  std::cout<<"************************* Test createAllBasicCases2 *************************"<<std::endl;
  createAllBasicCases2<Map,InsertVertex<Map> >(-1,-1);

  std::cout<<"************************* Test createAllBasicCases2<3> *************************"<<std::endl;
  createAllBasicCases2<Map,InsertVertex<Map> >(-1,3);

  std::cout<<"************************* Test different cases *************************"<<std::endl;
  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex1: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex2: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_edge(map);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex3: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_edge(map);
  map.template sew<1>(d1, d1); map.template sew<1>(d1->beta(2), d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex4: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_edge(map);
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex5: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,3);
  d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex6: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex7: " << flush; insert_cell_0_in_cell_1(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex8: " << flush; insert_cell_0_in_cell_1(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<3>(d1, d2); d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex9: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex10: " << flush; insert_cell_0_in_cell_1(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex11: " << flush; insert_cell_0_in_cell_1(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex12: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex13: " << flush; insert_cell_0_in_cell_1(map,d1->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex14: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_hexahedron(map);
  d2 = make_combinatorial_hexahedron(map);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "insert vertex15: " << flush; insert_cell_0_in_cell_1(map,d1->beta(1)->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex16: " << flush; insert_cell_0_in_cell_1(map,d1->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex17: " << flush; insert_cell_0_in_cell_1(map,d1->beta(0));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert vertex18: " << flush; insert_cell_0_in_cell_1(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
    
  cout << "***************************** TEST INSERT VERTEX 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT EDGE 3D:"
       << endl;

  /*  d1 = map.create_dart();
  d2 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge1: " << flush; insert_cell_1_in_cell_2(map,d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();*/

  /*  d1 = make_edge(map,, );
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge2: " << flush; insert_cell_1_in_cell_2(map,d1, d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();*/

  d1 = make_edge(map);
  d2 = make_edge(map);
  map.template sew<3>(d1, d2); map.template sew<3>(d1->beta(2), d2->beta(2)); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  /*  cout << "insert edge3: " << flush; insert_cell_1_in_cell_2(map,d1, d1);
      map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;*/
  map.clear();

  d1 = make_combinatorial_polygon(map,4 );
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge4: " << flush; insert_cell_1_in_cell_2(map,d1, d1->beta(1)->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,4);
  d2 = make_combinatorial_polygon(map,4);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge5: " << flush; insert_cell_1_in_cell_2(map,d1, d1->beta(1)->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,4);
  d2 = make_combinatorial_polygon(map,4);
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert edge6: " << flush; insert_cell_1_in_cell_2(map,d1, d1->beta(1)->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
  map.clear();
    
  cout << "***************************** TEST INSERT EDGE 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT DANGLING EDGE 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge1: " << flush; insert_dangling_cell_1_in_cell_2(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_edge(map);
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge2: " << flush; insert_dangling_cell_1_in_cell_2(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_edge(map);
  d2 = make_edge(map);
  map.template sew<3>(d1, d2); map.template sew<3>(d1->beta(2), d2->beta(2)); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge3: " << flush; insert_dangling_cell_1_in_cell_2(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,4);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge4: " << flush; insert_dangling_cell_1_in_cell_2(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_polygon(map,4);
  d2 = make_combinatorial_polygon(map,4);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert dangling edge5: " << flush; insert_dangling_cell_1_in_cell_2(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  cout << "***************************** TEST INSERT DANGLING EDGE 3D DONE."
       << endl;

  cout << "***************************** TEST INSERT FACET 3D:"
       << endl;

  std::vector<Dart_handle> v;

  d1 = make_combinatorial_polygon(map,4);
  v.push_back(d1); v.push_back(v[0]->beta(1));
  v.push_back(v[1]->beta(1)); v.push_back(v[2]->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert facet1: " << flush; insert_cell_2_in_cell_3(map,v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();

  d1 = make_combinatorial_polygon(map,3);
  d2 = make_combinatorial_polygon(map,3);
  map.template sew<2>(d1, d2);
  v.push_back(d1); v.push_back(v[0]->beta(1)); v.push_back(v[1]->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "insert facet2: " << flush; insert_cell_2_in_cell_3(map,v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();

  d1 = make_combinatorial_hexahedron(map);
  d2 = d1->beta(2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet3: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  v.push_back(d2); v.push_back(v[0]->beta(1)->beta(2)->beta(1)); 
  v.push_back(v[1]->beta(1)->beta(2)->beta(1)); v.push_back(v[2]->beta(1)->beta(2)->beta(1));
  cout << "insert facet3: " << flush; insert_cell_2_in_cell_3(map,v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();
    
  d1 = make_combinatorial_hexahedron(map);
  d2 = make_combinatorial_hexahedron(map);
  map.template sew<3>(d1,d2);   
  d3 = d1->beta(2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "remove facet4: " << flush; CGAL::remove_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  v.push_back(d3); v.push_back(v[0]->beta(1)->beta(2)->beta(1)); 
  v.push_back(v[1]->beta(1)->beta(2)->beta(1)); v.push_back(v[2]->beta(1)->beta(2)->beta(1));
  cout << "insert facet4: " << flush; insert_cell_2_in_cell_3(map,v.begin(),v.end());
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear(); v.clear();    
    
  cout << "***************************** TEST INSERT FACET 3D DONE."
       << endl;

  /*
  cout << "***************************** TEST EDGE CONTRACTION 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge1: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge2: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map,, );
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge3: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map,, );
  map.template sew<1>(d1, d1); map.template sew<1>(d1->beta(2), d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge4: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map,, );
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge5: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_triangle(map,, , );
  d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "contract edge6: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "contract edge7: " << flush; CGAL::contract_cell<Map,1>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "contract edge8: " << flush; CGAL::contract_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<3>(d1, d2); d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "contract edge9: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout)<<std::flush << ", valid=" << map.is_valid() << endl;

  cout << "contract edge10: " << flush; CGAL::contract_cell<Map,1>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  cout << "contract edge11: " << flush; CGAL::contract_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge12: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<2>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge13: " << flush; CGAL::contract_cell<Map,1>(map,d1->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge14: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_combinatorial_hexahedron(map,, 1);
  d2 = make_combinatorial_hexahedron(map,, 1);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d3 = insert_cell_0_in_cell_2(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  std::vector<Dart_handle> V;
  {
    for ( typename Map::Dart_iterator_of_vertex_2 it =
	    map.dart_iterator_of_vertex_2_begin(d3);
	  it!=map.dart_iterator_of_vertex_2_end(d3); ++it)
      V.push_back(*it);
  }

  for (typename std::vector<Dart_handle>::iterator it=V.begin(); it!=V.end(); ++it)
    {
      cout << "contract edge15: " << flush; CGAL::contract_cell<Map,1>(map,*it);
      map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
    }
  map.clear();
   
  cout << "***************************** TEST EDGE CONTRACTION 3D DONE."
       << endl;

  cout << "***************************** TEST FACET CONTRACTION 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet1: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet2: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map,, );
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet3: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_edge(map,, );
  map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet4: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_edge(map,, );
  map.template sew<1>(d1, d1);map.template sew<1>(d1->beta(2), d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet5: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_triangle(map,, , );
  d2 = d1->beta(0);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge6: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet6: " << flush; CGAL::contract_cell<Map,2>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<3>(d1, d2); d2 = d1->beta(0); d3 = d1->beta(1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge7: " << flush; CGAL::contract_cell<Map,1>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet7: " << flush; CGAL::contract_cell<Map,2>(map,d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
   
  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<2>(d1,d2); d3=d1->beta(0);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge8: " << flush; CGAL::contract_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet8: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  d3 = make_triangle(map,, , );
  map.template sew<2>(d1,d2); map.template sew<2>(d1->beta(0),d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge9: " << flush; CGAL::contract_cell<Map,1>(map,d1->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet9: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  d3 = make_triangle(map,, , );
  d4 = make_triangle(map,, , );
  map.template sew<2>(d1,d2); map.template sew<2>(d3,d4); map.template sew<3>(d1,d3); map.template sew<3>(d2,d4); 
  d3=d1->beta(0);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge10: " << flush; CGAL::contract_cell<Map,1>(map,d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet10: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_triangle(map,, , );
  d4 = make_triangle(map,, , );
  map.template sew<3>(d1,d4);
  d2 = make_triangle(map,, , );
  d4 = make_triangle(map,, , );
  map.template sew<3>(d2,d4);   
  d3 = make_triangle(map,, , );
  d4 = make_triangle(map,, , );
  map.template sew<3>(d3,d4);
  map.template sew<2>(d1,d2);          map.template sew<2>(d1->beta(3),d2->beta(3));  
  map.template sew<2>(d1->beta(0),d3); map.template sew<2>(d1->beta(0,3),d3->beta(3));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge11: " << flush; CGAL::contract_cell<Map,1>(map,d1->beta(1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet11: " << flush; CGAL::contract_cell<Map,2>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
     
  cout << "***************************** TEST FACET CONTRACTION 3D DONE."
       << endl;

  cout << "***************************** TEST VOLUME CONTRACTION 3D:"
       << endl;

  d1 = map.create_dart();
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume1: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = map.create_dart(); map.template sew<1>(d1, d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume2: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;

  d1 = make_edge(map,, );
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume3: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_edge(map,, );
  map.template sew<1>(d1, d1);map.template sew<1>(d1->beta(2), d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume4: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_triangle(map,, , );
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume5: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<2>(d1, d2);map.template sew<2>(d1->beta(0), d2->beta(1));
  map.template sew<2>(d1->beta(1), d2->beta(0));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume6: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_triangle(map,, , );
  d2 = make_triangle(map,, , );
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume7: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();

  d1 = make_combinatorial_tetrahedron(map,, , ,
			);
  d2 = make_combinatorial_tetrahedron(map,, , ,
			);
  map.template sew<3>(d1, d2);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge8: " << flush; CGAL::contract_cell<Map,1>(map,d1->beta(1,2,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet8: " << flush; CGAL::contract_cell<Map,2>(map,d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet8: " << flush; CGAL::contract_cell<Map,2>(map,d1->beta(1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume8: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  d1 = make_combinatorial_tetrahedron(map,, , ,
			);
  d2 = make_combinatorial_tetrahedron(map,, , ,
			);
  d3 = make_combinatorial_tetrahedron(map,, , ,
			);
  map.template sew<3>(d1, d2); map.template sew<3>(d1->beta(2), d3);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract edge9: " << flush; CGAL::contract_cell<Map,1>(map,d1->beta(1,2,1));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet9: " << flush; CGAL::contract_cell<Map,2>(map,d1->beta(2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract facet9: " << flush; CGAL::contract_cell<Map,2>(map,d1->beta(1,2));
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  cout << "contract volume8: " << flush; CGAL::contract_cell<Map,3>(map,d1);
  map.display_characteristics(cout) << ", valid=" << map.is_valid() << endl;
  map.clear();
   
  cout << "***************************** TEST VOLUME CONTRACTION 3D DONE."
       << endl;
  */
}

#endif // CGAL_COMBINATORIAL_MAP_3_TEST
