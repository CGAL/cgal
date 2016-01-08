// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
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

#ifndef COMBINATORIAL_MAP_TEST_ITERATORS_H
#define COMBINATORIAL_MAP_TEST_ITERATORS_H

// This function allows to test if the iterators iterate through the dart exactly once.
template<typename CMap>
bool test_iterators_2(CMap & cm)
{
  const int NBTESTS = 9;

  typename CMap::size_type marks[NBTESTS];
  typename CMap::size_type nbdarts[NBTESTS];
  typename CMap::size_type nb;

  for (int i=0; i<NBTESTS; ++i)
  {
    marks[i]=cm.get_new_mark();
    nbdarts[i]=0;
  }

  for (typename CMap::Dart_range::iterator it=cm.darts().begin(),
       itend=cm.darts().end(); it!=itend; ++it)
  {
    for (int i=0; i<NBTESTS; ++i)
    {
      if ( !cm.is_marked(it, marks[i]) )
      {
       switch(i)
       {
         case 0: nb = CGAL::mark_cell<CMap, 0>(cm, it, marks[i]); break;
         case 1: nb = CGAL::mark_cell<CMap, 1>(cm, it, marks[i]); break;
         case 2: nb = CGAL::mark_cell<CMap, 2>(cm, it, marks[i]); break;
         case 3: nb = CGAL::mark_cell<CMap, 3>(cm, it, marks[i]); break;
         case 4: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0> >(cm, it, marks[i]); break;
         case 5: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1> >(cm, it, marks[i]); break;
         case 6: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2> >(cm, it, marks[i]); break;
         case 7: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,2> >(cm, it, marks[i]); break;
         case 8: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,2> >(cm, it, marks[i]); break;
       };
       nbdarts[i] += nb;
      }
    }
  }

  bool res = true;
  for (int i=0; i<NBTESTS; ++i)
  {
    if(nbdarts[i]!=cm.number_of_darts()) res=false;
    cm.free_mark( marks[i] );
  }
  return res;
}

// This function allows to test if the iterators iterate through the dart exactly once.
template<typename CMap>
bool test_iterators_3(CMap & cm)
{
  const int NBTESTS = 18;

  typename CMap::size_type marks[NBTESTS];
  typename CMap::size_type nbdarts[NBTESTS];
  typename CMap::size_type nb;

  for (int i=0; i<NBTESTS; ++i)
  {
    marks[i]=cm.get_new_mark();
    nbdarts[i]=0;
  }
  
  for (typename CMap::Dart_range::iterator it=cm.darts().begin(),
       itend=cm.darts().end(); it!=itend; ++it)
  {
    for (int i=0; i<NBTESTS; ++i)
    {
      if ( !cm.is_marked(it, marks[i]) )
      {
       switch(i)
       {
         case 0: nb = CGAL::mark_cell<CMap, 0>(cm, it, marks[i]); break;
         case 1: nb = CGAL::mark_cell<CMap, 1>(cm, it, marks[i]); break;
         case 2: nb = CGAL::mark_cell<CMap, 2>(cm, it, marks[i]); break;
         case 3: nb = CGAL::mark_cell<CMap, 3>(cm, it, marks[i]); break;
         case 4: nb = CGAL::mark_cell<CMap, 4>(cm, it, marks[i]); break;
         case 5: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0> >(cm, it, marks[i]); break;
         case 6: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1> >(cm, it, marks[i]); break;
         case 7: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2> >(cm, it, marks[i]); break;
         case 8: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 3> >(cm, it, marks[i]); break;
         case 9: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,2> >(cm, it, marks[i]); break;
         case 10: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,2> >(cm, it, marks[i]); break;
         case 11: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,3> >(cm, it, marks[i]); break;
         case 12: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,3> >(cm, it, marks[i]); break;
         case 13: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2,3> >(cm, it, marks[i]); break;
         case 14: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,2,3> >(cm, it, marks[i]); break;
         case 15: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,2,3> >(cm, it, marks[i]); break;
         case 16: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_involution<CMap,3> >(cm, it, marks[i]); break;
         case 17: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_involution_inv<CMap,3> >(cm, it, marks[i]); break;
       };
       nbdarts[i] += nb;
      }
    }
  }

  bool res = true;  
  for (int i=0; i<NBTESTS; ++i)
  {
    if(nbdarts[i]!=cm.number_of_darts()) res=false;
    cm.free_mark( marks[i] );
  }
  return res;
}

// This function allows to test if the iterators iterate through the dart exactly once.
template<typename CMap>
bool test_iterators_4(CMap & cm)
{
  const int NBTESTS = 31;

  typename CMap::size_type marks[NBTESTS];
  typename CMap::size_type nbdarts[NBTESTS];
  typename CMap::size_type nb;

  for (int i=0; i<NBTESTS; ++i)
  {
    marks[i]=cm.get_new_mark();
    nbdarts[i]=0;
  }

  for (typename CMap::Dart_range::iterator it=cm.darts().begin(),
       itend=cm.darts().end(); it!=itend; ++it)
  {
    for (int i=0; i<NBTESTS; ++i)
    {
      if ( !cm.is_marked(it, marks[i]) )
      {
       switch(i)
       {
         case 0: nb = CGAL::mark_cell<CMap, 0>(cm, it, marks[i]); break;
         case 1: nb = CGAL::mark_cell<CMap, 1>(cm, it, marks[i]); break;
         case 2: nb = CGAL::mark_cell<CMap, 2>(cm, it, marks[i]); break;
         case 3: nb = CGAL::mark_cell<CMap, 3>(cm, it, marks[i]); break;
         case 4: nb = CGAL::mark_cell<CMap, 4>(cm, it, marks[i]); break;
         case 5: nb = CGAL::mark_cell<CMap, 5>(cm, it, marks[i]); break;
         case 6: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0> >(cm, it, marks[i]); break;
         case 7: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1> >(cm, it, marks[i]); break;
         case 8: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2> >(cm, it, marks[i]); break;
         case 9: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 3> >(cm, it, marks[i]); break;
         case 10: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 4> >(cm, it, marks[i]); break;
         case 11: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,2> >(cm, it, marks[i]); break;
         case 12: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,2> >(cm, it, marks[i]); break;
         case 13: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,3> >(cm, it, marks[i]); break;
         case 14: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,3> >(cm, it, marks[i]); break;
         case 15: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2,3> >(cm, it, marks[i]); break;
         case 16: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,4> >(cm, it, marks[i]); break;
         case 17: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,4> >(cm, it, marks[i]); break;
         case 18: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2,4> >(cm, it, marks[i]); break;
         case 19: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 3,4> >(cm, it, marks[i]); break;
         case 20: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,2,3> >(cm, it, marks[i]); break;
         case 21: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,2,3> >(cm, it, marks[i]); break;
         case 22: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,2,4> >(cm, it, marks[i]); break;
         case 23: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,2,4> >(cm, it, marks[i]); break;
         case 24: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 0,3,4> >(cm, it, marks[i]); break;
         case 25: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 1,3,4> >(cm, it, marks[i]); break;
         case 26: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap, 2,3,4> >(cm, it, marks[i]); break;
         case 27: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_involution<CMap,3> >(cm, it, marks[i]); break;
         case 28: nb = CGAL::mark_orbit<CMap,
               CGAL::CMap_dart_const_iterator_basic_of_involution_inv<CMap,3> >(cm, it, marks[i]); break;
         case 29: nb = CGAL::mark_orbit<CMap,
                CGAL::CMap_dart_const_iterator_basic_of_involution<CMap,4> >(cm, it, marks[i]); break;
         case 30: nb = CGAL::mark_orbit<CMap,
                CGAL::CMap_dart_const_iterator_basic_of_involution_inv<CMap,4> >(cm, it, marks[i]); break;
       };
       nbdarts[i] += nb;
      }
    }
  }

  bool res = true;
  for (int i=0; i<NBTESTS; ++i)
  {
    if(nbdarts[i]!=cm.number_of_darts()) res=false;
    cm.free_mark( marks[i] );
  }
  return res;
}

#endif // COMBINATORIAL_MAP_TEST_ITERATORS_H
