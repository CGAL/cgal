// Copyright (c) 2016 CNRS and LIRIS' Establishments (France).
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
#ifndef CGAL_GENERALIZED_MAP_SEWABLE_H
#define CGAL_GENERALIZED_MAP_SEWABLE_H

#include <CGAL/GMap_dart_const_iterators.h>
#include <CGAL/Unique_hash_map.h>

/* Definition of functor used to test if two darts are i-sewable
 * (we use functors as there are different specializations).
 * @todo Specializations ?
 */
namespace CGAL
{
namespace internal
{
// Generic case (and for the moment the only one).
template<typename GMap, unsigned int i, unsigned int dim=GMap::dimension>
struct GMap_is_sewable_functor
{
  static bool run( const GMap& amap,
                   typename GMap::Dart_const_handle adart1,
                   typename GMap::Dart_const_handle adart2 )
  {
    CGAL_assertion( i<=GMap::dimension );
    if ( !amap.template is_free<i>(adart1) ||
         !amap.template is_free<i>(adart2) )
      return false;

    if ( adart1==adart2 )
    {
      if ( i==1 ) return true;
      return false;
    }

    // hash map to build the isomorphism between the two i-cells.
    CGAL::Unique_hash_map<typename GMap::Dart_const_handle,
        typename GMap::Dart_const_handle,
        typename GMap::Hash_function> bijection;

    typename GMap::size_type m1 = amap.get_new_mark();
    typename GMap::size_type m2 = amap.get_new_mark();
    CGAL::GMap_dart_const_iterator_basic_of_involution<GMap,i>
        I1(amap, adart1, m1);
    CGAL::GMap_dart_const_iterator_basic_of_involution<GMap,i>
        I2(amap, adart2, m2);
    bool res = true;
    typename GMap::size_type mbijection = amap.get_new_mark();

    while ( res && I1.cont() && I2.cont() )
    {
      amap.mark(I1, mbijection);
      bijection[I1]=I2;

      CGAL_assertion( amap.template is_free<i>(I1) );
      CGAL_assertion( amap.template is_free<i>(I2) );

      // We can remove this constraint which is not required for
      // generalized map definition, but which is quite "normal"
      // Indeed in this case we try to i-sew an i-cell with itself (case
      // of folded cells).
      if ( I1==adart2 || I2==adart1 ) res=false;

      for ( unsigned int j=0; res && j<=GMap::dimension; ++j )
      {
        if ( j+1!=i && j!=i && j!=i+1 )
        {
          if ( amap.is_free(I1,j) )
          {
            if ( !amap.is_free(I2,j) ) res=false;
          }
          else
          {
            if ( amap.is_free(I2,j) ) res=false;
            else if ( amap.is_marked(amap.alpha(I1,j), mbijection) )
            {
              if ( bijection[amap.alpha(I1,j)]!=amap.alpha(I2,j) ) res=false;
            }
          }
        }
      }
      ++I1; ++I2;
    }
    if ( I1.cont()!=I2.cont() )
      res = false;

    amap.negate_mark(m1);
    amap.negate_mark(m2);
    I1.rewind(); I2.rewind();
    while ( amap.number_of_marked_darts(mbijection)>0 )
    {
      amap.unmark(I1, mbijection);
      ++I1; ++I2;
    }

    CGAL_assertion( amap.is_whole_map_marked(m1) );
    CGAL_assertion( amap.is_whole_map_marked(m2) );
    CGAL_assertion( amap.is_whole_map_unmarked(mbijection) );
    amap.free_mark(m1);
    amap.free_mark(m2);
    amap.free_mark(mbijection);

    return res;
  }
};

} //namespace internal

} //namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_SEWABLE_H
//******************************************************************************
