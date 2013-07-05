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
#ifndef CGAL_COMBINATORIAL_MAP_SEWABLE_H
#define CGAL_COMBINATORIAL_MAP_SEWABLE_H

#include <CGAL/Dart_const_iterators.h>
#include <CGAL/Unique_hash_map.h>

/* Definition of functor used to test if two darts are i-sewable
 * (we use functors as there are different specializations).
 */
namespace CGAL
{
#define CGAL_BETAINV(i) (i>1?i:(i==1?0:1))
namespace internal
{
// Generic case for 1<=i<=dimension, and 3<dim.
template<typename CMap, unsigned int i, unsigned int dim=CMap::dimension>
struct Is_sewable_functor
{
  static bool run( const CMap* amap,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion( 1<=i && i<=CMap::dimension );
    CGAL_assertion( 3<dim );
    CGAL_assertion( adart1!=NULL && adart2!=NULL );
    CGAL_assertion( adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle );

    if ( !adart1->template is_free<i>() ||
         !adart2->template is_free<CGAL_BETAINV(i)>() )
      return false;

    if ( adart1==adart2 )
    {
      if ( i==1 ) return true;
      return false;
    }

    // hash map to build the isomorphism between the two i-cells.
    CGAL::Unique_hash_map<typename CMap::Dart_const_handle,
        typename CMap::Dart_const_handle> bijection;

    int m1 = amap->get_new_mark();
    int m2 = amap->get_new_mark();
    CGAL::CMap_dart_const_iterator_basic_of_involution<CMap,i>
        I1(*amap, adart1, m1);
    CGAL::CMap_dart_const_iterator_basic_of_involution_inv<CMap,i>
        I2(*amap, adart2, m2);
    bool res = true;
    int mbijection = amap->get_new_mark();

    while ( res && I1.cont() && I2.cont() )
    {
      amap->mark(I1, mbijection);
      bijection[I1]=I2;

      CGAL_assertion( I1->template is_free<i>() );
      CGAL_assertion( I2->template is_free<CGAL_BETAINV(i)>() );

      // We can remove this constraint which is not required for
      // combinatorial map definition, but which is quite "normal"
      // Indeed in this case we try to i-sew an i-cell with itself (case
      // of folded cells).
      if ( i>1 && (I1==adart2 || I2==adart1) ) res=false;

      if ( i>2)
      {
        if ( I1->template is_free<1>() )
        {
          if ( !I2->template is_free<0>() ) res=false;
        }
        else
        {
          if ( I2->template is_free<0>() ) res=false;
          else if ( amap->is_marked(I1->template beta<1>(), mbijection) )
          {
            if ( bijection[I1->template beta<1>()]!=I2->template beta<0>() )
              res=false;
          }
        }
      }

      for ( unsigned int j=2; res && j<=CMap::dimension; ++j )
      {
        if ( j+1!=i && j!=i && j!=i+1 )
        {
          if ( I1->is_free(j) )
          {
            if ( !I2->is_free(j) ) res=false;
          }
          else
          {
            if ( I2->is_free(j) ) res=false;
            else if ( amap->is_marked(I1->beta(j), mbijection) )
            {
              if ( bijection[I1->beta(j)]!=I2->beta(j) ) res=false;
            }
          }
        }
      }
      ++I1; ++I2;
    }
    if ( I1.cont()!=I2.cont() )
      res = false;

    amap->negate_mark(m1);
    amap->negate_mark(m2);
    I1.rewind(); I2.rewind();
    while ( amap->number_of_marked_darts(mbijection)>0 )
    {
      amap->unmark(I1, mbijection);
      ++I1; ++I2;
    }

    CGAL_assertion( amap->is_whole_map_marked(m1) );
    CGAL_assertion( amap->is_whole_map_marked(m2) );
    CGAL_assertion( amap->is_whole_map_unmarked(mbijection) );
    amap->free_mark(m1);
    amap->free_mark(m2);
    amap->free_mark(mbijection);

    return res;
  }
};

// Specialization for i=0 and 3<dim.
template<typename CMap, unsigned int dim>
struct Is_sewable_functor<CMap, 0, dim>
{
  static bool run( const CMap* amap,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  { return Is_sewable_functor<CMap,1, dim>::run(amap, adart2, adart1); }
};

// Specialization for i=0 and dim=1.
template<typename CMap>
struct Is_sewable_functor<CMap, 0, 1>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<0>() ||
         !adart2->template is_free<1>() )
      return false;
    return true;
  }
};
// Specialization for i=1 and dim=1.
template<typename CMap>
struct Is_sewable_functor<CMap, 1, 1>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<1>() ||
         !adart2->template is_free<0>() )
      return false;
    return true;
  }
};

// Specialization for i=0 and dim=2.
template<typename CMap>
struct Is_sewable_functor<CMap, 0, 2>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<0>() ||
         !adart2->template is_free<1>() )
      return false;
    return true;
  }
};
// Specialization for i=1 and dim=2.
template<typename CMap>
struct Is_sewable_functor<CMap, 1, 2>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<1>() ||
         !adart2->template is_free<0>() )
      return false;
    return true;
  }
};
// Specialization for i=2 and dim=2.
template<typename CMap>
struct Is_sewable_functor<CMap, 2, 2>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<2>() ||
         !adart2->template is_free<2>() || adart1==adart2 )
      return false;
    return true;
  }
};

// Specialization for i=0 and dim=3.
template<typename CMap>
struct Is_sewable_functor<CMap, 0, 3>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<0>() ||
         !adart2->template is_free<1>() )
      return false;

    if ( adart1->template is_free<3>() )
    {
      if ( !adart2->template is_free<3>() ) return false;
      return true;
    }

    // Here adart1 is not 3-free
    if ( adart2->template is_free<3>() ) return false;

    CGAL_assertion( adart1->template beta<3>()->template is_free<1>() &&
                    adart2->template beta<3>()->template is_free<0>() );
    return true;
  }
};
// Specialization for i=1 and dim=3.
template<typename CMap>
struct Is_sewable_functor<CMap, 1, 3>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<1>() ||
         !adart2->template is_free<0>() )
      return false;

    if ( adart1->template is_free<3>() )
    {
      if ( !adart2->template is_free<3>() ) return false;
      return true;
    }

    // Here adart1 is not 3-free
    if ( adart2->template is_free<3>() ) return false;

    CGAL_assertion( adart1->template beta<3>()->template is_free<0>() &&
                    adart2->template beta<3>()->template is_free<1>() );
    return true;
  }
};
// Specialization for i=2 and dim=3.
template<typename CMap>
struct Is_sewable_functor<CMap, 2, 3>
{
  static bool run( const CMap* /*amap*/,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<2>() ||
         !adart2->template is_free<2>() || adart1==adart2 )
      return false;
    return true;
  }
};
// Specialization for i=3 and dim=3.
template<typename CMap>
struct Is_sewable_functor<CMap, 3, 3>
{
  static bool run( const CMap* amap,
                   typename CMap::Dart_const_handle adart1,
                   typename CMap::Dart_const_handle adart2 )
  {
    CGAL_assertion(adart1!=NULL && adart2!=NULL);
    CGAL_assertion(adart1!=CMap::null_dart_handle &&
        adart2!=CMap::null_dart_handle);

    if ( !adart1->template is_free<3>() ||
         !adart2->template is_free<3>() )
      return false;

    CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap,1> I1(*amap, adart1);
    CGAL::CMap_dart_const_iterator_basic_of_orbit<CMap,0> I2(*amap, adart2);
    bool res=true;
    while (res && I1.cont() && I2.cont())
    {
      CGAL_assertion( I1->template is_free<3>() ||
                      I2->template is_free<3>() );

      // We can remove this constraint which is not required for
      // combinatorial map definition, but which is quite "normal" as it avoid
      // fold cells.
      if ( I1==adart2 || I2==adart1 ) res=false;
      else if ( I1.prev_operation()!=I2.prev_operation() ) res=false;

      ++I1; ++I2;
    }
    if ( I1.cont()!=I2.cont() )
      res=false;

    return res;
  }
};

} //namespace internal

} //namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_SEWABLE_H
//******************************************************************************
