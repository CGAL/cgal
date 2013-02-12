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
#ifndef CGAL_COMBINATORIAL_MAP_FUNCTORS_H
#define CGAL_COMBINATORIAL_MAP_FUNCTORS_H

#include <CGAL/Dart_const_iterators.h>
#include <CGAL/Cell_const_iterators.h>
#include <CGAL/Combinatorial_map_basic_operations.h>
#include <CGAL/Unique_hash_map.h>
#include <stack>

/* Definition of functors used to manage attributes (we need functors as
 * attributes are stored in tuple, thus all the access must be done at
 * compiling time). Some of these functors are used with
 * Foreach_enabled_attributes to iterate through all the non void attribs.
 *
 * Call_split_functor<CMap,i> to call the OnSplit functors on two given
 *    i-attributes.
 *
 * Call_merge_functor<CMap,i> to call the OnMerge functors on two given
 *    i-attributes.
 *
 * Reserve_mark_functor<CMap> to reserve one mark for each non void attribute.
 *
 * Test_is_valid_attribute_functor<CMap> to test if a given i-cell is valid
 *   (all its darts are linked to the same attribute, no other dart is linked
 *    with this attribute).
 *
 * Count_cell_functor<CMap> to count the nuber of i-cells.
 *
 * Count_bytes_one_attribute_functor<CMap> to count the memory occupied by
 *  i-attributes.
 *
 * Decrease_attribute_functor<CMap> to decrease by one the ref counting of
 *  a given i-attribute.
 *
 * Beta_functor<Dart, i...> to call several beta on the given dart. Indices are
 *   given as parameter of the run function.
 *
 * Beta_functor_static<Dart, i...> to call several beta on the given dart.
 *   Indices are given as template arguments.
 *
 * Group_attribute_functor_of_dart<CMap> to group the <i>-attributes of two
 *    given darts (except for j-dim). Only the attributes of the two given
 *    darts are possibly modified.
 *
 * Group_attribute_functor_of_dart_run<CMap,i> same than
 *   Group_attribute_functor_of_dart<CMap>::run<i>, with i template argument
 *   given in the struct to enable specialization.
 *
 * Group_attribute_functor<CMap> to group the <i>-attributes of two
 *    given i-cells (except for j-adim). If one i-attribute is NULL, we set the
 *    darts of its i-cell to the second attribute. If both i-attributes are
 *    non NULL, we overide all the i-attribute of the second i-cell to the
 *    first i-attribute.
 *
 * Degroup_attribute_functor_run<CMap> to degroup one i-attributes in two
 *   (except for j-adim).
 *
 * Test_split_attribute_functor<CMap,i> to test if there is some i-attributes
 *   that are split after an operation. Modified darts are given in a
 *   std::deque.
 */

namespace CGAL
{
  namespace internal
  {

  /** @file Combinatorial_map_functors.h
   * Definition of functors used for dD Combinatorial map.
   */
  // **************************************************************************
  // Functor which call Functor::operator() on the two given cell_attributes
  template<typename Cell_attribute, typename Functor>
  struct Apply_cell_functor
  {
    static void run(Cell_attribute& acell1, Cell_attribute& acell2)
    {
      Functor() (acell1,acell2);
    }
  };
  //...except for Null_functor.
  template<typename Cell_attribute>
  struct Apply_cell_functor<Cell_attribute,Null_functor>
  {
    static void run(Cell_attribute&, Cell_attribute&)
    {}
  };
  // **************************************************************************
  // Functor used to call the On_split functor between the two given darts.
  template<typename CMap, unsigned int i,
           typename Enabled=typename CMap::Helper::
         #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
           template
         #endif
           Attribute_type<i>::type>
  struct Call_split_functor
  {
    static void run(typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    {
      Apply_cell_functor
          <typename CMap::Helper::template Attribute_type<i>::type,
          typename CMap::Helper::
          template Attribute_type<i>::type::On_split>::
          run(*(adart1->template attribute<i>()),
              *(adart2->template attribute<i>()));
    }
    static void
    run(typename CMap::Helper::template Attribute_handle<i>::type a1,
        typename CMap::Helper::template Attribute_handle<i>::type a2)
    {
      Apply_cell_functor
          <typename CMap::Helper::template Attribute_type<i>::type,
          typename CMap::Helper::
          template Attribute_type<i>::type::On_split>::
          run(*a1, *a2);
    }
  };
  // Specialization for disabled attributes.
  template<typename CMap,unsigned int i>
  struct Call_split_functor<CMap,i,CGAL::Void>
  {
    static void run(typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // **************************************************************************
  // Functor used to call the On_merge functor between the two given darts.
  template<typename CMap,unsigned int i,
           typename Enabled=typename CMap::Helper::
         #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
           template
         #endif
           Attribute_type<i>::type>
  struct Call_merge_functor
  {
    static void run(typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    {
      Apply_cell_functor
          <typename CMap::Helper::template Attribute_type<i>::type,
          typename CMap::Helper::template Attribute_type<i>::type::On_merge>::
          run(*(adart1->template attribute<i>()),
              *(adart2->template attribute<i>()));
    }
    static void
    run(typename CMap::Helper::template Attribute_handle<i>::type a1,
        typename CMap::Helper::template Attribute_handle<i>::type a2)
    {
      Apply_cell_functor
          <typename CMap::Helper::template Attribute_type<i>::type,
          typename CMap::Helper::template Attribute_type<i>::type::On_merge>::
          run(*a1, *a2);
    }
  };
  // Specialization for disabled attributes.
  template<typename CMap,unsigned int i>
  struct Call_merge_functor<CMap,i,CGAL::Void>
  {
    static void run(typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // **************************************************************************
  /// Functor used to reserve one mark for each enabled attribute.
  template<typename CMap>
  struct Reserve_mark_functor
  {
    template <unsigned int i>
    static void run(const CMap* amap, std::vector<int>* marks)
    { (*marks)[i] = amap->get_new_mark(); }
  };
  // **************************************************************************
  /// Functor used to test if a cell is valid
  template<typename CMap>
  struct Test_is_valid_attribute_functor
  {
    template <unsigned int i>
    static void run(const CMap* amap,
                    typename CMap::Dart_const_handle adart,
                    std::vector<int>* marks, bool *ares)
    {
      if (!amap->template is_valid_attribute<i>(adart,(*marks)[i]) )
      {
        (*ares)=false;
        std::cerr << "CMap not valid: a "<<i<<"-cell is not correctly "
          "associated with an attribute for " << &(*adart)<< std::endl;
      }
    }
  };
  // **************************************************************************
  /// Functor for counting i-cell
  template<typename CMap>
  struct Count_cell_functor
  {
    template <unsigned int i>
    static void run( const CMap* amap,
                     typename CMap::Dart_const_handle adart,
                     std::vector<int>* amarks,
                     std::vector<unsigned int>* ares )
    {
      if ( (*amarks)[i]!=-1 && !amap->is_marked(adart, (*amarks)[i]) )
      {
        ++ (*ares)[i];
        mark_cell<CMap,i>(*amap, adart, (*amarks)[i]);
      }
    }
  };
  // **************************************************************************
  /// Functor for counting the memory occupation of attributes
  /// Be careful not reentrant !!! TODO a  Foreach_enabled_attributes
  /// taking an instance of a functor as argument allowing to compute
  /// and return values.
  template<typename CMap>
  struct Count_bytes_one_attribute_functor
  {
    template <unsigned int i>
    static void run( const CMap* amap )
    {
      res += amap->template attributes<i>().capacity()*
        sizeof(typename CMap::template Attribute_type<i>::type);
    }

    static typename CMap::size_type res;
  };
  template<typename CMap>
  typename CMap::size_type Count_bytes_one_attribute_functor<CMap>::res = 0;

  template<typename CMap>
  struct Count_bytes_all_attributes_functor
  {
    static typename CMap::size_type run( const CMap& amap )
    {
      Count_bytes_one_attribute_functor<CMap>::res = 0;
      CMap::Helper::template Foreach_enabled_attributes
        <Count_bytes_one_attribute_functor<CMap> >::run(&amap);
      return Count_bytes_one_attribute_functor<CMap>::res;
    }
  };
  // **************************************************************************
  /// Functor used to call decrease_attribute_ref_counting<i>
  /// on each i-cell attribute enabled
  template<typename CMap>
  struct Decrease_attribute_functor
  {
    template <unsigned int i>
    static void run(CMap* amap, typename CMap::Dart_handle adart)
    { amap->template decrease_attribute_ref_counting<i>(adart); }
  };
  // **************************************************************************
  /// Functor used for link_beta to update the i-attributes of
  /// adart2 on the attributes of this dart, except if i=j.
  ///    (j is the dimension of the beta modified between adart1 and adart2,
  ///     so that after the modification we will have beta_j(adart1)==adart2)
  /// Only attributes of dh1 or dh2 can be modified.
  template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
           typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Group_attribute_functor_of_dart_run
  {
    /// Group the i-attribute of dh1 and dh2.
    static void run(CMap* amap,
                    typename CMap::Dart_handle dh1,
                    typename CMap::Dart_handle dh2)
    {
      CGAL_static_assertion( 1<=i && i<=CMap::dimension );
      CGAL_static_assertion( i!=j );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<i>::value>=0,
                                "Group_attribute_functor_of_dart_run<i> but "
                                "i-attributes are disabled");
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh1!=CMap::null_dart_handle );

      typename CMap::Helper::template Attribute_handle<i>::type
          a1=dh1->template attribute<i>();
      typename CMap::Helper::template Attribute_handle<i>::type
          a2=dh2->template attribute<i>();

      // If the two attributes are equal, nothing to do.
      if ( a1==a2 ) return;

      if ( a1==NULL ) amap->template set_attribute_of_dart<i>(dh1, a2);
      else            amap->template set_attribute_of_dart<i>(dh2, a1);
    }
  };
  // Specialization for i=0 and 2<=j.
  template<typename CMap, unsigned int j, typename T>
  struct Group_attribute_functor_of_dart_run<CMap, 0, j, T>
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle dh1,
                    typename CMap::Dart_handle dh2)
    {
      CGAL_static_assertion(j!=0 && j!=1);
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh1!=CMap::null_dart_handle );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<0>::value>=0,
                                "Group_attribute_functor_of_dart_run<0> but "
                                "0-attributes are disabled");
      typename CMap::Helper::template Attribute_handle<0>::type
          a1=NULL, a2=NULL;

      // First extremity
      typename CMap::Dart_handle od = dh2->other_extremity();
      if ( od!=NULL )
      {
        a1=dh1->template attribute<0>();
        a2=od->template attribute<0>();

        if ( a1==NULL && a2!=NULL )
        {
          amap->template set_attribute_of_dart<0>(dh1, a2);
        }
      }

      // Second extremity
      od = dh1->other_extremity();
      if ( od!=NULL )
      {
        a1=od->template attribute<0>();
        a2=dh2->template attribute<0>();

        if ( a1!=NULL )
        {
          amap->template set_attribute_of_dart<0>(dh2, a1);
        }
      }
    }
  };
  // Specialization for i=0 and j=0.
  template<typename CMap, typename T>
  struct Group_attribute_functor_of_dart_run<CMap, 0, 0, T>
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle dh1,
                    typename CMap::Dart_handle dh2)
    {
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh1!=CMap::null_dart_handle );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<0>::value>=0,
                                "Group_attribute_functor_of_dart_run<0> but "
                                "0-attributes are disabled");
      typename CMap::Dart_handle od = dh2->other_extremity();
      if ( od!=NULL )
      {
        typename CMap::Helper::template Attribute_handle<0>::type
            a1=dh1->template attribute<0>();
        typename CMap::Helper::template Attribute_handle<0>::type
            a2=od->template attribute<0>();

        if ( a1==NULL && a2!=NULL )
        {
          amap->template set_attribute_of_dart<0>(dh1, a2);
        }
      }
    }
  };
  // Specialization for i=0 and j=1.
  template<typename CMap, typename T>
  struct Group_attribute_functor_of_dart_run<CMap, 0, 1, T>
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle dh1,
                    typename CMap::Dart_handle dh2)
    {
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh1!=CMap::null_dart_handle );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<0>::value>=0,
                                "Group_attribute_functor_of_dart_run<0> but "
                                "0-attributes are disabled");
      typename CMap::Dart_handle od = dh1->other_extremity();
      if ( od!=NULL )
      {
        typename CMap::Helper::template Attribute_handle<0>::type
            a1=od->template attribute<0>();
        typename CMap::Helper::template Attribute_handle<0>::type
            a2=dh2->template attribute<0>();

        if ( a1!=NULL )
        {
          amap->template set_attribute_of_dart<0>(dh2, a1);
        }
      }
    }
  };
  // Specialization for void attributes.
  template<typename CMap, unsigned int i, unsigned int j>
  struct Group_attribute_functor_of_dart_run<CMap,i,j,CGAL::Void>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // Specialization for i=j. Do nothing as j is the dimension to not consider.
  template<typename CMap, unsigned int i, typename T>
  struct Group_attribute_functor_of_dart_run<CMap,i,i,T>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // **************************************************************************
  /// Functor used for link_beta to update the attributes of
  /// adart2 on the attributes of this dart, except for j-attributes.
  ///    (j is the dimension of the beta modified between adart1 and adart2,
  ///     so that after the modification we will have beta_j(adart1)==adart2)
  /// We define run<i> to allows to use this functor with
  /// Foreach_enabled_attributes.
  ///   If you know i at compiling time, use directly
  ///   Group_attribute_functor_of_dart_run.
  template<typename CMap, unsigned int j=CMap::dimension+1>
  struct Group_attribute_functor_of_dart
  {
    template <unsigned int i>
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    {
      Group_attribute_functor_of_dart_run<CMap,i,j>::run(amap,adart1,adart2);
    }
  };
  // **************************************************************************
  // Functor used to group the two i-attributes of the two i-cells, except the
  // attribute of j
  //    (j is the dimension of the beta modified between adart1 and adart2).
  template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
           typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Group_attribute_functor_run
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    {
      CGAL_static_assertion( 1<=i && i<=CMap::dimension );
      CGAL_static_assertion( i!=j );
      CGAL_static_assertion_msg
          ( CMap::Helper::template Dimension_index<i>::value>=0,
           "Group_attribute_functor_run<i> but i-attributes are disabled" );
      CGAL_assertion( adart1!=NULL && adart2!=NULL );
      CGAL_assertion( adart1!=CMap::null_dart_handle &&
          adart2!=CMap::null_dart_handle );

      typename CMap::Helper::template Attribute_handle<i>::type
          a1=adart1->template attribute<i>();
      typename CMap::Helper::template Attribute_handle<i>::type
          a2=adart2->template attribute<i>();

      // If the two attributes are equal, nothing to do.
      if ( a1 == a2 ) return;

      typename CMap::Dart_handle toSet = NULL;

      // If the attribute associated to adart1 is NULL, set it with
      // the attribute associated to adart2 (necessarily != NULL)
      if (a1 == NULL)
      { toSet  = adart1; a1 = a2; }
      else
      {
        toSet = adart2;
        if (a2 != NULL)
        {
          Call_merge_functor<CMap, i>::run(a1, a2);
        }
      }
      amap->template set_attribute<i>(toSet, a1);
    }
  };
  // Specialization for i=0 and 2<=j.
  template<typename CMap, unsigned int j, typename T>
  struct Group_attribute_functor_run<CMap, 0, j, T>
  {
    static void run( CMap* amap,
                     typename CMap::Dart_handle dh1,
                     typename CMap::Dart_handle dh2 )
    {
      CGAL_static_assertion_msg
          ( CMap::Helper::template Dimension_index<0>::value>=0,
           "Group_attribute_functor_run<0> but 0-attributes are disabled" );
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh2!=CMap::null_dart_handle );
      CGAL_static_assertion(j!=0 && j!=1);

      typename CMap::Helper::template Attribute_handle<0>::type
          a1=NULL, a2=NULL;
      typename CMap::Dart_handle toSet=NULL;
      // First extremity
      typename CMap::Dart_handle od=dh2->other_extremity();
      if ( od!=NULL )
      {
        a1=dh1->template attribute<0>();
        a2=od->template attribute<0>();
        if ( a1!=a2 )
        {
          if ( a1==NULL )
          { toSet=dh1; a1=a2; }
          else
          {
            toSet=od;
            if ( a2!=NULL )
            {
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
      // Second extremity
      od = dh1->other_extremity();
      if ( od!=NULL )
      {
        a1=od->template attribute<0>();
        a2=dh2->template attribute<0>();
        if ( a1!=a2 )
        {
          if ( a1==NULL )
          { toSet=od; a1=a2; }
          else
          {
            toSet=dh2;
            if ( a2!=NULL )
            {
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
    }
  };
  // Specialization for i=0 and j=0.
  template<typename CMap, typename T>
  struct Group_attribute_functor_run<CMap, 0, 0, T>
  {
    static void run( CMap* amap,
                     typename CMap::Dart_handle dh1,
                     typename CMap::Dart_handle dh2 )
    {
      CGAL_static_assertion_msg
          ( CMap::Helper::template Dimension_index<0>::value>=0,
           "Group_attribute_functor_run<0> but 0-attributes are disabled" );
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh2!=CMap::null_dart_handle );

      typename CMap::Dart_handle od=dh2->other_extremity();
      if ( od!=NULL )
      {
        typename CMap::Dart_handle toSet=NULL;
        typename CMap::Helper::template Attribute_handle<0>::type
            a1=dh1->template attribute<0>();
        typename CMap::Helper::template Attribute_handle<0>::type
            a2=od->template attribute<0>();
        if ( a1!=a2 )
        {
          if ( a1==NULL )
          { toSet=dh1; a1=a2; }
          else
          {
            toSet=od;
            if ( a2!=NULL )
            {
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
    }
  };
  // Specialization for i=0 and j=1.
  template<typename CMap, typename T>
  struct Group_attribute_functor_run<CMap, 0, 1, T>
  {
    static void run( CMap* amap,
                     typename CMap::Dart_handle dh1,
                     typename CMap::Dart_handle dh2 )
    {
      CGAL_static_assertion_msg
          ( CMap::Helper::template Dimension_index<0>::value>=0,
           "Group_attribute_functor_run<0> but 0-attributes are disabled" );
      CGAL_assertion( dh1!=NULL && dh2!=NULL );
      CGAL_assertion( dh1!=CMap::null_dart_handle &&
          dh2!=CMap::null_dart_handle );

      typename CMap::Dart_handle od = dh1->other_extremity();
      if ( od!=NULL )
      {
        typename CMap::Dart_handle toSet=NULL;
        typename CMap::Helper::template Attribute_handle<0>::type
            a1=od->template attribute<0>();
        typename CMap::Helper::template Attribute_handle<0>::type
            a2=dh2->template attribute<0>();
        if ( a1!=a2 )
        {
          if ( a1==NULL )
          { toSet=od; a1=a2; }
          else
          {
            toSet=dh2;
            if ( a2!=NULL )
            {
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
    }
  };
  // Specialization for void attributes.
  template<typename CMap, unsigned int i, unsigned int j>
  struct Group_attribute_functor_run<CMap, i, j, CGAL::Void>
  {
    static void run( CMap*,
                     typename CMap::Dart_handle,
                     typename CMap::Dart_handle )
    {}
  };
  // Specialization for i=j. Do nothing as j is the dimension to not consider.
  template<typename CMap, unsigned int i, typename T>
  struct Group_attribute_functor_run<CMap,i,i,T>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // **************************************************************************
  /// Functor used for sew to update the attributes of
  /// adart2 on the attributes of this dart, except for j-attributes.
  ///    (j is the dimension of the beta modified between adart1 and adart2,
  ///     so that after the modification we will have beta_j(adart1)==adart2)
  /// We define run<i> to allows to use this functor with
  /// Foreach_enabled_attributes.
  ///   If you know i at compiling time, use directly
  ///   Group_attribute_functor_run.
  template<typename CMap, unsigned int j=CMap::dimension+1>
  struct Group_attribute_functor
  {
    template <unsigned int i>
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    { Group_attribute_functor_run<CMap,i,j>::run(amap,adart1,adart2); }
  };

  // **************************************************************************
  // Beta functor, used to combine several beta.
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template<typename Dart_handle, typename ... Betas>
  struct Beta_functor;
  template<typename Dart_handle, typename ... Betas>
  struct Beta_functor<Dart_handle, int, Betas...>
  {
    static Dart_handle run(Dart_handle ADart, int B, Betas... betas)
    { return Beta_functor<Dart_handle, Betas...>::run(ADart->beta(B),
                                                      betas...); }
  };
  template<typename Dart_handle>
  struct Beta_functor<Dart_handle, int>
  {
    static Dart_handle run(Dart_handle ADart, int B)
    {
      CGAL_assertion( ADart!=NULL );
      return ADart->beta(B);
    }
  };
  // **************************************************************************
  template<typename Dart_handle, typename ... Betas>
  struct Beta_functor_static;
  template<typename Dart_handle, int B, typename ... Betas>
  struct Beta_functor_static<Dart_handle, B, Betas...>
  {
    static Dart_handle run(Dart_handle ADart)
    { return Beta_functor_static<Dart_handle, Betas...>::
          run(ADart->beta<B>()); }
  };
  template<typename Dart_handle, int B>
  struct Beta_functor_static<Dart_handle, B>
  {
    static Dart_handle run(Dart_handle ADart)
    {
      CGAL_assertion( ADart!=NULL );
      return ADart->beta<B>();
    }
  };
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  // **************************************************************************
  /// Functor used to call update_dart_of_attribute<i>
  /// on each i-cell attribute enabled
  template<typename CMap>
  struct Update_dart_of_attribute_functor
  {
    template <unsigned int i>
    static void run(CMap* amap, typename CMap::Dart_handle ah, int amark)
    { amap->template update_dart_of_attribute<i>(ah,amark); }
  };

  template<typename CMap, unsigned int i, typename Enabled=
           typename CMap::Helper::
         #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
           template
         #endif
           Attribute_type<i>::type>
  struct Update_dart_of_one_attribute_functor
  {
    static void run(CMap* amap, typename CMap::Dart_handle ah, int amark)
    { amap->template update_dart_of_attribute<i>(ah,amark); }
  };
  template<typename CMap, unsigned int i>
  struct Update_dart_of_one_attribute_functor<CMap, i, CGAL::Void>
  {
    static void run(CMap*, typename CMap::Dart_handle, int)
    {}
  };
  // **************************************************************************
  // Functor used to degroup one i-attribute of one i-cell in two, except the
  // attribute of j.
  template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
           typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Degroup_attribute_functor_run
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    {
      CGAL_static_assertion( i<=CMap::dimension );
      CGAL_static_assertion( i!=j );
      CGAL_static_assertion_msg
          ( CMap::Helper::template Dimension_index<i>::value>=0,
           "Degroup_attribute_functor_run<i> but i-attributes are disabled" );
      CGAL_assertion( adart1!=NULL && adart2!=NULL );
      CGAL_assertion( adart1!=CMap::null_dart_handle &&
          adart2!=CMap::null_dart_handle );

      typename CMap::Helper::template Attribute_handle<i>::type
          a1=adart1->template attribute<i>();

      // If the two attributes are different, nothing to do.
      if ( a1!=adart2->template attribute<i>() || a1==NULL ) return;

      CGAL_assertion( (!belong_to_same_cell<CMap,i,CMap::dimension>
                      (*amap, adart1, adart2)) );

      typename CMap::Helper::template Attribute_handle<i>::type
          a2 = amap->template create_attribute<i>(*a1);

      amap->template set_attribute<i>(adart2, a2);
      Call_split_functor<CMap, i>::run(a1, a2);
    }
  };
  // Specialization for void attributes.
  template<typename CMap, unsigned int i, unsigned int j>
  struct Degroup_attribute_functor_run<CMap, i, j, CGAL::Void>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // Specialization for i==j.
  template<typename CMap, unsigned int i, typename T>
  struct Degroup_attribute_functor_run<CMap, i, i, T>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  // **************************************************************************
  /// Functor used by Test_split_attribute_functor_run to process one dart.
  template<typename CMap>
  struct Test_split_attribute_functor_one_dart
  {
    typedef typename CMap::Dart_handle Dart_handle;

    // Test the split of the i-cell containing the given dart adart.
    // When we process a dart, we search in the Unique_hash_map if its
    // i-attribute was already found. If yes, it means that we already
    // found an i-cell with this attribute, thus this attribute is split.
    // We mark (with mark) all the darts of the i-cell containing adart to
    // process them exactly once.
    template<unsigned int i>
    static void run( CMap* amap,
                     Dart_handle adart,
                     Unique_hash_map<typename CMap::Helper::template
                     Attribute_handle<i>::type, int> & found_attributes,
                     int mark )
    {
      CGAL_assertion( amap!=NULL );
      CGAL_assertion( adart!=NULL );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<i>::value>=0,
                                "Test_split_attribute_functor_one_dart<i> but "
                                "i-attributes are disabled");

      typedef typename CMap::Helper::template Attribute_handle<i>::type
      Attribute_handle_i;

      // If the current dart has no attribute, or if it is aldready marked,
      // nothing to do.
      if ( adart->template attribute<i>()==NULL ||
           amap->is_marked(adart, mark) )
        return;

      Attribute_handle_i a1 = adart->template attribute<i>();
      if ( found_attributes.is_defined(a1) )
      {  // Here the attribute was already present in the hash_map
        Attribute_handle_i a2 = amap->template create_attribute<i>(*a1);
        a2->set_dart(adart);

        for ( CMap_dart_iterator_basic_of_cell<CMap, i>
              itj(*amap, adart, mark); itj.cont(); ++itj )
        {
          amap->template set_attribute_of_dart<i>(itj, a2);
          amap->mark(itj, mark);
        }
        Call_split_functor<CMap, i>::run(a1, a2);
      }
      else
      {
        // Here the attribute was not in the hash_map.
        found_attributes[a1]=1;
        a1->set_dart(adart);

        for ( CMap_dart_iterator_basic_of_cell<CMap, i>
              itj(*amap, adart, mark); itj.cont(); ++itj )
        {
          CGAL_assertion( itj->template attribute<i>()==a1 );
          amap->mark(itj, mark);
        }
      }
    }
  };
  // **************************************************************************
  /// Functor used for unsew to test if i-attributes are split after an
  /// operation, except for j-attributes.
  ///   (j is the dimension of the beta modified for darts in modified_darts,
  ///    if j==0 modified_darts_2 are the darts modified for beta_1).
  template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
           typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Test_split_attribute_functor_run
  {
    static void run( CMap* amap,
                     std::deque<typename CMap::Dart_handle>
                     *modified_darts,
                     std::deque<typename CMap::Dart_handle>
                     */*modified_darts2*/,
                     int mark_modified_darts=-1)
    {
      CGAL_static_assertion( 1<=i && i<=CMap::dimension );
      CGAL_assertion( i!=j );
      CGAL_assertion( amap!=NULL );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<i>::value>=0,
                                "Test_split_attribute_functor_run<i> but "
                                "i-attributes are disabled");

      typedef typename CMap::Helper::template Attribute_handle<i>::type
      Attribute_handle_i;

      Unique_hash_map<Attribute_handle_i, unsigned int> found_attributes;

      int mark = amap->get_new_mark(); // to mark incident cells.
      typename std::deque<typename CMap::Dart_handle>::iterator
          it=modified_darts->begin();
      for ( ; it!=modified_darts->end(); ++it )
      {
        Test_split_attribute_functor_one_dart<CMap>::
            run<i>(amap, *it, found_attributes, mark);
      }

      // Now we unmark all the marked darts.
      amap->negate_mark(mark);
      for ( it=modified_darts->begin(); it!=modified_darts->end(); ++it )
      {
        if ( mark_modified_darts!=-1 )
          amap->unmark(*it, mark_modified_darts);

        if ( !amap->is_marked(*it, mark) )
          mark_cell<i>(*amap, *it, mark);
      }

      CGAL_assertion( amap->is_whole_map_marked(mark) );
      amap->free_mark(mark);
    }
  };
  // Specialization for i=0 and 2<=j.
  template<typename CMap, unsigned int j, typename T>
  struct Test_split_attribute_functor_run<CMap, 0, j, T>
  {
    static void run( CMap* amap,
                     std::deque<typename CMap::Dart_handle>
                     *modified_darts,
                     std::deque<typename CMap::Dart_handle>
                     */*modified_darts2*/,
                     int mark_modified_darts=-1)
    {
      CGAL_assertion( j!=0 && j!=1 );
      CGAL_assertion( amap!=NULL );
      CGAL_static_assertion_msg(CMap::Helper::template
                                Dimension_index<0>::value>=0,
                                "Test_split_attribute_functor_run<0> but "
                                "0-attributes are disabled");

      typedef typename CMap::Helper::template Attribute_handle<0>::type
          Attribute_handle_0;

      Unique_hash_map<Attribute_handle_0, unsigned int> found_attributes;
      typename CMap::Dart_handle od=NULL;

      int mark = amap->get_new_mark(); // to mark incident cells.
      typename std::deque<typename CMap::Dart_handle>::iterator
          it=modified_darts->begin();
      for ( ; it!=modified_darts->end(); ++it )
      {
        Test_split_attribute_functor_one_dart<CMap>::
            run<0>(amap, *it, found_attributes, mark);

        od=(*it)->other_extremity();
        if ( od!=NULL )
          Test_split_attribute_functor_one_dart<CMap>::
              run<0>(amap, od, found_attributes, mark);
      }

      // Now we unmark all the marked darts.
      amap->negate_mark(mark);
      for ( it=modified_darts->begin(); it!=modified_darts->end(); ++it )
      {
        if ( mark_modified_darts!=-1 )
          amap->unmark(*it, mark_modified_darts);

        if ( !amap->is_marked(*it, mark) )
          mark_cell<0>(*amap, *it, mark);

        od=(*it)->other_extremity();
        if ( od!=NULL && !amap->is_marked(od, mark) )
          mark_cell<0>(*amap, od, mark);
      }

      CGAL_assertion( amap->is_whole_map_marked(mark) );
      amap->free_mark(mark);
    }
  };
  // Specialization for i=0 and j=0.
  template<typename CMap, typename T>
  struct Test_split_attribute_functor_run<CMap, 0, 0, T>
  {
  };
  // Specialization for i=1 and j=0.
  template<typename CMap, typename T>
  struct Test_split_attribute_functor_run<CMap, 1, 0, T>
  {
    // No sense as we use always <0,0> ?
  };
  // Specialization for void attributes.
  template<typename CMap, unsigned int i, unsigned int j>
  struct Test_split_attribute_functor_run<CMap, i, j, CGAL::Void>
  {
  };
  /// Functor used for unsew to test if i-attributes are split after an
  /// operation, except for j-attributes.
  /// We define run<i> to allows to use this functor with
  /// Foreach_enabled_attributes.
  template<typename CMap, unsigned int j=CMap::dimension+1>
  struct Test_split_attribute_functor
  {
    typedef typename CMap::Dart_handle Dart_handle;

    // Test the split of the i-cell containing dart adart.
    template<unsigned int i>
    static void test_one_dart( CMap* amap,
                               Dart_handle adart,
                               Unique_hash_map<typename CMap::Helper::template
                               Attribute_handle<i>::type, int> &
                               found_attributes, int mark )
    {
      typedef typename CMap::Helper::template Attribute_handle<i>::type
      Attribute_handle_i;

      if ( adart->template attribute<i>()!=NULL &&
           !amap->is_marked(adart, mark) )
      {
        Attribute_handle_i a1 = adart->template attribute<i>();
        if ( found_attributes.is_defined(a1) )
        {  // Here the attribute was already present in the hash_map
          Attribute_handle_i a2 = amap->template create_attribute<i>(*a1);

          for ( CMap_dart_iterator_basic_of_cell<CMap, i>
                itj(*amap, adart, mark); itj.cont(); ++itj )
          {
            amap->template set_attribute_of_dart<i>(itj, a2);
            amap->mark(itj, mark);
          }
          a2->set_dart(adart);
          Call_split_functor<CMap, i>::run(a1, a2);
        }
        else
        {
          // Here the attribute was not in the hash_map.
          found_attributes[a1]=1;
          a1->set_dart(adart);

          for ( CMap_dart_iterator_basic_of_cell<CMap, i>
                itj(*amap, adart, mark); itj.cont(); ++itj )
          {
            CGAL_assertion( itj->template attribute<i>()==a1 );
            amap->mark(itj, mark);
          }
        }
      }
    }

    // Test the split of i-attributes, for all modified darts given in
    // modified_darts, and marked with mark_modified_darts.
    // For each split attribute, create a new i-attribute, associate
    // it with the new i-cell and call onsplit functors.
    template <unsigned int i>
    static void run( CMap* amap,
                     std::deque<typename CMap::Dart_handle>
                     *modified_darts,
                     int mark_modified_darts=-1)
    {
      CGAL_assertion( i!=j );
      CGAL_assertion( amap!=NULL );

      typedef typename CMap::Helper::template Attribute_handle<i>::type
      Attribute_handle_i;

      Unique_hash_map<Attribute_handle_i, unsigned int> found_attributes;

      int mark = amap->get_new_mark(); // to mark incident cells.
      for ( typename std::deque<Dart_handle>::iterator
            it=modified_darts->begin(); it!=modified_darts->end(); ++it )
      {
        test_one_dart<i>(amap, *it, found_attributes, mark);

        if ( /*j!=1 &&*/ i==0 ) // TODO Here we test too many cases, but more works to do
        {
          typename CMap::Dart_handle od = (*it)->other_extremity();
          if ( od!=NULL )
            test_one_dart<i>(amap, od, found_attributes, mark);
        }
      }

      // Now we unmark all the marked darts.
      amap->negate_mark(mark);
      for ( typename std::deque<Dart_handle>::iterator
            it=modified_darts->begin(); it!=modified_darts->end(); ++it )
      {
        if ( mark_modified_darts!=-1 )
          amap->unmark(*it, mark_modified_darts);

        if ( !amap->is_marked(*it, mark) )
          for ( CMap_dart_iterator_basic_of_cell<CMap,i>
                itj(*amap, *it, mark); itj.cont(); ++itj )
          {
            amap->mark(itj, mark);
          }

        if ( /*j!=1 &&*/ i==0 )
        {
          Dart_handle od = (*it)->other_extremity();
          if ( od!=NULL && !amap->is_marked(od, mark) )
            for ( CMap_dart_iterator_basic_of_cell<CMap,i>
                  itj(*amap, od, mark); itj.cont(); ++itj )
            {
              amap->mark(itj, mark);
            }
        }
      }

      CGAL_assertion( amap->is_whole_map_marked(mark) );
      amap->free_mark(mark);
    }
  };

  } // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTORS_H //
// EOF //
