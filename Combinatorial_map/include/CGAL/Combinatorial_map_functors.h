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
#include <vector>

/* Definition of functors used to manage attributes (we need functors as
 * attributes are stored in tuple, thus all the access must be done at
 * compiling time). Some of these functors are used with
 * Foreach_enabled_attributes to iterate through all the non void attribs.
 * Functors allowing to group/ungroup attributes are defined in
 * Combinatorial_map_group_functors.h (included at the end of this file)
 *
 * Reserve_mark_functor<CMap> to reserve one mark for each non void attribute.
 *
 * Set_i_attribute_functor<CMap, i> to set the i-attribute of a given
 *   i-cell.
 *
 * internal::Call_split_functor<CMap,i> to call the OnSplit functors on two
 *    given i-attributes.
 *
 * internal::Call_merge_functor<CMap,i> to call the OnMerge functors on two
 *    given i-attributes.
 *
 * internal::Test_is_valid_attribute_functor<CMap> to test if a given i-cell is
 *    valid (all its darts are linked to the same attribute, no other dart is
 *    linked with this attribute).
 *
 * internal::Count_cell_functor<CMap> to count the nuber of i-cells.
 *
 * internal::Count_bytes_one_attribute_functor<CMap> to count the memory
 *    occupied by i-attributes.
 *
 * internal::Decrease_attribute_functor<CMap> to decrease by one the ref
 *    counting of a given i-attribute.
 *
 * internal::Beta_functor<Dart, i...> to call several beta on the given dart.
 *   Indices are given as parameter of the run function.
 *
 * internal::Beta_functor_static<Dart, i...> to call several beta on the given
 *   dart. Indices are given as template arguments.
 *
 * internal::Set_i_attribute_of_dart_functor<CMap, i> to set the i-attribute
 *   of a given dart.
 */

namespace CGAL
{
  /** @file Combinatorial_map_functors.h
   * Definition of functors used for dD Combinatorial map.
   */
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
  /// Functor used to display the address of the i-cell attribute
  template<typename CMap>
  struct Display_attribute_functor
  {
    template <unsigned int i>
    static void run(const CMap* amap,
                    typename CMap::Dart_const_handle adart)
    {
      if ( adart->template attribute<i>()==NULL )
        std::cout<<"NULL";
      else
        std::cout<<&*(adart->template attribute<i>());
    }
  };
  // **************************************************************************
  namespace internal
  {
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
  /// Decrease the cell attribute reference counting of the given dart.
  /// The attribute is removed if there is no more darts linked with it.
  template<typename CMap, unsigned int i, typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Decrease_attribute_functor_run
  {
    static void run(CMap* amap, typename CMap::Dart_handle adart)
    {
      if ( adart->template attribute<i>()!=NULL )
      {
        adart->template attribute<i>()->dec_nb_refs();
        if ( adart->template attribute<i>()->get_nb_refs()==0 )
          amap->template erase_attribute<i>(adart->template attribute<i>());
      }
    }
  };
  /// Specialization for void attributes.
  template<typename CMap, unsigned int i>
  struct Decrease_attribute_functor_run<CMap,i,CGAL::Void>
  {
    static void run(CMap*, typename CMap::Dart_handle)
    {}
  };
  // **************************************************************************
  /// Functor used to call decrease_attribute_ref_counting<i>
  /// on each i-cell attribute enabled
  template<typename CMap>
  struct Decrease_attribute_functor
  {
    template <unsigned int i>
    static void run(CMap* amap, typename CMap::Dart_handle adart)
    { Decrease_attribute_functor_run<CMap,i>::run(amap, adart); }
  };
  // **************************************************************************
  /// Functor used to set the i-attribute of a given dart.
  template<typename CMap, unsigned int i, typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Set_i_attribute_of_dart_functor
  {
    static void run( CMap* amap, typename CMap::Dart_handle dh,
                     typename CMap::Helper::template Attribute_handle<i>::type
                     ah )
    {
      CGAL_static_assertion(i<=CMap::dimension);
      CGAL_assertion( dh!=NULL && dh!=CMap::null_dart_handle );

      if ( dh->template attribute<i>()==ah ) return;

      Decrease_attribute_functor_run<CMap, i>::run(amap, dh);
      dh->template set_attribute<i>(ah);
    }
  };
  /// Specialization for void attributes.
  template<typename CMap, unsigned int i>
  struct Set_i_attribute_of_dart_functor<CMap,i,CGAL::Void>
  {
    static void run( CMap*, typename CMap::Dart_handle,
                     typename CMap::Helper::template Attribute_handle<i>::type)
    {}
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
  template<typename Dart_handle, int ... Betas>
  struct Beta_functor_static;
  template<typename Dart_handle, int B, int ... Betas>
  struct Beta_functor_static<Dart_handle, B, Betas...>
  {
    static Dart_handle run(Dart_handle ADart)
    { return Beta_functor_static<Dart_handle, Betas...>::
          run(ADart->template beta<B>()); }
  };
  template<typename Dart_handle, int B>
  struct Beta_functor_static<Dart_handle, B>
  {
    static Dart_handle run(Dart_handle ADart)
    {
      CGAL_assertion( ADart!=NULL );
      return ADart->template beta<B>();
    }
  };
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  // **************************************************************************
  /// Functor used to call update_dart_of_attribute<i>
  /// on each i-cell attribute enabled
  // TODO REMOVE ?
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
  } // namespace internal
  // **************************************************************************
  /// Functor used to set the i-attribute of a given i-cell.
  /// We can use any range as Range type, by default we use
  /// Dart_of_cell_range<i>
  template<typename CMap, unsigned int i,
           typename Range=typename CMap::template Dart_of_cell_range<i>,
           typename T=typename CMap::Helper::template Attribute_type<i>::type>
  struct Set_i_attribute_functor
  {
    static void run( CMap* amap, typename CMap::Dart_handle dh,
                     typename CMap::Helper::template Attribute_handle<i>::type
                     ah )
    {
      CGAL_static_assertion(i<=CMap::dimension);
      CGAL_assertion( dh!=NULL && dh!=CMap::null_dart_handle && ah!=NULL );

      for ( typename Range::iterator it(*amap, dh); it.cont(); ++it)
      {
        if ( it->template attribute<i>()!=ah )
        {
          internal::Decrease_attribute_functor_run<CMap, i>::run(amap, it);
          it->template set_attribute<i>(ah);
        }
      }
      ah->set_dart(dh);
    }
  };
  /// Specialization for void attributes.
  template<typename CMap, unsigned int i>
  struct Set_i_attribute_functor<CMap,i,CGAL::Void>
  {
    static void run( CMap*, typename CMap::Dart_handle,
                     typename CMap::Helper::template Attribute_handle<i>::type)
    {}
  };
  // **************************************************************************
} // namespace CGAL

#include <CGAL/internal/Combinatorial_map_group_functors.h>

#endif // CGAL_COMBINATORIAL_MAP_FUNCTORS_H //
// EOF //
