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
#include <stack>
#include <set>

/* Definition of functors used to manage attributes (we need functors as
 * attributes are stored in tuple, thus all the access must be done at
 * compiling time.
 *
 * Group_attribute_functor_of_dart<CMap> to group the <i>-attributes of two
 *    given darts (except for j-dim). Only the attributes of the two given
 *    darts are possibly modified.
 * Group_attribute_functor_of_dart_run<CMap,i> same than
 *   Group_attribute_functor_of_dart<CMap>::run<i>, with i template argument
 *   given in the struct to enable specialization.
 *
 * Group_attribute_functor<CMap> to group the <i>-attributes of two
 *    given i-cells (except forj-adim). If one i-attribute is NULL, we set the
 *    darts of its i-cell to the second attribute. If both i-attributes are
 *    non NULL, we overide all the i-attribute of the second i-cell to the
 *    first i-attribute.
 *
 * Call_split_functor<CMap,i> to call the OnSplit functor on two given
 *    i-attributes.
 *
 ** Call_merge_functor<CMap,i> to call the OnMerge functor on two given
 *    i-attributes.
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

  // Functor used to call the On_split functor between the two given darts.
  template<typename Map,unsigned int i,
           typename Enabled=typename Map::Helper::
         #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
           template
         #endif
           Attribute_type<i>::type>
  struct Call_split_functor
  {
    static void run(typename Map::Dart_handle adart1,
                    typename Map::Dart_handle adart2)
    {
      Apply_cell_functor
          <typename Map::Helper::template Attribute_type<i>::type,
          typename Map::Helper::template Attribute_type<i>::type::On_split>::
          run(*(adart1->template attribute<i>()),
              *(adart2->template attribute<i>()));
    }
    static void
    run(typename Map::Helper::template Attribute_handle<i>::type a1,
        typename Map::Helper::template Attribute_handle<i>::type a2)
    {
      Apply_cell_functor
          <typename Map::Helper::template Attribute_type<i>::type,
          typename Map::Helper::template Attribute_type<i>::type::On_split>::
          run(*a1, *a2);
    }
  };

  // Specialization for disabled attributes.
  template<typename Map,unsigned int i>
  struct Call_split_functor<Map,i,CGAL::Void>
  {
    static void run(typename Map::Dart_handle,
                    typename Map::Dart_handle)
    {}
  };

  // Functor used to call the On_merge functor between the two given darts.
  template<typename Map,unsigned int i,
           typename Enabled=typename Map::Helper::
         #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
           template
         #endif
           Attribute_type<i>::type>
  struct Call_merge_functor
  {
    static void run(typename Map::Dart_handle adart1,
                    typename Map::Dart_handle adart2)
    {
      Apply_cell_functor
          <typename Map::Helper::template Attribute_type<i>::type,
          typename Map::Helper::template Attribute_type<i>::type::On_merge>::
          run(*(adart1->template attribute<i>()),
              *(adart2->template attribute<i>()));
    }
    static void
    run(typename Map::Helper::template Attribute_handle<i>::type a1,
        typename Map::Helper::template Attribute_handle<i>::type a2)
    {
      Apply_cell_functor
          <typename Map::Helper::template Attribute_type<i>::type,
          typename Map::Helper::template Attribute_type<i>::type::On_merge>::
          run(*a1, *a2);
    }
  };

  // Specialization for disabled attributes.
  template<typename Map,unsigned int i>
  struct Call_merge_functor<Map,i,CGAL::Void>
  {
    static void run(typename Map::Dart_handle,
                    typename Map::Dart_handle)
    {}
  };

  // **************************************************************************
  /// Functor used for link_beta to update the i-attributes of
  /// adart2 on the attributes of this dart, except if i=j.
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
      CGAL_static_assertion(i<=CMap::dimension);
      CGAL_static_assertion(i!=j);
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
  template<typename CMap, unsigned int i, unsigned int j>
  struct Group_attribute_functor_of_dart_run<CMap,i,j,CGAL::Void>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };
  template<typename CMap, unsigned int i, typename T>
  struct Group_attribute_functor_of_dart_run<CMap,i,i,T>
  {
    static void run(CMap*,
                    typename CMap::Dart_handle,
                    typename CMap::Dart_handle)
    {}
  };

  /// Functor used for link_beta to update the attributes of
  /// adart2 on the attributes of this dart, except for j-attributes.
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
  template<typename CMap, unsigned int i, unsigned int j=CMap::dimension+1,
           typename T=
           typename CMap::Helper::template Attribute_type<i>::type>
  struct Group_attribute_functor_run
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2)
    {
      CGAL_static_assertion( i<=CMap::dimension );
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
            { // TODO REMOVE ?
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
            { // TODO REMOVE ?
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
    }
  };
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
            { // TODO REMOVE ?
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
    }
  };
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
            { // TODO REMOVE ?
              Call_merge_functor<CMap, 0>::run(a1, a2);
            }
          }
          // TODO set_attribute templated by a range
          amap->template set_attribute<0>(toSet, a1);
        }
      }
    }
  };
  template<typename CMap, unsigned int i, unsigned int j>
  struct Group_attribute_functor_run<CMap, i, j, CGAL::Void>
  {
    static void run( CMap*,
                     typename CMap::Dart_handle,
                     typename CMap::Dart_handle )
    {}
  };

  /// Functor used for sew to update the attributes of
  /// adart2 on the attributes of this dart, except for j-attributes.
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
    {
      Group_attribute_functor_run<CMap,i,j>::run(amap,adart1,adart2);
    }
  };

  // Functor used to degroup the two n-attributes of the two darts, except the
  // attribute of adim
  template<typename CMap,unsigned int i>
  struct Degroup_attribute_functor_run
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2, int adim)
    {
      CGAL_assertion( adim==-1 ||
                      (1<=adim && (unsigned int)adim<=CMap::dimension) );
      if (i!=adim )
      {
        amap->template degroup_enabled_attribute
          <i, typename CMap::Helper::template Attribute_type<i>::type>
          (adart1, adart2);
      }
    }
  };
  template<typename CMap>
  struct Degroup_attribute_functor_run<CMap, 0>
  {
    static void run(CMap* amap,
                    typename CMap::Dart_handle adart1,
                    typename CMap::Dart_handle adart2, int adim)
    {
      CGAL_assertion( adim==-1 ||
                      (1<=adim && (unsigned int)adim<=CMap::dimension) );
      typename CMap::Dart_handle od = adart1->other_extremity();
      if ( od!=NULL )
        amap->template degroup_enabled_attribute
          <0, typename CMap::Helper::template Attribute_type<0>::type >
          (od, adart2);

      if ( adim!=1 )
      {
        od = adart2->other_extremity();
        if ( od!=NULL )
        {
          amap->template degroup_enabled_attribute
            <0, typename CMap::Helper::template Attribute_type<0>::type>
            (adart1, od);
        }
      }
    }
  };
  template<typename Map>
  struct Degroup_attribute_functor
  {
    template <unsigned int i>
    static void run(Map* amap,typename Map::Dart_handle adart1,
                    typename Map::Dart_handle adart2, int adim)
    {
      Degroup_attribute_functor_run<Map,i>::run(amap,adart1,adart2,adim);
    }
  };


  ////////////////////////////////////
    template<typename Dart_handle>
    struct Couple_dart_and_dim
    {
      Couple_dart_and_dim(Dart_handle ad1,Dart_handle ad2,int adim) :
        d1(ad1), d2(ad2), dim(adim)
      {}
      Dart_handle d1,d2;
      int dim;
    };

    // Functor used to degroup one attribute of two given darts
    template <typename CMap, unsigned int i, typename Type_attr>
    struct Degroup_one_attribute_functor
    {
      static bool run(CMap* amap,
                      typename CMap::Dart_handle adart1,
                      typename CMap::Dart_handle adart2)
      {
        CGAL_assertion(amap!=NULL);
        return amap->template degroup_enabled_attribute<i, Type_attr>
          (adart1,adart2);
      }
    };

    // Specialization for i-attributes disabled.
    template <typename CMap, unsigned int i>
    struct Degroup_one_attribute_functor<CMap,i,CGAL::Void>
    {
      static bool run(CMap*,
                      typename CMap::Dart_handle,
                      typename CMap::Dart_handle)
      { return false; }
    };

/*    // Functor used to degroup one attribute of one dart
    template <typename CMap, unsigned int i, typename Type_attr, typename Range>
    struct Degroup_one_attribute_of_dart_functor
    {
      static bool run(CMap* amap,
                      typename CMap::Dart_handle adart1,
                      typename CMap::Dart_handle adart2)
      {
        CGAL_assertion(amap!=NULL);
        return amap->template degroup_enabled_attribute_of_dart
          <i, Type_attr, Range>(adart1,adart2);
      }
    };

    // Specialization for i-attributes disabled.
    template <typename CMap, unsigned int i, typename Range>
    struct Degroup_one_attribute_of_dart_functor<CMap, i, CGAL::Void, Range>
    {
      static bool run(CMap*,
                      typename CMap::Dart_handle,
                      typename CMap::Dart_handle)
      { return false; }
    };*/

    /// Functor used to call decrease_attribute_ref_counting<i>
    /// on each i-cell attribute enabled
    template<typename Map>
    struct Decrease_attribute_functor
    {
      template <unsigned int i>
      static void run(Map* amap, typename Map::Dart_handle adart)
      { amap->template
          decrease_attribute_ref_counting<i>(adart/*,Tag_true()*/); }
    };

    /// Functor used to call update_dart_of_attribute<i>
    /// on each i-cell attribute enabled
    template<typename Map>
    struct Update_dart_of_attribute_functor
    {
      template <unsigned int i>
      static void run(Map* amap, typename Map::Dart_handle ah, int amark)
      { amap->template update_dart_of_attribute<i>(ah,amark); }
    };

    template<typename Map, unsigned int i, typename Enabled=
             typename Map::Helper::
         #ifndef CGAL_CFG_TEMPLATE_IN_DEFAULT_PARAMETER_BUG
             template
         #endif
             Attribute_type<i>::type>
    struct Update_dart_of_one_attribute_functor
    {
      static void run(Map* amap, typename Map::Dart_handle ah, int amark)
      { amap->template update_dart_of_attribute<i>(ah,amark); }
    };
    template<typename Map, unsigned int i>
    struct Update_dart_of_one_attribute_functor<Map, i, CGAL::Void>
    {
      static void run(Map*, typename Map::Dart_handle, int)
      {}
    };

    /// Functor used to reserve one mark for each enabled attribute.
    template<typename Map>
    struct Reserve_mark_functor
    {
      template <unsigned int i>
      static void run(const Map* amap, std::vector<int>* marks)
      { (*marks)[i] = amap->get_new_mark(); }
    };

    /// Functor used to test if a cell is valid
    template<typename Map>
    struct Test_is_valid_attribute_functor
    {
      template <unsigned int i>
      static void run(const Map* amap,
                      typename Map::Dart_const_handle adart,
                      std::vector<int>* marks, bool *ares)
      {
        if (!amap->template is_valid_attribute<i>(adart,(*marks)[i]) )
        {
          (*ares)=false;
          std::cerr << "Map not valid: a "<<i<<"-cell is not correctly "
            "associated with an attribute for " << &(*adart)<< std::endl;
        }
      }
    };

    /// Functor for counting i-cell
    template<typename Map>
    struct Count_cell_functor
    {
      template <unsigned int i>
      static void run( const Map* amap,
                       typename Map::Dart_const_handle adart,
                       std::vector<int>* amarks,
                       std::vector<unsigned int>* ares )
      {
        if ( (*amarks)[i]!=-1 && !amap->is_marked(adart, (*amarks)[i]) )
        {
          ++ (*ares)[i];
          mark_cell<Map,i>(*amap, adart, (*amarks)[i]);
        }
      }
    };


    /// Functor for counting the memory occupation of attributes
    /// Be careful not reentrant !!! TODO a  Foreach_enabled_attributes
    /// taking an instance of a functor as argument allowing to compute
    /// and return values.
    template<typename Map>
    struct Count_bytes_one_attribute_functor
    {
      template <unsigned int i>
      static void run( const Map* amap )
      {
        res += amap->template attributes<i>().capacity()*
          sizeof(typename Map::template Attribute_type<i>::type);
      }

      static typename Map::size_type res;
    };
    template<typename Map>
    typename Map::size_type Count_bytes_one_attribute_functor<Map>::res = 0;

    template<typename Map>
    struct Count_bytes_all_attributes_functor
    {
      static typename Map::size_type run( const Map& amap )
      {
        Count_bytes_one_attribute_functor<Map>::res = 0;
        Map::Helper::template Foreach_enabled_attributes
          <Count_bytes_one_attribute_functor<Map> >::run(&amap);
        return Count_bytes_one_attribute_functor<Map>::res;
      }
    };

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

    template<typename Map, unsigned int i>
    struct Store_incident_cells
    {
      template <unsigned int j>
      static void run( Map* amap, typename Map::Dart_handle adart,
                       int  mark_for_icell,
                       int* mark_for_incident_cells,
                       std::deque<std::deque<typename Map::Dart_handle> >
                        *store )
      {
        if ( i==j ) return;

        const int mark_for_jcells = mark_for_incident_cells
            [Map::Helper::template Dimension_index<j>::value];

        std::deque<std::deque<typename Map::Dart_handle> >& jcells =
            store[Map::Helper::template Dimension_index<j>::value];

        CGAL_assertion( amap!=NULL );
        CGAL_assertion( adart!=NULL );
        CGAL_assertion( amap->is_reserved(mark_for_icell) );
        CGAL_assertion( amap->is_reserved(mark_for_jcells) );

        if ( !amap->is_marked(adart, mark_for_jcells) &&
             adart->template attribute<j>()!=NULL )
        {
          jcells.push_back(std::deque<typename Map::Dart_handle>());
          for ( CMap_dart_iterator_basic_of_cell<Map,j>
                itj(*amap, adart, mark_for_jcells); itj.cont(); ++itj )
          {
            if ( !amap->is_marked(itj, mark_for_icell) )
            {
              jcells.back().push_back(itj);
            }
            amap->mark(itj, mark_for_jcells);
          }
          if ( jcells.back().empty() ) jcells.pop_back();
        }

        if ( i!=1 && j==0 )
        {
          typename Map::Dart_handle od = adart->other_extremity();

          if ( od!=NULL && !amap->is_marked(od, mark_for_jcells) &&
               od->template attribute<j>()!=NULL )
          {
            jcells.push_back(std::deque<typename Map::Dart_handle>());
            for ( CMap_dart_iterator_basic_of_cell<Map,j>
                  itj(*amap, od, mark_for_jcells); itj.cont(); ++itj )
            {
              if ( !amap->is_marked(itj, mark_for_icell) )
              {
                jcells.back().push_back(itj);
              }
              amap->mark(itj, mark_for_jcells);
            }
            if ( jcells.back().empty() ) jcells.pop_back();
          }
        }
      }
    };

    template<typename Map, unsigned int i>
    struct Test_split_with_deque
    {
      template <unsigned int j>
      static void run( Map* amap,
                       int* mark_for_incident_cells,
                       std::deque<std::deque<typename Map::Dart_handle> >
                        *store )
      {
        const int mark_for_jcells = mark_for_incident_cells
            [Map::Helper::template Dimension_index<j>::value];
        amap->negate_mark( mark_for_jcells );

        if ( i==j ) return;

        std::deque<std::deque<typename Map::Dart_handle> >& jcells =
            store[Map::Helper::template Dimension_index<j>::value];

        CGAL_assertion( amap!=NULL );
        CGAL_assertion( amap->is_reserved(mark_for_jcells) );

        int nbofjcell = 0;
        typename Map::Helper::template Attribute_handle<j>::type
            a1 = NULL;
        typename Map::Helper::template Attribute_handle<j>::type
            a2=NULL;

        int nb=0;
        for ( typename std::deque<std::deque<typename Map::Dart_handle> >::
              iterator it=jcells.begin(); it!=jcells.end(); ++it )
        {
          nbofjcell = 0;
          for ( typename std::deque<typename Map::Dart_handle>::iterator
                itj=it->begin(); itj!=it->end(); ++itj )
          {
            if ( !amap->is_marked(*itj, mark_for_jcells) )
            {
              ++nbofjcell;
              if ( nbofjcell>1 )
              {
                a2 = amap->template create_attribute<j>(*a1);
                // std::cout<<"A2 "<<&*a2<<"  "<<&**itj<<": ";
                // We call the on_split functor
              }
              else
              {
                a1=(*itj)->template attribute<j>();
                a1->set_dart(*itj);
                // std::cout<<"A1 "<<&*a1<<"  "<<&**itj<<": ";
              }

              for ( CMap_dart_iterator_basic_of_cell<Map,j>
                  itj2(*amap, *itj, mark_for_jcells);
                  itj2.cont(); ++itj2 )
              {
                // std::cout<<&*itj2<<", ";
                if ( nbofjcell>1 )
                  amap->template set_attribute_of_dart<j>(itj2, a2);
                ++nb;
                amap->mark(itj2, mark_for_jcells);
              }
              // std::cout<<std::endl;

              if ( nbofjcell>1 )
                Apply_cell_functor
                    <typename Map::Helper::template Attribute_type<j>::type,
                    typename Map::Helper::template Attribute_type<j>::type::
                    On_split>::run(*a1, *a2);
            }
          }
        }
        //std::cout<<"number of marked darts for <"<<j<<"> : "<<amap->number_of_marked_darts(mark_for_jcells)<<std::endl;
        //std::cout<<"number of iterated darts : "<<nb<<std::endl;

      }
    };

    template<typename Map, unsigned int i>
    struct Test2_split_with_deque
    {
      template<unsigned int j>
      static void test_one_dart( Map* amap,
                                 typename Map::Dart_handle adart,
                                 std::set<typename Map::Helper::template
                                 Attribute_handle<j>::type>& found_attributes,
                                 int mark)
      {
        if ( adart->template attribute<j>()!=NULL &&
             !amap->is_marked(adart, mark) )
        {
          typename Map::Helper::template Attribute_handle<j>::type
              a1 = adart->template attribute<j>();
          if ( !found_attributes.insert(a1).second )
          {  // Here the attribute was already present in the set
            typename Map::Helper::template Attribute_handle<j>::type
                a2 = amap->template create_attribute<j>(*a1);
            // std::cout<<"A2 "<<&*a2<<"  "<<&**itj<<": ";

            for ( CMap_dart_iterator_basic_of_cell<Map,j>
                  itj(*amap, adart, mark);
                  itj.cont(); ++itj )
            {
              // std::cout<<&*itj<<", ";
              amap->template set_attribute_of_dart<j>(itj, a2);
              amap->mark(itj, mark);
            }

            Apply_cell_functor
                <typename Map::Helper::template Attribute_type<j>::type,
                typename Map::Helper::template Attribute_type<j>::type::
                On_split>::run(*a1, *a2);
          }
          else
          {
            // Here the attribute was not in the set as we are able
            // to insert it.
            a1->set_dart(adart);
            // std::cout<<"A1 "<<&*a1<<"  "<<&**itj<<": ";

            for ( CMap_dart_iterator_basic_of_cell<Map,j>
                  itj(*amap, adart, mark);
                  itj.cont(); ++itj )
            {
              // std::cout<<&*itj<<", ";
              CGAL_assertion( itj->template attribute<j>()==a1 );
              amap->mark(itj, mark);
            }
          }
          // std::cout<<std::endl;
        }
      }

      template <unsigned int j>
      static void run( Map* amap,
                       std::deque<typename Map::Dart_handle>
                        *modified_darts,
                       int mark_modified_darts
                       /*,
                       std::deque<typename Map::Dart_handle>
                        *modified_darts2*/)
      {
        if ( i==j ) return;

        CGAL_assertion( amap!=NULL );

        std::set<typename Map::Helper::template
            Attribute_handle<j>::type> found_attributes;

        int mark = amap->get_new_mark();
        int nb=0;
        for ( typename std::deque<typename Map::Dart_handle>::
              iterator it=modified_darts->begin();
              it!=modified_darts->end(); ++it )
        {
          test_one_dart<j>(amap, *it, found_attributes, mark);

          if ( i!=1 && j==0 )
          {
            typename Map::Dart_handle od = (*it)->other_extremity();
            if ( od!=NULL )
              test_one_dart<j>(amap, od, found_attributes, mark);
          }
        }

/*        if ( i+1==j )
        {
          for ( typename std::deque<typename Map::Dart_handle>::
                iterator it=modified_darts2->begin();
                it!=modified_darts2->end(); ++it )
          {
            test_one_dart<j>(amap, *it, found_attributes, mark);

            if ( i!=1 && j==0 )
            {
              typename Map::Dart_handle od = (*it)->other_extremity();
              if ( od!=NULL )
                test_one_dart<j>(amap, od, found_attributes, mark);
            }
          }
        }*/
        // std::cout<<"Size of set for dim "<<j<<": "<<found_attributes.size()<<std::endl;
        // std::cout<<"Size of deque for dim "<<j<<": "<<modified_darts->size()<<std::endl;
        // std::cout<<"number of marked darts : "<<amap->number_of_marked_darts(mark)<<std::endl;
        // std::cout<<"number of iterated darts : "<<nb<<std::endl;

        // Now we unmark all the marked darts.
        amap->negate_mark(mark);

        for ( typename std::deque<typename Map::Dart_handle>::
              iterator it=modified_darts->begin();
              it!=modified_darts->end(); ++it )
        {
          if ( mark_modified_darts!=-1 )
            amap->unmark(*it, mark_modified_darts);

          if ( !amap->is_marked(*it, mark) )
            for ( CMap_dart_iterator_basic_of_cell<Map,j>
                itj(*amap, *it, mark);
                itj.cont(); ++itj )
            {
              amap->mark(itj, mark);
            }

          if ( i!=1 && j==0 )
          {
            typename Map::Dart_handle od = (*it)->other_extremity();
            if ( od!=NULL && !amap->is_marked(od, mark) )
              for ( CMap_dart_iterator_basic_of_cell<Map,j>
                  itj(*amap, od, mark);
                    itj.cont(); ++itj )
              {
                amap->mark(itj, mark);
              }
          }
        }

 /*       if ( i+1==j )
        {
          for ( typename std::deque<typename Map::Dart_handle>::
                iterator it=modified_darts2->begin();
                it!=modified_darts2->end(); ++it )
          {
            if ( mark_modified_darts!=-1 )
              amap->unmark(*it, mark_modified_darts);

            if ( !amap->is_marked(*it, mark) )
              for ( CMap_dart_iterator_basic_of_cell<Map,j>
                  itj(*amap, *it, mark);
                  itj.cont(); ++itj )
              {
                amap->mark(itj, mark);
              }

            if ( i!=1 && j==0 )
            {
              typename Map::Dart_handle od = (*it)->other_extremity();
              if ( od!=NULL && !amap->is_marked(od, mark) )
                for ( CMap_dart_iterator_basic_of_cell<Map,j>
                    itj(*amap, od, mark);
                      itj.cont(); ++itj )
                {
                  amap->mark(itj, mark);
                }
            }
          }
        }*/

        CGAL_assertion( amap->is_whole_map_marked(mark) );
        amap->free_mark(mark);
      }
    };

  } // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTORS_H //
// EOF //
