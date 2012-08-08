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

namespace CGAL {
  namespace internal {

    /** @file Combinatorial_map_functors.h
     * Definition of functors used for dD Combinatorial map.
     */
    
    template<typename Dart_handle>
    struct Couple_dart_and_dim
    {
      Couple_dart_and_dim(Dart_handle ad1,Dart_handle ad2,int adim) :
        d1(ad1), d2(ad2), dim(adim)
      {}
      Dart_handle d1,d2;
      int dim;
    };

    // Functor used to group one attribute of two given darts
    template <typename CMap, unsigned int i, typename Type_attr>
    struct Group_one_attribute_functor
    {
      static void run(CMap* amap,
                      typename CMap::Dart_handle adart1,
                      typename CMap::Dart_handle adart2)
      { 
        CGAL_assertion(amap!=NULL);
        amap->template group_enabled_attribute<i, Type_attr>(adart1,adart2);
      }
    };

    // Specialization for i-attributes disabled.
    template <typename CMap, unsigned int i>
    struct Group_one_attribute_functor<CMap,i,CGAL::Void>
    {
      static void run(CMap*,
                      typename CMap::Dart_handle,
                      typename CMap::Dart_handle)
      {}
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

    // Functor used to degroup one attribute of one dart
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
    };

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

    // Functor used to group the two n-attributes of the two darts, except the
    // attribute of adim (adim==-1 || 1<=adim<=dimension)
    template<typename Map,unsigned int i>
    struct Group_attribute_functor_run
    {
      static void run(Map* amap,
                      typename Map::Dart_handle adart1, 
                      typename Map::Dart_handle adart2, int adim)
      {
        if ( i!=adim )
        {
          amap->template group_enabled_attribute
            <i, typename Map::Helper::template Attribute_type<i>::type>
            (adart1, adart2);
        }    
      }
    };

    template<typename CMap>
    struct Group_attribute_functor_run<CMap,0>
    {
      static void run(CMap* amap,
                      typename CMap::Dart_handle adart1, 
                      typename CMap::Dart_handle adart2, int adim)
      {
        typename CMap::Dart_handle od = adart1->other_extremity();
        if ( od!=NULL )
        {
          amap->template group_enabled_attribute
            <0, typename CMap::Helper::template Attribute_type<0>::type>
            (od, adart2);
        }
      
        if ( adim!=1 )
        {
          od = adart2->other_extremity();
          if ( od!=NULL )
            amap->template group_enabled_attribute
              <0, typename CMap::Helper::template Attribute_type<0>::type>
              (adart1, od);
        }
      }
    };

    template<typename Map>
    struct Group_attribute_functor
    {
      template <unsigned int i>
      static void run(Map* amap,
                      typename Map::Dart_handle adart1, 
                      typename Map::Dart_handle adart2, int adim)
      {
        CGAL_assertion( adim==-1 || 
                        (1<=adim && (unsigned int)adim<=Map::dimension) );
        Group_attribute_functor_run<Map,i>::run(amap,adart1,adart2,adim);
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
    
    // Functor which call operator() on the cell_attribute...
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

    /// Functor used for link_beta to update the attributes of 
    /// adart2 on the attributes of this dart, except for adimension-attributes.
    template<typename CMap, unsigned int i>
    struct Group_attribute_functor_of_dart_run
    {
      static void run(CMap* amap,
                      typename CMap::Dart_handle dh1,
                      typename CMap::Dart_handle dh2,
                      int adim)
      {
        CGAL_assertion( adim==-1 || 
                        (0<=adim && (unsigned int)adim<=CMap::dimension) );
        if ( adim!=i ) 
        {
          amap->template group_enabled_attribute_of_dart
            <i, typename CMap::Helper::template Attribute_type<i>::type>
            (dh1, dh2);
        }
      }
    };
    template<typename CMap>
    struct Group_attribute_functor_of_dart_run<CMap, 0>
    {
      static void run(CMap* amap,
                      typename CMap::Dart_handle dh1,
                      typename CMap::Dart_handle dh2,
                      int adim)
      {
        CGAL_assertion( adim==-1 || 
                        (0<=adim && (unsigned int)adim<=CMap::dimension) );
        // todo ASSERT (1<=adim && ...) ???
        if ( adim!=0 )
        {
          typename CMap::Dart_handle od = dh1->other_extremity();
          if ( od!=NULL )
          {
            typename CMap::Helper::template  Attribute_handle<0>::type
              a1=od->template attribute<0>();
            if ( a1!=NULL && a1!=dh2->template attribute<0>() )
              amap->template set_attribute_of_dart<0>(dh2, a1);
          }
        }

        if ( adim!=1 )
        {
          typename CMap::Dart_handle od = dh2->other_extremity();
          if ( od!=NULL && dh1->template attribute<0>()==NULL )
          {
            typename CMap::Helper::template Attribute_handle<0>::type
              a2=od->template attribute<0>();
            if ( a2!=NULL )
              amap->template set_attribute_of_dart<0>(dh1, a2);
          }            
        }          
      }
    };
      
    template<typename CMap>
    struct Group_attribute_functor_of_dart
    {
      template <unsigned int i>
      static void run(CMap* amap,
                      typename CMap::Dart_handle adart1, 
                      typename CMap::Dart_handle adart2, int adim)
      {
        CGAL_assertion( adim==-1 || 
                        (0<=adim && (unsigned int)adim<=CMap::dimension) );
        Group_attribute_functor_of_dart_run<CMap,i>::
          run(amap,adart1,adart2,adim);
      }
    };

    // Functor used to call the On_split functor between the two given darts.
    template<typename Map,unsigned int i,
             typename Enabled=
             typename Map::Helper::
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
    };

    // Specialization for disabled attributes.
    template<typename Map,unsigned int i>
    struct Call_split_functor<Map,i,CGAL::Void>
    {
      static void run(typename Map::Dart_handle, 
                      typename Map::Dart_handle)
      {}    
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
      { return Beta_functor<Dart_handle, Betas...>::run(ADart->beta(B), betas...); }
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
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
 
  } // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTORS_H //
// EOF //
