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

    
    /// Functor used to test if it is possible to i-sew two darts,
    /// 2<=i<=dimension
    template <typename Map,unsigned int i>
    struct is_sewable_functor{
      static bool run(const Map& amap,
                      typename Map::Dart_const_handle adart1,
                      typename Map::Dart_const_handle adart2)
      {
        CGAL_static_assertion(2<=i && i<=Map::dimension);
        CGAL_assertion(adart1!=NULL && adart2!=NULL);
      
        if ( !adart1->is_free(i) || !adart2->is_free(i) || adart1==adart2 )
          return false;
      
        CMap_dart_const_iterator_of_involution<Map,i>     I1(amap, adart1);
        CMap_dart_const_iterator_of_involution_inv<Map,i> I2(amap, adart2);
        bool res = true;
        while (res && I1.cont() && I2.cont())
        {
          // We can remove this constraint which is not required for 
          // combinatorial map definition, but which is quite "normal"
          if ( I1==adart2 || I2==adart1 ) res=false;
          
          // Special case to consider beta0 and beta1
          if ( i>2 )
          {
            if ( I1->is_free(0)!=I2->is_free(1) )      res = false;
            else if ( I1->is_free(1)!=I2->is_free(0) ) res = false;
          }
          
          // General case
          for (unsigned int j=2;res && j<=Map::dimension; ++j)
          {
            if ( j+1!=i && j!=i && j!=i+1 && 
                 I1->is_free(j)!=I2->is_free(j) )
            { res = false; }
          }
          ++I1; ++I2;
        }
        if (I1.cont() != I2.cont()) 
          res = false;
      
        return res;
      }
    };
    
    /// Functor used to test if it is possible to 1-sew two darts.
    template <typename Map>
    struct is_sewable_functor<Map,1>{
      static bool run(const Map& amap,
                      typename Map::Dart_const_handle adart1,
                      typename Map::Dart_const_handle adart2)
      {
        CGAL_assertion(adart1!=NULL && adart2!=NULL);
      
        if ( !adart1->is_free(1) || !adart2->is_free(0) )
          return false;
      
        if ( adart1 == adart2 ) return true;
      
        CMap_dart_const_iterator_of_involution    <Map,1> I1(amap, adart1);
        CMap_dart_const_iterator_of_involution_inv<Map,1> I2(amap, adart2);
        bool res = true;
        while (res && I1.cont() && I2.cont())
        {
          // We can remove this constraint which is not required for 
          // combinatorial map definition, but which imposes quite "normal"
          // configurations
          if ( I1==adart2 || I2==adart1 ) res=false;
          
          for (unsigned int j=3;res && j<=Map::dimension; ++j)
          {
            if ( I1->is_free(j)!=I2->is_free(j) )
            {
              res = false;
            }
          }
          ++I1; ++I2;
        }
        if (I1.cont() != I2.cont()) 
          res = false;
      
        return res;
      }
    };

    /// Functor used to test if it is possible to 0-sew two darts.
    template <typename Map>
    struct is_sewable_functor<Map,0>{
      static bool run(const Map& amap,
                      typename Map::Dart_const_handle adart1,
                      typename Map::Dart_const_handle adart2)
      { return is_sewable_functor<Map,1>::run(amap,adart2,adart1); }
    };

    /// Functor used to i-topo_sew two darts, 2<=i<=dimension.
    template <typename Map,unsigned int i>
    struct topo_sew_functor{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_static_assertion(2<=i && i<=Map::dimension);
        CGAL_assertion( (is_sewable_functor<Map,i>::run(amap,adart1,adart2)) );

        CMap_dart_iterator_of_involution<Map,i>     I1(amap, adart1);
        CMap_dart_iterator_of_involution_inv<Map,i> I2(amap, adart2);
        while ( I1.cont() )        
        {
          amap.basic_link_beta(I1, I2, i);
          ++I1; ++I2;
        }
      }
    };

    /// Functor used to 1-topo_sew two darts.
    template <typename Map>
    struct topo_sew_functor<Map,1>{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_assertion( (is_sewable_functor<Map,1>::run(amap,adart1,adart2)) );

        int mark = amap.get_new_mark();
        std::vector<typename Map::Dart_handle> dartv;
        for (CMap_dart_iterator_basic_of_cell<Map,0> it(amap,adart1,mark); 
             it.cont(); ++it)
        {
          amap.mark(*it,mark);
          dartv.push_back(*it);
        }

        CMap_dart_iterator_of_involution<Map,1>     I1(amap, adart1);
        CMap_dart_iterator_of_involution_inv<Map,1> I2(amap, adart2);
        while ( I1.cont() )        
        {
          if ( amap.is_marked(*I1,mark) )
            amap.template basic_link_beta<1>(*I1, *I2);
          else
            amap.template basic_link_beta<0>(*I1, *I2);
          ++I1; ++I2;
        }

        for (typename std::vector<typename Map::Dart_handle>::iterator 
               it=dartv.begin(); it!=dartv.end(); ++it)
        { amap.unmark(*it,mark); }
        CGAL_assertion( amap.is_whole_map_unmarked(mark) );
        amap.free_mark(mark);
      }
    };

    /// Functor used to 0-topo_sew two darts.
    template <typename Map>
    struct topo_sew_functor<Map,0>{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      { topo_sew_functor<Map,1>::run(amap,adart2,adart1); }
    };

    /// Functor used to i-sew two darts, 2<=i<=dimension.
    template <typename Map,unsigned int i>
    struct sew_functor{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_static_assertion(2<=i && i<=Map::dimension);
        CGAL_assertion( (is_sewable_functor<Map,i>::run(amap,adart1,adart2)) );
      
        CMap_dart_iterator_of_involution<Map,i>     I1(amap, adart1);
        CMap_dart_iterator_of_involution_inv<Map,i> I2(amap, adart2);
        while ( I1.cont() )        
        {
          amap.group_all_attributes_except(I1,I2,i);
          ++I1; ++I2;
        }

        I1.rewind(); I2.rewind();      
        while ( I1.cont() )
        {
          amap.basic_link_beta(I1, I2, i);
          ++I1; ++I2;
        }
      }
    };

    /// Functor used to 0-sew two darts.
    template <typename Map>
    struct sew_functor<Map,0>{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {
        CGAL_assertion( (is_sewable_functor<Map,0>::run(amap,adart1,adart2)) );

        int mark = amap.get_new_mark();
        std::vector<typename Map::Dart_handle> dartv;
        for (CMap_dart_iterator_basic_of_cell<Map,0> it(amap,adart1,mark);
             it.cont(); ++it)
        {
          amap.mark(it,mark);
          dartv.push_back(it);
        }

        CMap_dart_iterator_of_involution<Map,1>     I1(amap, adart1);
        CMap_dart_iterator_of_involution_inv<Map,1> I2(amap, adart2);
        while ( I1.cont() )
        {
          typename Map::Dart_handle od1=I1->other_extremity();
          typename Map::Dart_handle od2=I2->other_extremity();
          if (od1!=NULL && od2!=NULL)
            amap.group_all_attributes_except(od1, od2, 1);
          ++I1; ++I2;	  
        }

        I1.rewind(); I2.rewind();      
        while ( I1.cont() )
        {
          if ( amap.is_marked(I1,mark) )
            amap.template basic_link_beta<0>(I1, I2);
          else
            amap.template basic_link_beta<1>(I1, I2);
          ++I1; ++I2;
        }

        for (typename std::vector<typename Map::Dart_handle>::iterator 
               it=dartv.begin(); it!=dartv.end(); ++it)
        { amap.unmark(*it,mark); }
        CGAL_assertion( amap.is_whole_map_unmarked(mark) );
        amap.free_mark(mark);
      }
    };

    /// Functor used to 1-sew two darts.
    template <typename Map>
    struct sew_functor<Map,1>{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {
        CGAL_assertion( (is_sewable_functor<Map,1>::run(amap,adart1,adart2)) );
        int mark = amap.get_new_mark();
        std::vector<typename Map::Dart_handle> dartv;
        for (CMap_dart_iterator_basic_of_cell<Map,0> it(amap,adart1,mark); 
             it.cont(); ++it)
        {
          amap.mark(it,mark);
          dartv.push_back(it);
        }

        CMap_dart_iterator_of_involution<Map,1>     I1(amap, adart1);
        CMap_dart_iterator_of_involution_inv<Map,1> I2(amap, adart2);
        while ( I1.cont() )
        {
          amap.group_all_attributes_except(I1,I2,1);
          ++I1; ++I2;
        }

        I1.rewind(); I2.rewind();      
        while ( I1.cont() )
        {
          if ( amap.is_marked(I1,mark) )
            amap.template basic_link_beta<1>(I1, I2);
          else
            amap.template basic_link_beta<0>(I1, I2);
          ++I1; ++I2;
        }

        for (typename std::vector<typename Map::Dart_handle>::iterator 
               it=dartv.begin(); it!=dartv.end(); ++it)
        { amap.unmark(*it,mark); }
        CGAL_assertion( amap.is_whole_map_unmarked(mark) );
        amap.free_mark(mark);
      }
    };

    /// Functor used to i-topo_unsew one dart, 2<=i<=dimension.
    template <typename Map,unsigned int i>
    struct topo_unsew_functor{
      static void run(Map& amap,typename Map::Dart_handle adart)
      {      
        CGAL_assertion( adart!=NULL && !adart->is_free(i) );
        CGAL_static_assertion(2<=i && i<=Map::dimension);

        CMap_dart_iterator_of_involution<Map,i> it(amap, adart);
        while ( it.cont() )
        {
          amap.unlink_beta(*it, i);
          ++it;
        }
      }
    };

    /// Functor used to 1-topo_unsew one dart.
    template <typename Map>
    struct topo_unsew_functor<Map,1>{
      static void run(Map& amap,typename Map::Dart_handle adart)
      {      
        CGAL_assertion( adart!=NULL && !adart->is_free(1) );

        int mark = amap.get_new_mark();
        std::vector<typename Map::Dart_handle> dartv;
        for (CMap_dart_iterator_basic_of_cell<Map,0> it(amap,adart,mark);
             it.cont(); ++it)
        {
          amap.mark(*it,mark);
          dartv.push_back(*it);
        }

        {
          CMap_dart_iterator_of_involution<Map,1> it(amap, adart);
          while ( it.cont() )        
          {
            if ( amap.is_marked(*it,mark) ) amap.unlink_beta<1>(*it);
            else amap.unlink_beta<0>(*it);
            ++it;
          }
        }

        for (typename std::vector<typename Map::Dart_handle>::iterator 
               it=dartv.begin(); it!=dartv.end(); ++it)
        { amap.unmark(*it,mark); }
        CGAL_assertion( amap.is_whole_map_unmarked(mark) );
        amap.free_mark(mark);
      }
    };

    /// Functor used to 1-topo_unsew one dart.
    template <typename Map>
    struct topo_unsew_functor<Map,0>{
      static void run(Map& amap,typename Map::Dart_handle adart)
      {   
        CGAL_assertion( adart!=NULL && !adart->is_free(0) );
        topo_unsew_functor<Map,1>::run(adart->beta(0));  
      }
    };

    /// Functor used to i-unsew one dart, 2<=i<=dimension.
    template <typename Map,unsigned int i>
    struct unsew_functor{
      static void run(Map& amap,typename Map::Dart_handle adart)
      {      
        CGAL_static_assertion(2<=i && i<=Map::dimension);
        CGAL_assertion( adart!=NULL && !adart->is_free(i) );

        std::stack<Couple_dart_and_dim<typename Map::Dart_handle> > todegroup;
      
        CMap_dart_iterator_of_involution<Map,i> it(amap, adart);
        while ( it.cont() )        
        {
          todegroup.push(Couple_dart_and_dim<typename Map::Dart_handle>
                         (it,it->beta(i),i));
          amap.unlink_beta(it, i);
          ++it;
        }

        while (!todegroup.empty() )
        {
          Couple_dart_and_dim<typename Map::Dart_handle> c=todegroup.top();
          todegroup.pop();
                amap.degroup_all_attributes_except(c.d1,c.d2,c.dim);
        }
      }
    };

    /// Functor used to 1-unsew one dart.
    template <typename Map>
    struct unsew_functor<Map,1>{
      static void run(Map& amap,typename Map::Dart_handle adart)
      {
        CGAL_assertion( adart!=NULL && !adart->is_free(1) );
        typename Map::Dart_handle d2 = NULL;

        int mark = amap.get_new_mark();
        std::vector<typename Map::Dart_handle> dartv;
        for (CMap_dart_iterator_basic_of_cell<Map,0> it(amap,adart,mark);
             it.cont(); ++it)
        {
          amap.mark(it,mark);
          dartv.push_back(it);
        }

        {
          CMap_dart_iterator_of_involution<Map,1> it(amap, adart);
          while ( it.cont() )        
          {
            if ( amap.is_marked(it,mark) ) 
            { d2 = it->beta(1); amap.template unlink_beta<1>(it); }
            else 
            { d2 = it->beta(0); amap.template unlink_beta<0>(it); }
            amap.degroup_all_attributes_except(it,d2,1);
            ++it;
          }
        }
      
        for (typename std::vector<typename Map::Dart_handle>::iterator 
               it=dartv.begin(); it!=dartv.end(); ++it)
        { amap.unmark(*it,mark); }
        CGAL_assertion( amap.is_whole_map_unmarked(mark) );
        amap.free_mark(mark);
      }
    };

    /// Functor used to 0-unsew one dart.
    template <typename Map>
    struct unsew_functor<Map,0>{
      static void run(Map& amap,typename Map::Dart_handle adart)
      {
        CGAL_assertion( adart!=NULL && !adart->is_free(0) );
        typename Map::Dart_handle d2 = NULL;

        int mark = amap.get_new_mark();
        std::vector<typename Map::Dart_handle> dartv;
        for (CMap_dart_iterator_basic_of_cell<Map,0> it(amap,adart,mark);
             it.cont(); ++it)
        {
          amap.mark(it,mark);
          dartv.push_back(it);
        }

        {
          CMap_dart_iterator_of_involution<Map,1> it(amap, adart);
          while ( it.cont() )        
          {
            if ( amap.is_marked(it,mark) ) 
            { d2 = it->beta(0); amap.template unlink_beta<0>(it); }
            else 
            { d2 = it->beta(1); amap.template unlink_beta<1>(it); }

            typename Map::Dart_handle od1=it->other_extremity();
            typename Map::Dart_handle od2=d2->other_extremity();
            if ( od1!=NULL && od2!=NULL )
              amap.degroup_all_attributes_except(od1,od2,1);

            ++it;
          }
        }

        for (typename std::vector<typename Map::Dart_handle>::iterator 
               it=dartv.begin(); it!=dartv.end(); ++it)
        { amap.unmark(*it,mark); }
        CGAL_assertion( amap.is_whole_map_unmarked(mark) );
        amap.free_mark(mark);
      }
    };

    /// Functor used to i-link two darts, 2<=i<=dimension.
    template <typename Map,unsigned int i>
    struct basic_link_beta_functor{
      static void run(Map&,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_assertion(adart1 != NULL && adart2 != NULL && adart1!=adart2);
        CGAL_static_assertion( i>=2 && i<=Map::dimension );
        adart1->basic_link_beta(adart2, i);
        adart2->basic_link_beta(adart1, i);
      }
    };

    /// Functor used to 0-link two darts.
    template <typename Map>
    struct basic_link_beta_functor<Map,0>{
      static void run(Map&,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_assertion(adart1 != NULL && adart2 != NULL );
        adart1->basic_link_beta(adart2, 0);
        adart2->basic_link_beta(adart1, 1);
      }
    };

    /// Functor used to 1-link two darts.
    template <typename Map>
    struct basic_link_beta_functor<Map,1>{
      static void run(Map& ,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_assertion(adart1 != NULL && adart2 != NULL);
        adart1->basic_link_beta(adart2, 1);
        adart2->basic_link_beta(adart1, 0);
      }
    };

    /// Functor used to i-unlink one dart.
    template <typename Map,unsigned int i>
    struct unlink_beta_functor{
      static void run(Map&,typename Map::Dart_handle adart)
      {      
        CGAL_assertion(adart != NULL && !adart->is_free(i));
        CGAL_static_assertion(2<=i && i<=Map::dimension);
        adart->beta(i)->unlink_beta(i);
        adart->unlink_beta(i);
      }
    };

    /// Functor used to 0-unlink one dart.
    template <typename Map>
    struct unlink_beta_functor<Map,0>{
      static void run(Map&,typename Map::Dart_handle adart)
      {      
        CGAL_assertion(adart != NULL && !adart->is_free(0));
        adart->beta(0)->unlink_beta(1);
        adart->unlink_beta(0);
      }
    };

    /// Functor used to 1-unlink one dart.
    template <typename Map>
    struct unlink_beta_functor<Map,1>{
      static void run(Map&,typename Map::Dart_handle adart)
      {
        CGAL_assertion(adart != NULL && !adart->is_free(1));
        adart->beta(1)->unlink_beta(0);
        adart->unlink_beta(1);      
      }
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

    /// Functor used to i-link two darts, 2<=i<=dimension.
    template <typename CMap,unsigned int i>
    struct link_beta_functor{
      static void run(CMap& amap,typename CMap::Dart_handle adart1,
                      typename CMap::Dart_handle adart2)
      {      
        CGAL_assertion(adart1 != NULL && adart2 != NULL && adart1!=adart2 );
        CGAL_static_assertion( 2<=i && i<=CMap::dimension );
        adart1->basic_link_beta(adart2, i);
        adart2->basic_link_beta(adart1, i);
        CMap::Helper::template Foreach_enabled_attributes
          <Group_attribute_functor_of_dart<CMap> >::run(&amap,adart1,adart2,i);
      }
    };
      
    /// Functor used to 0-link two darts.
    template <typename Map>
    struct link_beta_functor<Map,0>{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_assertion(adart1 != NULL && adart2 != NULL);
        adart1->basic_link_beta(adart2,0);
        adart2->basic_link_beta(adart1, 1);
        Map::Helper::template Foreach_enabled_attributes
          <Group_attribute_functor_of_dart<Map> >::run(&amap,adart1,adart2,0);
      }
    };
      
    /// Functor used to 1-link two darts.
    template <typename Map>
    struct link_beta_functor<Map,1>{
      static void run(Map& amap,typename Map::Dart_handle adart1,
                      typename Map::Dart_handle adart2)
      {      
        CGAL_assertion(adart1 != NULL && adart2 != NULL);
        adart1->basic_link_beta(adart2,1);
        adart2->basic_link_beta(adart1,0);
        Map::Helper::template Foreach_enabled_attributes
          <Group_attribute_functor_of_dart<Map> >::run(&amap,adart1,adart2,1);
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

  } // namespace internal
} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_FUNCTORS_H //
// EOF //
