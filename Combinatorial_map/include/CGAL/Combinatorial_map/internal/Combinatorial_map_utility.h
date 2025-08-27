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
#ifndef CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H
#define CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H 1

#include <CGAL/tuple.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Concurrent_compact_container.h>
#include <iostream>
#include <cstdint>
#include <type_traits>

#include <boost/function.hpp>
#include <boost/mpl/has_xxx.hpp>

/** Some utilities allowing to manage attributes. Indeed, as they as stores
 *  in tuples, we need to define functors with variadic templated arguments
 *  to deal with these attributes.
 *
 *  The class Combinatorial_map_helper<CMap> defines:
 *
 */
namespace CGAL
{
  namespace internal
  {
    // There is a problem on windows to handle tuple containing void.
    // To solve this, we transform such a tuple in tuple containing Void.
    template<typename T>
    struct Convert_void
    { typedef T type; };

    template<>
    struct Convert_void<void>
    { typedef CGAL::Void type; };

    // Get the type Dart_info defined as inner type of T.
    // If T::Dart_info is not defined or if T::Dart_info is void, defined
    // CGAL::Void as type.
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_dart_info,Dart_info,false)
    template<typename T, bool typedefined=Has_dart_info<T>::value >
    struct Get_dart_info
    { typedef CGAL::Void type; };
    template<typename T>
    struct Get_dart_info<T, true>
    { typedef typename Convert_void<typename T::Dart_info>::type type; };

    // Get the type Darts_with_id as inner type of T.
    // If T::Darts_with_id is not defined or if T::Darts_widh_id is Tag_false
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_darts_with_id,Darts_with_id,false)
    template<typename T, bool typedefined=Has_darts_with_id<T>::value >
    struct Get_darts_with_id
    { typedef CGAL::Tag_false type; };
    template<typename T>
    struct Get_darts_with_id<T, true>
    { typedef CGAL::Tag_true type; };

    // Get the type Attributes defined as inner type of T.
    // If T::Attributes is not defined, defined std::tuple<> as type.
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_attributes_tuple,Attributes,false)
    template<typename T, bool typedefined=Has_attributes_tuple<T>::value >
    struct Get_attributes_tuple
    { typedef std::tuple<> type; };
    template<typename T>
    struct Get_attributes_tuple<T, true>
    { typedef typename T::Attributes type; };

    // Get the Index_type
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_index_type,Index_type,false)
    template<typename T, typename I=typename T::Use_index,
             bool typedefined=Has_index_type<T>::value >
    struct Get_index_type
    { typedef std::uint32_t type; }; // By default use uint32_t for index type
    template<typename T>
    struct Get_index_type<T, CGAL::Tag_true, true>
    { typedef typename T::Index_type type; };

    // Get the Concurrent_tag
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_concurrent_tag,Use_concurrent_container,false)
    template<typename T, bool typedefined=Has_concurrent_tag<T>::value >
    struct Get_concurrent_tag
    { typedef CGAL::Tag_false type; };
    template<typename T>
    struct Get_concurrent_tag<T, true>
    { typedef CGAL::Tag_true type; };

    // Convert a tuple in a same tuple where each void type was replaced into
    // CGAL::Void.
    template<typename ... Items>
    struct Convert_tuple_with_void;
    template<typename ... Items>
    struct Convert_tuple_with_void<std::tuple<Items...> >
    {
      typedef std::tuple<typename Convert_void<Items>::type... > type;
    };

    // Length of a variadic template
    template<typename ... T>
    struct My_length;
    template<typename T1, typename ... T>
    struct My_length<std::tuple<T1, T...> >
    {
      static const int value = My_length<std::tuple<T...> >::value + 1;
    };
    template<>
    struct My_length<std::tuple<> >
    {
      static const int value = 0;
    };

    //count the number of time a given type is present in a tuple
    template<class Type,class Tuple>
    struct Number_of_type_in_tuple;
    template<class Type,typename ... Items>
    struct Number_of_type_in_tuple<Type,std::tuple<Type,Items...> >{
      static const int value=Number_of_type_in_tuple
        <Type,std::tuple<Items...> >::value+1;
    };
    template<class Type,class Other, typename ... Items>
    struct Number_of_type_in_tuple<Type,std::tuple<Other,Items...> >{
      static const int value=Number_of_type_in_tuple
        <Type,std::tuple<Items...> >::value;
    };
    template<class Type>
    struct Number_of_type_in_tuple<Type,std::tuple<> >{
      static const int value=0;
    };

    //count the number of different types from Type is present in a tuple
    template<class Type, class Tuple>
    struct Number_of_different_type_in_tuple;
    template<class Type, typename Other, typename ... Items>
    struct Number_of_different_type_in_tuple<Type,std::tuple
                                             <Other, Items...> >
    {
      static const int value=Number_of_different_type_in_tuple
        <Type,std::tuple<Items...> >::value+1;
    };
    template<class Type, typename ... Items>
    struct Number_of_different_type_in_tuple<Type, std::tuple
                                             <Type,Items...> >
    {
      static const int value=Number_of_different_type_in_tuple
        <Type,std::tuple<Items...> >::value;
    };
    template<class Type>
    struct Number_of_different_type_in_tuple<Type, std::tuple<> >
    {
      static const int value=0;
    };

    //count the number of time a given type have been found
    //within a tuple, until reaching position the k-th type of the tuple.
    //dim is the total size of the tuple
    template <class Type,int k,class T,
              int dim=CGAL::internal::My_length<T>::value-1>
    struct Nb_type_in_tuple_up_to_k;

    template <class Type,int dim,int k,class T1,class ... T>
    struct Nb_type_in_tuple_up_to_k<Type,k,std::tuple<T1,T...>,dim>
    {
      static const int pos= Nb_type_in_tuple_up_to_k
        <Type,k,std::tuple<T...>,dim>::pos - 1;

      static const int value =
        ( pos==k  ) ?  ( std::is_same<T1,Type>::value ? 0:-dim-1 )
        :  ( ( pos<k ) ? ( ( std::is_same<T1,Type>::value ? 1:0 )
                           + Nb_type_in_tuple_up_to_k
                           <Type,k,std::tuple
                           <T...>,dim >::value)
             :0
             );
    };

    template <class Type,int dim,int k,class T1>
    struct Nb_type_in_tuple_up_to_k<Type,k,std::tuple<T1>,dim >
    {
      static const int pos=dim;
      static const int value=(pos==k?
                              (std::is_same<T1,Type>::value?0:-dim-1) :
                              0);
    };

    //count the number of time a type different from Type have been found
    //within a tuple, until reaching position the k-th type of the tuple.
    //dim is the total size of the tuple
    template <class Type, int k,class T,
              int dim=CGAL::internal::My_length<T>::value-1>
    struct Nb_type_different_in_tuple_up_to_k;

    template <class Type,int dim,int k,class T1,class ... T>
    struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                              std::tuple<T1,T...>,dim>
    {
      static const int pos = Nb_type_different_in_tuple_up_to_k
        <Type,k,std::tuple<T...>,dim >::pos - 1;

      static const int value =
        ( pos==k  ) ?  ( std::is_same<T1,Type>::value ? -dim-1 : 0 )
        :  ( ( pos<k ) ? ( ( std::is_same<T1,Type>::value ? 0:1 )
                           + Nb_type_different_in_tuple_up_to_k
                           <Type,k,std::tuple<T...>,dim >::value)
             :0
             );
    };

    template <class Type,int dim,int k,class T1>
    struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                              std::tuple<T1>,dim >
    {
      static const int pos=dim;
      static const int value=(pos==k?
                              (std::is_same<T1,Type>::value?-dim-1:0) :
                              0);
    };

    //Convert a tuple of T... to a tuple of Functor<T>::type...
    template <template <class D> class Functor,class T>
    struct Tuple_converter;
    template <template <class D> class Functor,class ...T>
    struct Tuple_converter<Functor,std::tuple<T...> >{
      typedef std::tuple<typename Functor<T>::type... > type;
    };

    // To scan a given tuple, and keep only type different from Type
    // to build the tuple Attribute_type.
    template <class Type,class Res, class Tuple=std::tuple<> >
    struct Keep_type_different_of;

    template < class Type,class ... Res >
    struct Keep_type_different_of<Type,std::tuple<>,
                                  std::tuple<Res...> >
    {
      typedef std::tuple<Res...> type;
    };

    template < class Type,class ... T, class ... Res >
    struct Keep_type_different_of<Type,
                                  std::tuple<Type,T ...>,
                                  std::tuple<Res...> >
    {
      typedef typename Keep_type_different_of
      <Type,std::tuple<T ...>,std::tuple<Res...> >::type type;
    };

    template < class Type, class Other, class ... T, class ... Res >
    struct Keep_type_different_of<Type,std::tuple<Other,T...>,
                                  std::tuple<Res...> >
    {
      typedef typename Keep_type_different_of
      <Type, std::tuple<T...>,
       std::tuple<Res...,Other> >::type type;
    };

    //Helper class to statically call a functor
    // for (int i=n;i>=0;--i) Functor::run<i>(....)
    //Usage:
    //
    // struct Functor{
    //   template <int n>
    //   static void run(){std::cout << n << std::endl;}
    // };
    //
    // Foreach_static<Functor,5>::run();
    //
    template <class Functor,int n>
    struct Foreach_static{
      template <class  ... T>
      static void run(T& ... t){
        Functor:: template run<n>(t...);
        Foreach_static<Functor,n-1>::run(t...);
      }
    };

    template <class Functor>
    struct Foreach_static<Functor,0>{
      template <class  ... T>
      static void run(T& ... t)
      {
        Functor:: template run<0>( t... );
      }
    };

    //Helper function that is calling
    //Functor if TAG is different from Void
    template <class Functor,int n,class Type>
    struct Conditionnal_run{
      template <class ... T>
      static void run(T& ... t){
        Functor:: template run<n>(t...);
      }
    };

    template <class Functor,int n>
    struct Conditionnal_run<Functor,n,Void>
    {
      template <class ... T>
      static void run(T& ...){}
    };

    //Helper function that is calling
    //Functor if TAG is different from Void and n!=j
    template <class Functor,int n,int j,class Type>
    struct Conditionnal_run_except{
      template <class ... T>
      static void run(T& ... t){
        Functor:: template run<n>(t...);
      }
    };

    template <class Functor,int n,int j>
    struct Conditionnal_run_except<Functor,n,j,Void>
    {
      template <class ... T>
      static void run(T& ...){}
    };

    template <class Functor,int n,class Type>
    struct Conditionnal_run_except<Functor,n,n,Type>
    {
      template <class ... T>
      static void run(T& ...){}
    };

    template <class Functor,int n>
    struct Conditionnal_run_except<Functor,n,n,Void>
    {
      template <class ... T>
      static void run(T& ...){}
    };

    //Same as Foreach_static excepted that Functor
    //is called for case k only if the k-th type in the tuple
    //is different from Void. Note that to the converse of Foreach_static
    //Functor are called from n =0 to k
    template <class Functor,class T,int n=0, int startn=0>
    struct Foreach_static_restricted;

    template <class Functor,class Head, class ... Items,int n, int startn>
    struct Foreach_static_restricted<Functor,
                                     std::tuple<Head,Items...>,n, startn>
    {
      template <class  ... T>
      static void run(T& ... t){
        if(n>=startn)
        { Conditionnal_run<Functor,n,Head>::run(t...); }
        Foreach_static_restricted
          <Functor,std::tuple<Items...>, n+1, startn>::run(t...);
      }
    };

    template <class Functor,int n, int startn>
    struct Foreach_static_restricted<Functor,std::tuple<>,n, startn>{
      template <class  ... T>
      static void run(T& ... ){}
    };

    //Same as Foreach_static_restricted excepted that Functor
    //is called for case k only if the k-th type in the tuple
    //is different from Void and k!=j.
    template <class Functor,int j,class T,int n=0>
    struct Foreach_static_restricted_except;

    template <class Functor,int j,class Head, class ... Items,int n>
    struct Foreach_static_restricted_except<Functor, j,
        std::tuple<Head,Items...>,n>
    {
      template <class  ... T>
      static void run(T& ... t){
        Conditionnal_run_except<Functor,n,j,Head>::run(t...);
        Foreach_static_restricted_except
          <Functor,j,std::tuple<Items...>,n+1>::run(t...);
      }
    };

    template <class Functor,int j,int n>
    struct Foreach_static_restricted_except<Functor,j,std::tuple<>,n>
    {
      template <class  ... T>
      static void run(T& ... ){}
    };

    //Apply a functor to each element of a tuple
    template<class Functor,class Tuple,
             int pos=CGAL::internal::My_length<Tuple>::value-1>
    struct Apply_functor_to_each_tuple_element
    {
      static void run(Tuple& t){
        Functor() ( std::get<pos>(t) );
        Apply_functor_to_each_tuple_element<Functor,Tuple,pos-1>::run(t);
      }
    };

    template<class Functor,class Tuple>
    struct Apply_functor_to_each_tuple_element<Functor,Tuple,-1>
    {
      static void run(Tuple&){}
    };

    struct Clear_functor
    {
      template<class T>
      void operator()(T&t)
      { t.clear(); }
    };

    struct Clear_all
    {
      template<class Tuple>
      static void run(Tuple& t)
      {
        Apply_functor_to_each_tuple_element<Clear_functor,Tuple>::run(t);
      }
    };

    // Helper class, templated by a given combinatorial map.
    template <class CMap>
    struct Combinatorial_map_helper
    {
      // defines as type Compact_container<T>
      template <class T>
      struct Add_compact_container{
        typedef std::allocator_traits<typename CMap::Alloc> Allocator_traits;
        typedef typename Allocator_traits::template rebind_alloc<T> Attr_allocator;
        typedef typename CMap::template Container_for_attributes<T> type;
      };

      template<class Container, class UseIndex>
      struct Iterator_type
      {
        using type=typename Container::iterator;
        using const_type=typename Container::const_iterator;
      };

      template<class Container>
      struct Iterator_type<Container, CGAL::Tag_true>
      {
        using type=typename Container::Index;
        using const_type=typename Container::Index;
      };

      // defines as type Compact_container<T>::iterator
      template <class T>
      struct Add_compact_container_iterator{
        typedef std::allocator_traits<typename CMap::Alloc> Allocator_traits;
        typedef typename Allocator_traits::template rebind_alloc<T> Attr_allocator;

        // TODO? case when there is no Use_index typedef in CMap
        using type=typename Iterator_type<typename CMap::template
        Container_for_attributes<T>, typename CMap::Use_index>::type;
      };

      // defines as type Compact_container<T>::const_iterator
      template <class T>
      struct Add_compact_container_const_iterator{
        typedef std::allocator_traits<typename CMap::Alloc> Allocator_traits;
        typedef typename Allocator_traits::template rebind_alloc<T> Attr_allocator;

        using type=typename Iterator_type<typename CMap::template
        Container_for_attributes<T>, typename CMap::Use_index>::const_type;
      };

      // All the attributes (with CGAL::Void)
      typedef typename CGAL::internal::Convert_tuple_with_void
      <typename CMap::Attributes>::type Attributes;

      // defines as type Cell_attribute_binary_functor<T>
      template <class T>
      struct Define_cell_attribute_binary_functor{
        typedef typename boost::function<void(T&, T&)> type;
      };

      // Enabled attributes (without CGAL::Void)
      typedef typename CGAL::internal::Keep_type_different_of
      <CGAL::Void,Attributes>::type Enabled_attributes;

      // Number of all attributes
      /*  Does not compile on windows !!
          static const unsigned int number_of_attributes =
          CGAL::internal::My_length<Attributes>::value; */

      // Number of enabled attributes
      static const unsigned int nb_attribs =
          Number_of_different_type_in_tuple<Void,Enabled_attributes>::value;

      // Given a dimension of the cell, return the index of
      // corresponding attribute
      template <int d>
      struct Dimension_index
      { static const int value=
            Nb_type_different_in_tuple_up_to_k<Void,d,Attributes>::value; };

      // All these type contains as many entries than the number of
      // enabled attributes
      typedef typename Tuple_converter
      < Add_compact_container, Enabled_attributes >::type Attribute_containers;

      typedef typename Tuple_converter< Add_compact_container_iterator,
                                        Enabled_attributes >::type
      Attribute_iterators;
      typedef typename Tuple_converter< Add_compact_container_const_iterator,
                                        Enabled_attributes >::type
      Attribute_const_iterators;

      typedef Attribute_containers      Attribute_ranges;
      typedef Attribute_iterators       Attribute_descriptors;
      typedef Attribute_const_iterators Attribute_const_descriptors;

      typedef typename Tuple_converter< Define_cell_attribute_binary_functor,
                                        Enabled_attributes >::type
      Merge_functors;
      typedef typename Tuple_converter< Define_cell_attribute_binary_functor,
                                        Enabled_attributes >::type
      Split_functors;

      //Helper class allowing to retrieve the type of the
      // attribute of dimension d
      template<int d, int in_tuple=(d<CGAL::internal::My_length
                                    <Attributes>::value)>
      struct Attribute_type
      { typedef typename std::tuple_element<d,Attributes>::type type; };

      template<int d>
      struct Attribute_type<d,0>
      { typedef Void type; };

      // Helper class allowing to retrieve the d-cell-descriptor attribute
      template<int d, class Type=typename Attribute_type<d>::type,
               typename WithIndex=typename CMap::Use_index>
      struct Attribute_descriptor
      {
         typedef typename std::tuple_element
               <Dimension_index<d>::value,Attribute_descriptors>::type type;
      };

      template<int d>
      struct Attribute_descriptor<d, CGAL::Void, CGAL::Tag_false>
      { typedef CGAL::Void* type; };

      template<int d>
      struct Attribute_descriptor<d, CGAL::Void, CGAL::Tag_true>
      { typedef typename CMap::Dart_index type; };

      // Helper class allowing to retrieve the d-cell-const descriptor attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_const_descriptor
      {
        typedef typename std::tuple_element
          <Dimension_index<d>::value, Attribute_const_descriptors>::type type;
      };

      template<int d>
      struct Attribute_const_descriptor<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retrieve the d-cell-iterator attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_iterator
      {
        typedef typename std::tuple_element
           <Dimension_index<d>::value, Attribute_iterators>::type type;
      };

      template<int d>
      struct Attribute_iterator<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retrieve the d-cell-const descriptor attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_const_iterator
      {
        typedef typename std::tuple_element
           <Dimension_index<d>::value, Attribute_const_iterators>::type type;
      };

      template<int d>
      struct Attribute_const_iterator<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retrieve the d-cell-attribute range
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_range
      {
        typedef typename std::tuple_element
             <Dimension_index<d>::value, Attribute_ranges>::type type;
      };

      template<int d>
      struct Attribute_range<d, CGAL::Void>
      { typedef CGAL::Void type; };

      // Helper class allowing to retrieve the d-cell-attribute const range
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_const_range
      {
        typedef const typename std::tuple_element
             <Dimension_index<d>::value, Attribute_ranges >::type type;
      };

      template<int d>
      struct Attribute_const_range<d, CGAL::Void>
      { typedef CGAL::Void type; };

      // To iterate onto each enabled attributes, starting from startn-attributes (0 by default)
      template <class Functor, int startn=0>
      struct Foreach_enabled_attributes
      {
        template <class ...Ts>
        static void run(Ts& ... t)
        { Foreach_static_restricted<Functor, Attributes, 0, startn>::run(t...); }
      };
      // To iterate onto each enabled attributes, except j-attributes
      template <class Functor, unsigned int j>
      struct Foreach_enabled_attributes_except
      {
        template <class ...Ts>
        static void run(Ts& ... t)
        { Foreach_static_restricted_except<Functor, j, Attributes>::run(t...); }
      };
  };

  // Helper class to define container type depending on the Concurrent_tag
  template<typename Concurrent_tag, class T, class Alloc_>
  struct Container_type
  {
    typedef CGAL::Compact_container<T, Alloc_> type;
  };
  template<class T, class Alloc_>
  struct Container_type<CGAL::Tag_true, T, Alloc_>
  {
    typedef CGAL::Concurrent_compact_container<T, Alloc_> type;
  };

} //namespace internal

} //namespace CGAL

#endif //CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H
