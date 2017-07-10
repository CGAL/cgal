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
#ifndef CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H
#define CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H 1

#include <CGAL/tuple.h>
#include <CGAL/Compact_container.h>
#include <iostream>

#include <boost/type_traits/is_same.hpp>
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
    // If T::Attributes is not defined, defined CGAL::cpp11::tuple<> as type.
    BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_attributes_tuple,Attributes,false)
    template<typename T, bool typedefined=Has_attributes_tuple<T>::value >
    struct Get_attributes_tuple
    { typedef CGAL::cpp11::tuple<> type; };
    template<typename T>
    struct Get_attributes_tuple<T, true>
    { typedef typename T::Attributes type; };

#if ! defined(CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES) &&  \
    ! defined(CGAL_CFG_NO_CPP0X_TUPLE)
    // Convert a tuple in a same tuple where each void type was replaced into
    // CGAL::Void.
    template<typename ... Items>
    struct Convert_tuple_with_void;
    template<typename ... Items>
    struct Convert_tuple_with_void<CGAL::cpp11::tuple<Items...> >
    {
      typedef CGAL::cpp11::tuple<typename Convert_void<Items>::type... > type;
    };

    // Length of a variadic template
    template<typename ... T>
    struct My_length;
    template<typename T1, typename ... T>
    struct My_length<CGAL::cpp11::tuple<T1, T...> >
    {
      static const int value = My_length<CGAL::cpp11::tuple<T...> >::value + 1;
    };
    template<>
    struct My_length<CGAL::cpp11::tuple<> >
    {
      static const int value = 0;
    };

    //count the number of time a given type is present in a tuple
    template<class Type,class Tuple>
    struct Number_of_type_in_tuple;
    template<class Type,typename ... Items>
    struct Number_of_type_in_tuple<Type,CGAL::cpp11::tuple<Type,Items...> >{
      static const int value=Number_of_type_in_tuple
        <Type,CGAL::cpp11::tuple<Items...> >::value+1;
    };
    template<class Type,class Other, typename ... Items>
    struct Number_of_type_in_tuple<Type,CGAL::cpp11::tuple<Other,Items...> >{
      static const int value=Number_of_type_in_tuple
        <Type,CGAL::cpp11::tuple<Items...> >::value;
    };
    template<class Type>
    struct Number_of_type_in_tuple<Type,CGAL::cpp11::tuple<> >{
      static const int value=0;
    };

    //count the number of different types from Type is present in a tuple
    template<class Type, class Tuple>
    struct Number_of_different_type_in_tuple;
    template<class Type, typename Other, typename ... Items>
    struct Number_of_different_type_in_tuple<Type,CGAL::cpp11::tuple
                                             <Other, Items...> >
    {
      static const int value=Number_of_different_type_in_tuple
        <Type,CGAL::cpp11::tuple<Items...> >::value+1;
    };
    template<class Type, typename ... Items>
    struct Number_of_different_type_in_tuple<Type, CGAL::cpp11::tuple
                                             <Type,Items...> >
    {
      static const int value=Number_of_different_type_in_tuple
        <Type,CGAL::cpp11::tuple<Items...> >::value;
    };
    template<class Type>
    struct Number_of_different_type_in_tuple<Type, CGAL::cpp11::tuple<> >
    {
      static const int value=0;
    };

    //count the number of time a given type have been found
    //within a tuple, until reaching position the k'th type of the tuple.
    //dim is the total size of the tuple
    template <class Type,int k,class T,
              int dim=CGAL::internal::My_length<T>::value-1>
    struct Nb_type_in_tuple_up_to_k;

    template <class Type,int dim,int k,class T1,class ... T>
    struct Nb_type_in_tuple_up_to_k<Type,k,CGAL::cpp11::tuple<T1,T...>,dim>
    {
      static const int pos= Nb_type_in_tuple_up_to_k
        <Type,k,CGAL::cpp11::tuple<T...>,dim>::pos - 1;

      static const int value =
        ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? 0:-dim-1 )
        :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 1:0 )
                           + Nb_type_in_tuple_up_to_k
                           <Type,k,CGAL::cpp11::tuple
                           <T...>,dim >::value)
             :0
             );
    };

    template <class Type,int dim,int k,class T1>
    struct Nb_type_in_tuple_up_to_k<Type,k,CGAL::cpp11::tuple<T1>,dim >
    {
      static const int pos=dim;
      static const int value=(pos==k?
                              (boost::is_same<T1,Type>::value?0:-dim-1) :
                              0);
    };

    //count the number of time a type different from Type have been found
    //within a tuple, until reaching position the k'th type of the tuple.
    //dim is the total size of the tuple
    template <class Type, int k,class T,
              int dim=CGAL::internal::My_length<T>::value-1>
    struct Nb_type_different_in_tuple_up_to_k;

    template <class Type,int dim,int k,class T1,class ... T>
    struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                              CGAL::cpp11::tuple<T1,T...>,dim>
    {
      static const int pos = Nb_type_different_in_tuple_up_to_k
        <Type,k,CGAL::cpp11::tuple<T...>,dim >::pos - 1;

      static const int value =
        ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
        :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                           + Nb_type_different_in_tuple_up_to_k
                           <Type,k,CGAL::cpp11::tuple<T...>,dim >::value)
             :0
             );
    };

    template <class Type,int dim,int k,class T1>
    struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                              CGAL::cpp11::tuple<T1>,dim >
    {
      static const int pos=dim;
      static const int value=(pos==k?
                              (boost::is_same<T1,Type>::value?-dim-1:0) :
                              0);
    };

    //Convert a tuple of T... to a tuple of Functor<T>::type...
    template <template <class D> class Functor,class T>
    struct Tuple_converter;
    template <template <class D> class Functor,class ...T>
    struct Tuple_converter<Functor,CGAL::cpp11::tuple<T...> >{
      typedef CGAL::cpp11::tuple<typename Functor<T>::type... > type;
    };

    // To scan a given tuple, and keep only type different from Type
    // to build the tuple Attribute_type.
    template <class Type,class Res, class Tuple=CGAL::cpp11::tuple<> >
    struct Keep_type_different_of;

    template < class Type,class ... Res >
    struct Keep_type_different_of<Type,CGAL::cpp11::tuple<>,
                                  CGAL::cpp11::tuple<Res...> >
    {
      typedef CGAL::cpp11::tuple<Res...> type;
    };

    template < class Type,class ... T, class ... Res >
    struct Keep_type_different_of<Type,
                                  CGAL::cpp11::tuple<Type,T ...>,
                                  CGAL::cpp11::tuple<Res...> >
    {
      typedef typename Keep_type_different_of
      <Type,CGAL::cpp11::tuple<T ...>,CGAL::cpp11::tuple<Res...> >::type type;
    };

    template < class Type, class Other, class ... T, class ... Res >
    struct Keep_type_different_of<Type,CGAL::cpp11::tuple<Other,T...>,
                                  CGAL::cpp11::tuple<Res...> >
    {
      typedef typename Keep_type_different_of
      <Type, CGAL::cpp11::tuple<T...>,
       CGAL::cpp11::tuple<Res...,Other> >::type type;
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
    //is called for case k only if the k'th type in the tuple
    //is different from Void. Note that to the converse of Foreach_static
    //Functor are called from n =0 to k
    template <class Functor,class T,int n=0>
    struct Foreach_static_restricted;

    template <class Functor,class Head, class ... Items,int n>
    struct Foreach_static_restricted<Functor,
                                     CGAL::cpp11::tuple<Head,Items...>,n>
    {
      template <class  ... T>
      static void run(T& ... t){
        Conditionnal_run<Functor,n,Head>::run(t...);
        Foreach_static_restricted
          <Functor,CGAL::cpp11::tuple<Items...>,n+1>::run(t...);
      }
    };

    template <class Functor,int n>
    struct Foreach_static_restricted<Functor,CGAL::cpp11::tuple<>,n>{
      template <class  ... T>
      static void run(T& ... ){}
    };

    //Same as Foreach_static_restricted excepted that Functor
    //is called for case k only if the k'th type in the tuple
    //is different from Void and k!=j.
    template <class Functor,int j,class T,int n=0>
    struct Foreach_static_restricted_except;

    template <class Functor,int j,class Head, class ... Items,int n>
    struct Foreach_static_restricted_except<Functor, j,
        CGAL::cpp11::tuple<Head,Items...>,n>
    {
      template <class  ... T>
      static void run(T& ... t){
        Conditionnal_run_except<Functor,n,j,Head>::run(t...);
        Foreach_static_restricted_except
          <Functor,j,CGAL::cpp11::tuple<Items...>,n+1>::run(t...);
      }
    };

    template <class Functor,int j,int n>
    struct Foreach_static_restricted_except<Functor,j,CGAL::cpp11::tuple<>,n>
    {
      template <class  ... T>
      static void run(T& ... ){}
    };
#else
    // Definitions of structs are moved to another file.
#include <CGAL/internal/Combinatorial_map_utility_novariadic.h>
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

    //Apply a functor to each element of a tuple
    template<class Functor,class Tuple,
             int pos=CGAL::internal::My_length<Tuple>::value-1>
    struct Apply_functor_to_each_tuple_element
    {
      static void run(Tuple& t){
        Functor() ( CGAL::cpp11::get<pos>(t) );
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
        typedef typename CMap::Alloc::template rebind<T>::other Attr_allocator;
        typedef typename CMap::template Container_for_attributes<T> type;
      };

      // defines as type Compact_container<T>::iterator
      template <class T>
      struct Add_compact_container_iterator{
        typedef typename CMap::Alloc::template rebind<T>::other Attr_allocator;
        typedef typename CMap::template Container_for_attributes<T>::iterator
        iterator_type;

        // TODO case when there is no Use_index typedef in CMap
        typedef typename boost::mpl::if_
        < typename boost::is_same<typename CMap::Use_index,Tag_true>::type,
          typename CMap::Dart_handle, iterator_type >::type type;
      };

      // defines as type Compact_container<T>::const_iterator
      template <class T>
      struct Add_compact_container_const_iterator{
        typedef typename CMap::Alloc::template rebind<T>::other Attr_allocator;
        typedef typename CMap::template Container_for_attributes<T>::
        const_iterator iterator_type;

        typedef typename boost::mpl::if_
             < typename boost::is_same<typename CMap::Use_index,Tag_true>::type,
               typename CMap::Dart_handle, iterator_type >::type type;
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
      typedef Attribute_iterators       Attribute_handles;
      typedef Attribute_const_iterators Attribute_const_handles;

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
      { typedef typename CGAL::cpp11::tuple_element<d,Attributes>::type type; };

      template<int d>
      struct Attribute_type<d,0>
      { typedef Void type; };

      // Helper class allowing to retreive the d-cell-handle attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_handle
      {
         typedef typename CGAL::cpp11::tuple_element
               <Dimension_index<d>::value,Attribute_handles>::type type;
      };

      template<int d>
      struct Attribute_handle<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retreive the d-cell-const handle attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_const_handle
      {
        typedef typename CGAL::cpp11::tuple_element
          <Dimension_index<d>::value, Attribute_const_handles>::type type;
      };

      template<int d>
      struct Attribute_const_handle<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retreive the d-cell-iterator attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_iterator
      {
        typedef typename CGAL::cpp11::tuple_element
           <Dimension_index<d>::value, Attribute_iterators>::type type;
      };

      template<int d>
      struct Attribute_iterator<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retreive the d-cell-const handle attribute
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_const_iterator
      {
        typedef typename CGAL::cpp11::tuple_element
           <Dimension_index<d>::value, Attribute_const_iterators>::type type;
      };

      template<int d>
      struct Attribute_const_iterator<d, CGAL::Void>
      { typedef CGAL::Void* type; };

      // Helper class allowing to retreive the d-cell-attribute range
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_range
      {
        typedef typename CGAL::cpp11::tuple_element
             <Dimension_index<d>::value, Attribute_ranges>::type type;
      };

      template<int d>
      struct Attribute_range<d, CGAL::Void>
      { typedef CGAL::Void type; };

      // Helper class allowing to retreive the d-cell-attribute const range
      template<int d, class Type=typename Attribute_type<d>::type>
      struct Attribute_const_range
      {
        typedef const typename CGAL::cpp11::tuple_element
             <Dimension_index<d>::value, Attribute_ranges >::type type;
      };

      template<int d>
      struct Attribute_const_range<d, CGAL::Void>
      { typedef CGAL::Void type; };

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
      // To iterate onto each enabled attributes
      template <class Functor>
      struct Foreach_enabled_attributes
      {
        template <class ...Ts>
        static void run(Ts& ... t)
        { Foreach_static_restricted<Functor, Attributes>::run(t...); }
      };
      // To iterate onto each enabled attributes, except j-attributes
      template <class Functor, unsigned int j>
      struct Foreach_enabled_attributes_except
      {
        template <class ...Ts>
        static void run(Ts& ... t)
        { Foreach_static_restricted_except<Functor, j, Attributes>::run(t...); }
      };
#else
    // This one cannot be moved in Combinatorial_map_utility_novariadic.h
    // because this is an inner struct which uses inner type Attributes.
    template <class Functor>
    struct Foreach_enabled_attributes
    {
      static void run()
      {Foreach_static_restricted<Functor,Attributes >::run();}

      template <class T1>
      static void run(T1& t1)
      {Foreach_static_restricted<Functor,Attributes >::run(t1);}

      template <class T1,class T2>
      static void run(T1& t1,T2& t2)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2);}

      template <class T1,class T2,class T3>
      static void run(T1& t1,T2& t2,T3& t3)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3);}

      template <class T1,class T2,class T3,class T4>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4);}

      template <class T1,class T2,class T3,class T4,class T5>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5);}

      template <class T1,class T2,class T3,class T4,class T5,class T6>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,t6);}

      template <class T1,class T2,class T3,class T4,class T5,class T6,class T7>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6,T7& t7)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,
                                                           t6,t7);}

      template <class T1,class T2,class T3,class T4,class T5,class T6,
                class T7,class T8>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6,T7& t7,T8& t8)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,t6,
                                                           t7,t8);}

      template <class T1,class T2,class T3,class T4,class T5,class T6,
                class T7,class T8,class T9>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6,T7& t7,T8& t8,
                      T9& t9)
      {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,
                                                           t5,t6,t7,t8,t9);}
    };
    // This one cannot be moved in Combinatorial_map_utility_novariadic.h
    // because this is an inner struct which uses inner type Attributes.
    template <class Functor, unsigned int j>
    struct Foreach_enabled_attributes_except
    {
      static void run()
      {Foreach_static_restricted_except<Functor,j,Attributes>::run();}

      template <class T1>
      static void run(T1& t1)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1);}

      template <class T1,class T2>
      static void run(T1& t1,T2& t2)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2);}

      template <class T1,class T2,class T3>
      static void run(T1& t1,T2& t2,T3& t3)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3);}

      template <class T1,class T2,class T3,class T4>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3,t4);}

      template <class T1,class T2,class T3,class T4,class T5>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3,t4,t5);}

      template <class T1,class T2,class T3,class T4,class T5,class T6>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3,t4,t5,t6);}

      template <class T1,class T2,class T3,class T4,class T5,class T6,class T7>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6,T7& t7)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3,t4,t5,
                                                           t6,t7);}

      template <class T1,class T2,class T3,class T4,class T5,class T6,
                class T7,class T8>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6,T7& t7,T8& t8)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3,t4,t5,t6,
                                                           t7,t8);}

      template <class T1,class T2,class T3,class T4,class T5,class T6,
                class T7,class T8,class T9>
      static void run(T1& t1,T2& t2,T3& t3,T4& t4,
                      T5& t5,T6& t6,T7& t7,T8& t8,
                      T9& t9)
      {Foreach_static_restricted_except<Functor,j,Attributes>::run(t1,t2,t3,t4,
                                                           t5,t6,t7,t8,t9);}
    };
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  };

} //namespace internal

} //namespace CGAL

#endif //CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H
