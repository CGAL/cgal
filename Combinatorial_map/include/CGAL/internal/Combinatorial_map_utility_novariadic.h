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
#ifndef CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_NOVARIADIC_H
#define CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_NOVARIADIC_H 1

#ifdef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

// This file is included in Combinatorial_map_utility.h, in the namespace
// CGAL::internal
//------------------------------------------------------------------------------
template <class T1>
struct Convert_tuple_with_void;
    
template <>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<> >
{
  typedef CGAL::cpp0x::tuple<> type;
};
template <class T1>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type > type;
};
template <class T1, class T2>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type > type;
};
template <class T1, class T2, class T3>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type > type;
};
template <class T1, class T2, class T3, class T4>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type,
                             typename Convert_void<T4>::type > type;
};
template <class T1, class T2, class T3, class T4, class T5>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4, T5> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type,
                             typename Convert_void<T4>::type,
                             typename Convert_void<T5>::type > type;
};
template <class T1, class T2, class T3, class T4, class T5, class T6>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type,
                             typename Convert_void<T4>::type,
                             typename Convert_void<T5>::type,
                             typename Convert_void<T6>::type > type;
};
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type,
                             typename Convert_void<T4>::type,
                             typename Convert_void<T5>::type,
                             typename Convert_void<T6>::type,
                             typename Convert_void<T7>::type > type;
};
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7, class T8>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type,
                             typename Convert_void<T4>::type,
                             typename Convert_void<T5>::type,
                             typename Convert_void<T6>::type,
                             typename Convert_void<T7>::type,
                             typename Convert_void<T8>::type > type;
};
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7, class T8, class T9>
struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,
                                                  T7,T8,T9> >
{
  typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
                             typename Convert_void<T2>::type,
                             typename Convert_void<T3>::type,
                             typename Convert_void<T4>::type,
                             typename Convert_void<T5>::type,
                             typename Convert_void<T6>::type,
                             typename Convert_void<T7>::type,
                             typename Convert_void<T8>::type,
                             typename Convert_void<T9>::type > type;
};
//------------------------------------------------------------------------------
template <class T>
struct My_length;
    
template <>
struct My_length<CGAL::cpp0x::tuple<> >
{
  static const int value = 0;
};
template <class T1>
struct My_length<CGAL::cpp0x::tuple<T1> >
{
  static const int value = 1;
};
template <class T1, class T2>
struct My_length<CGAL::cpp0x::tuple<T1,T2> >
{
  static const int value = 2;
};
template <class T1, class T2, class T3>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3> >
{
  static const int value = 3;
};
template <class T1, class T2, class T3, class T4>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4> >
{
  static const int value = 4;
};
template <class T1, class T2, class T3, class T4, class T5>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4, T5> >
{
  static const int value = 5;
};
template <class T1, class T2, class T3, class T4, class T5, class T6>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >
{
  static const int value = 6;
};
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >
{
  static const int value = 7;
};
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7, class T8>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
{
  static const int value = 8;
};
template <class T1, class T2, class T3, class T4, class T5, class T6,
          class T7, class T8, class T9>
struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >
{
  static const int value = 9;
};
//------------------------------------------------------------------------------
template<class Type, class Tuple>
struct Number_of_different_type_in_tuple;
    
template <class Type>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<> >
{
  static const int value=0;
};
    
template <class Type, class T1>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1> >
{
  static const int value=boost::is_same<Type,T1>::value?0:1;
};
    
template <class Type, class T1,class T2>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2> >::value;
};
    
template <class Type, class T1,class T2,class T3>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                      <T2,T3> >::value;
};
    
template <class Type, class T1,class T2,class T3,class T4>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                         <T1,T2,T3,T4> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple
    <Type,CGAL::cpp0x::tuple<T2,T3,T4> >::value;
};
    
template <class Type, class T1,class T2,class T3,class T4,class T5>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,
                                                                 T4,T5> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple
    <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5> >::value;
};

template <class Type, class T1,class T2,class T3,class T4,class T5,class T6>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                         <T1,T2,T3,T4,T5,T6> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple
    <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6> >::value;
};
    
template <class Type, class T1,class T2,class T3,class T4,class T5,class T6,
          class T7>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                         <T1,T2,T3,T4,T5,T6,T7> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                      <T2,T3,T4,T5,T6,T7> >::value;
};
    
template <class Type, class T1,class T2,class T3,class T4,class T5,class T6,
          class T7,class T8>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                         <T1,T2,T3,T4,T5,T6,T7,T8> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple
    <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8> >::value;
};
    
template <class Type, class T1,class T2,class T3,class T4,class T5,class T6,
          class T7,class T8,class T9>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple
                                         <T1,T2,T3,T4,T5,T6,T7,T8,T9> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) +
    Number_of_different_type_in_tuple
    <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9> >::value;
};
//------------------------------------------------------------------------------
template <class Type, int k,class T,
          int dim=CGAL::internal::My_length<T>::value-1>
struct Nb_type_different_in_tuple_up_to_k;
    
template <class Type,int dim,int k,class T1>
struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                          CGAL::cpp0x::tuple<T1>,dim >
{
  static const int pos=dim;
  static const int value= pos==k?
    boost::is_same<T1,Type>::value?-dim-1:0 : 0;
};
    
template <class Type,int dim,int k,class T1,class T2>
struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                          CGAL::cpp0x::tuple<T1,T2>,dim >
{
  static const int pos= Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2>,dim >::pos - 1;
      
  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple<T2>,dim >::value)
         :0
         );
};

template <class Type,int dim,int k,class T1,class T2,class T3>
struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                          CGAL::cpp0x::tuple<T1,T2,T3>,dim >
{
  static const int pos = Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k<Type,k,
                       CGAL::cpp0x::tuple<T2,T3>,dim >::value)
         :0
         );
};


template <class Type,int dim,int k,class T1,class T2,class T3,class T4>
struct Nb_type_different_in_tuple_up_to_k<Type,k,
                                          CGAL::cpp0x::tuple<T1,T2,T3,T4>,dim >
{
  static const int pos= Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3,T4>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple<T2,T3,T4>,dim >::value)
         :0
         );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,
          class T5>
struct Nb_type_different_in_tuple_up_to_k
<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5>,dim >
{
  static const int pos=Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5>,dim >::value)
         :0
         );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,
          class T5,class T6>
struct Nb_type_different_in_tuple_up_to_k
<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6>,dim >
{
  static const int pos=Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple
                       <T2,T3,T4,T5,T6>,dim >::value)
         :0
         );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7>
struct Nb_type_different_in_tuple_up_to_k
<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7>,dim >{
  static const int pos=Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7>,
                       dim >::value)
         :0
         );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8>
struct Nb_type_different_in_tuple_up_to_k
<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8>,dim >
{
  static const int pos=Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8>,
                       dim >::value)
         :0
         );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9>
struct Nb_type_different_in_tuple_up_to_k
<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9>,dim >
{
  static const int pos=Nb_type_different_in_tuple_up_to_k
    <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
    :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
                       + Nb_type_different_in_tuple_up_to_k
                       <Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9>,
                       dim >::value)
         :0
         );
};
//------------------------------------------------------------------------------
//Convert a tuple of T... to a tuple of Functor<T>::type...
template <template <class D> class Functor,class T>
struct Tuple_converter;
    
template <template <class D> class Functor>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<> >{
  typedef CGAL::cpp0x::tuple<> type;
};

template <template <class D> class Functor,class T1>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type > type;
};

template <template <class D> class Functor,class T1,class T2>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type > type;
};

template <template <class D> class Functor,class T1,class T2,class T3>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type ,
                             typename Functor<T3>::type > type;
};

template <template <class D> class Functor,class T1,class T2,class T3,
          class T4>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4> >
{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type,
                             typename Functor<T3>::type,
                             typename Functor<T4>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,
          class T4,class T5>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type,
                             typename Functor<T3>::type,
                             typename Functor<T4>::type,
                             typename Functor<T5>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,
          class T4,class T5,class T6>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type,
                             typename Functor<T3>::type,
                             typename Functor<T4>::type,
                             typename Functor<T5>::type,
                             typename Functor<T6>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,
          class T4,class T5,class T6,class T7>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type,
                             typename Functor<T3>::type,
                             typename Functor<T4>::type,
                             typename Functor<T5>::type,
                             typename Functor<T6>::type,
                             typename Functor<T7>::type> type;
};
    
template <template <class D> class Functor,class T1,class T2,class T3,
          class T4,class T5,class T6,class T7,class T8>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type,
                             typename Functor<T3>::type,
                             typename Functor<T4>::type,
                             typename Functor<T5>::type,
                             typename Functor<T6>::type,
                             typename Functor<T7>::type,
                             typename Functor<T8>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,
          class T4,class T5,class T6,class T7,class T8,class T9>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,
                                                  T8,T9> >
{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type,
                             typename Functor<T2>::type,
                             typename Functor<T3>::type,
                             typename Functor<T4>::type,
                             typename Functor<T5>::type,
                             typename Functor<T6>::type,
                             typename Functor<T7>::type,
                             typename Functor<T8>::type,
                             typename Functor<T9>::type> type;
};
//------------------------------------------------------------------------------
template <class Type,class Tuple>
struct Append_to_tuple;
  
template <class Type>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<> >
{ typedef CGAL::cpp0x::tuple<Type> type; };

template <class Type,class T1>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1> >
{ typedef CGAL::cpp0x::tuple<T1,Type> type; };

template <class Type,class T1,class T2>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2> >
{ typedef CGAL::cpp0x::tuple<T1,T2,Type> type; };

template <class Type,class T1,class T2,class T3>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,Type> type; };

template <class Type,class T1,class T2,class T3,class T4>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6,
          class T7>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6,
          class T7,class T8>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6,
          class T7,class T8,class T9>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >
{ typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,Type> type; };
//------------------------------------------------------------------------------
template <class Type,class Tuple,class Result=CGAL::cpp0x::tuple<> >
struct Keep_type_different_of;

template <class Type,class Result>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<>,Result> {
  typedef Result type;
};  
  
template <class Type,class T1,class Result>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1>,Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type type;
};

template <class Type,class Result,class T1,class T2>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2>,Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2>,New_result>::type type;
};
    
template <class Type,class Result,class T1,class T2,class T3>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3>,Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4>,Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3,T4>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,
          class T5>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5>,
                              Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,
          class T5,class T6>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6>,
                              Result>
{
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6>,New_result>::type type;
};
    
template <class Type,class Result,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7>,
                              Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,
                                                      T7,T8>,Result> {
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,
          class T5,class T6,class T7,class T8,class T9>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,
                                                      T8,T9>,Result>{
  typedef typename boost::mpl::if_
  < typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple
    <T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of
  <Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9>,New_result>::type type;
};
//------------------------------------------------------------------------------
template <class Functor,int n>
struct Foreach_static
{
  static void run() 
  {
    Functor:: template run<n>();
    Foreach_static<Functor,n-1>::run();
  }

  template<class T1>
  static void run(const T1& t1) 
  {
    Functor:: template run<n>(t1);
    Foreach_static<Functor,n-1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2) 
  {
    Functor:: template run<n>(t1,t2);
    Foreach_static<Functor,n-1>::run(t1,t2);
  }
  
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3)
  {
    Functor:: template run<n>(t1,t2,t3);
    Foreach_static<Functor,n-1>::run(t1,t2,t3);
  }
  
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4)
  {
    Functor:: template run<n>(t1,t2,t3,t4);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4);
  }
  
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5);
  }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6,t7);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor>
struct Foreach_static<Functor,0>
{
  static void run() 
  {
    Functor:: template run<0>();
  }

  template<class T1>
  static void run(const T1& t1) 
  {
    Functor:: template run<0>(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2) 
  {
    Functor:: template run<0>(t1,t2);
  }
  
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3)
  {
    Functor:: template run<0>(t1,t2,t3);
  }
  
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4)
  {
    Functor:: template run<0>(t1,t2,t3,t4);
  }
  
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5);
  }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6,t7);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6,t7,t8);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,
           class T7,class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }  
};
//------------------------------------------------------------------------------
template <class Functor,int n,class Type>
struct Conditionnal_run{
  static void run() 
  {
    Functor:: template run<n>();
  }

  template<class T1>
  static void run(const T1& t1) 
  {
    Functor:: template run<n>(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2) 
  {
    Functor:: template run<n>(t1,t2);
  }
  
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3)
  {
    Functor:: template run<n>(t1,t2,t3);
  }
  
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4)
  {
    Functor:: template run<n>(t1,t2,t3,t4);
  }
  
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5);
  }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }  
};
//------------------------------------------------------------------------------
template <class Functor,int n>
struct Conditionnal_run<Functor,n,Void>
{
  static void run(){}

  template<class T1>
  static void run(const T1&){}

  template<class T1,class T2>
  static void run(const T1&,const T2&){}
  
  template<class T1,class T2,class T3>
  static void run(const T1&,const T2&,const T3&){}
  
  template<class T1,class T2,class T3,class T4>
  static void run(const T1&,const T2&,const T3&,const T4&){}
  
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&){}
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,
                  const T6&){}
      
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,
                  const T6&,const T7&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,
                  const T6&,const T7&,const T8&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,
                  const T6&,const T7&,const T8&,const T9&){}
};
//------------------------------------------------------------------------------
template <class Functor,class T,int n=0>
struct Foreach_static_restricted;

template <class Functor,class TT1,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::
      run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::
      run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::
      run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,class TT4,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
                                 <TT1,TT2,TT3,TT4>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7);
  }  
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,class TT4,
          class TT5,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
                                 <TT1,TT2,TT3,TT4,TT5>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,
          class TT6,int n>
struct Foreach_static_restricted
<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5,TT6>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,
          class TT6,class TT7,int n>
struct Foreach_static_restricted
<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,
      n+1>::run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,
      CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7>,
      n+1>::run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7>,
      n+1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,
          class TT6,class TT7,class TT8,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
                                 <TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run();
  }
  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1);
  }
  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1,t2);
  }
  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1,t2,t3);
  }
  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::
      run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::
      run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::
      run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,
          class TT6,class TT7,class TT8,class TT9,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
                                 <TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3,t4);
  }
  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3,t4,t5);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3,t4,t5,t6);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
      <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8);
  }
  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,
           class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,
                  const T5& t5,const T6& t6,const T7& t7,const T8& t8,
                  const T9& t9){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple
                              <TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::
      run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};
//------------------------------------------------------------------------------
template <class Functor,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n>{
  static void run(){}

  template<class T1>
  static void run(const T1&){}

  template<class T1,class T2>
  static void run(const T1&,const T2&){}
  
  template<class T1,class T2,class T3>
  static void run(const T1&,const T2&,const T3&){}
  
  template<class T1,class T2,class T3,class T4>
  static void run(const T1&,const T2&,const T3&,const T4&){}

  template<class T1,class T2,class T3,class T4,class T5>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6, class T7>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&,
                  const T7&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6, class T7,
           class T8>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&,
                  const T7&,const T8&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6, class T7,
           class T8, class T9>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&,
                  const T7&,const T8&,const T9&){}
};

template<typename Dart_handle>
struct Beta_functor
{
  static Dart_handle run(Dart_handle ADart, int B)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3, int B4)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3)->beta(B4);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3, int B4,
                         int B5)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3)->beta(B4)->beta(B5);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3, int B4,
                         int B5, int B6)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3)->beta(B4)->beta(B5)->
      beta(B6);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3, int B4,
                         int B5, int B6, int B7)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3)->beta(B4)->beta(B5)->
      beta(B6)->beta(B7);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3, int B4,
                         int B5, int B6, int B7, int B8)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3)->beta(B4)->beta(B5)->
      beta(B6)->beta(B7)->beta(B8);
  }

  static Dart_handle run(Dart_handle ADart, int B1, int B2, int B3, int B4, int B5, int B6,
                         int B7, int B8, int B9)
  {
    CGAL_assertion( ADart!=NULL );
    return ADart->beta(B1)->beta(B2)->beta(B3)->beta(B4)->beta(B5)->
      beta(B6)->beta(B7)->beta(B8)->beta(B9);
  }
};


#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

#endif //CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_NOVARIADIC_H

