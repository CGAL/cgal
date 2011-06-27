// Copyright (c) 2010-2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#include <boost/static_assert.hpp>

#ifdef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>
#endif


namespace CGAL 
{
  // struct Disabled {}; // If we want to use CGAL::Disabled for disabled attributes
  // typedef CGAL::Void Disabled; 
  // typedef void Disabled; // If we wand to use void, does not compile on windows

  namespace internal
  {    
	  // There is a problem on windows to handle tuple containing void. To solve this, we
	  // transform such a tuple in tuple containing Disabled.
	template<typename T>
	struct Convert_void
	{ typedef T type; };

	template<>
	struct Convert_void<void>
	{ typedef CGAL::Void type; };

	#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
	template<typename ... Items>
	struct Convert_tuple_with_void;
    
	template<typename ... Items>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<Items...> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<Items>::type... > type;
	};
    
        template<typename ... T>
        struct My_length;    
        template<typename T1, typename ... T>
	struct My_length<CGAL::cpp0x::tuple<T1, T...> >
	{
	  static const int value = My_length<CGAL::cpp0x::tuple<T...> >::value + 1;
	};
    
        template<>
	struct My_length<CGAL::cpp0x::tuple<> >
	{
	  static const int value = 0;
	};
    

#else //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
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
			     typename Convert_void<T2>::type> type;
	};
	template <class T1, class T2, class T3>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type> type;
	};
	template <class T1, class T2, class T3, class T4>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type,
			     typename Convert_void<T4>::type> type;
	};
	template <class T1, class T2, class T3, class T4, class T5>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4, T5> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type,
			     typename Convert_void<T4>::type,
			     typename Convert_void<T5>::type> type;
	};
	template <class T1, class T2, class T3, class T4, class T5, class T6>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type,
			     typename Convert_void<T4>::type,
			     typename Convert_void<T5>::type,
			     typename Convert_void<T6>::type> type;
	};
	template <class T1, class T2, class T3, class T4, class T5, class T6, class T7>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type,
			     typename Convert_void<T4>::type,
			     typename Convert_void<T5>::type,
			     typename Convert_void<T6>::type,
			     typename Convert_void<T7>::type> type;
	};
	template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
			  class T8>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type,
			     typename Convert_void<T4>::type,
			     typename Convert_void<T5>::type,
			     typename Convert_void<T6>::type,
			     typename Convert_void<T7>::type,
			     typename Convert_void<T8>::type> type;
	};
	template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
			  class T8, class T9>
	struct Convert_tuple_with_void<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >
	{
	typedef CGAL::cpp0x::tuple<typename Convert_void<T1>::type,
			     typename Convert_void<T2>::type,
			     typename Convert_void<T3>::type,
			     typename Convert_void<T4>::type,
			     typename Convert_void<T5>::type,
			     typename Convert_void<T6>::type,
			     typename Convert_void<T7>::type,
			     typename Convert_void<T8>::type,
			     typename Convert_void<T9>::type> type;
	};
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
	template <class T1, class T2, class T3, class T4, class T5, class T6, class T7>
	struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >
	{
	  static const int value = 7;
	};
	template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
			  class T8>
	struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
	{
	  static const int value = 8;
	};
	template <class T1, class T2, class T3, class T4, class T5, class T6, class T7,
			  class T8, class T9>
	struct My_length<CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >
	{
	  static const int value = 9;
	};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
    
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
//count the number of time a given type is present in a tuple
template<class Type,class Tuple>
struct Number_of_type_in_tuple;
  
template<class Type,typename ... Items>
struct Number_of_type_in_tuple<Type,CGAL::cpp0x::tuple<Type,Items...> >{
  static const int value=Number_of_type_in_tuple<Type,CGAL::cpp0x::tuple<Items...> >::value+1;
};

template<class Type,class Other, typename ... Items>
struct Number_of_type_in_tuple<Type,CGAL::cpp0x::tuple<Other,Items...> >{
  static const int value=Number_of_type_in_tuple<Type,CGAL::cpp0x::tuple<Items...> >::value;
};

template<class Type>
struct Number_of_type_in_tuple<Type,CGAL::cpp0x::tuple<> >{
  static const int value=0;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//count the number of different types from Type is present in a tuple
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template<class Type, class Tuple>
struct Number_of_different_type_in_tuple;
  
template<class Type, typename Other, typename ... Items>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<Other, Items...> >
{
  static const int value=Number_of_different_type_in_tuple<Type,
         CGAL::cpp0x::tuple<Items...> >::value+1;
};

template<class Type, typename ... Items>
struct Number_of_different_type_in_tuple<Type, CGAL::cpp0x::tuple<Type,Items...> >
{
  static const int value=Number_of_different_type_in_tuple<Type,
      CGAL::cpp0x::tuple<Items...> >::value;
};

template<class Type>
struct Number_of_different_type_in_tuple<Type, CGAL::cpp0x::tuple<> >
{
  static const int value=0;
};
#else
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
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2> >::value;
};

template <class Type, class T1,class T2,class T3>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3> >::value;
};

template <class Type, class T1,class T2,class T3,class T4>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3,T4> >::value;
};

template <class Type, class T1,class T2,class T3,class T4,class T5>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5> >::value;
};

template <class Type, class T1,class T2,class T3,class T4,class T5,class T6>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6> >::value;
};

template <class Type, class T1,class T2,class T3,class T4,class T5,class T6,class T7>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7> >::value;
};

template <class Type, class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8> >::value;
};

template <class Type, class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
struct Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >
{
  static const int value= (boost::is_same<Type,T1>::value?0:1) + Number_of_different_type_in_tuple<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9> >::value;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES


#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
//count the number of time a given type have been found
//within a tuple, until reaching position the k'th type of the tuple.
//dim is the total size of the tuple
template <class Type,int k,class T,int dim=CGAL::internal::My_length<T>::value-1>
struct Nb_type_in_tuple_up_to_k;

template <class Type,int dim,int k,class T1,class ... T>
struct Nb_type_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T...>,dim>{
  static const int pos= Nb_type_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T...>,dim>::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? 0:-dim-1 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 1:0 )
                                     + Nb_type_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T...>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1>
struct Nb_type_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1>,dim >{
  static const int pos=dim;
  static const int value= pos==k? boost::is_same<T1,Type>::value?0:-dim-1 : 0;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//count the number of time a type different from Type have been found
//within a tuple, until reaching position the k'th type of the tuple.
//dim is the total size of the tuple
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class Type, int k,class T,int dim=CGAL::internal::My_length<T>::value-1>
struct Nb_type_different_in_tuple_up_to_k;

template <class Type,int dim,int k,class T1,class ... T>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T...>,dim>{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T...>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T...>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1>,dim >{
  static const int pos=dim;
  static const int value= pos==k? boost::is_same<T1,Type>::value?-dim-1:0 : 0;
};
#else
template <class Type, int k,class T,int dim=CGAL::internal::My_length<T>::value-1>
struct Nb_type_different_in_tuple_up_to_k;

template <class Type,int dim,int k,class T1>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1>,dim >{
  static const int pos=dim;
  static const int value= pos==k? boost::is_same<T1,Type>::value?-dim-1:0 : 0;
};

template <class Type,int dim,int k,class T1,class T2>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1,class T2,class T3>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3>,dim >::value)
                               :0
                   );
};


template <class Type,int dim,int k,class T1,class T2,class T3,class T4>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3,T4>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3,T4>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,class T5>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3,T4,T5>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,class T5,class T6>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3,T4,T5,T6>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8>,dim >::value)
                               :0
                   );
};

template <class Type,int dim,int k,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
struct Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9>,dim >{
  static const int pos= 
    Nb_type_different_in_tuple_up_to_k<Type,k,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9>,dim >::pos - 1;

  static const int value =
    ( pos==k  ) ?  ( boost::is_same<T1,Type>::value ? -dim-1 : 0 )
                :  ( ( pos<k ) ? ( ( boost::is_same<T1,Type>::value ? 0:1 )
		   + Nb_type_different_in_tuple_up_to_k<Type,k,
                                   CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9>,dim >::value)
                               :0
                   );
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//Convert a tuple of T... to a tuple of Functor<T>::type...
template <template <class D> class Functor,class T>
struct Tuple_converter;
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <template <class D> class Functor,class ...T>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T...> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T>::type... > type;
};
#else
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
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type > type;
};

template <template <class D> class Functor,class T1,class T2,class T3>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type , typename Functor<T3>::type > type;
};

template <template <class D> class Functor,class T1,class T2,class T3,class T4>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type, typename Functor<T3>::type, typename Functor<T4>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,class T4,class T5>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type, typename Functor<T3>::type, typename Functor<T4>::type, typename Functor<T5>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,class T4,class T5,class T6>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type, typename Functor<T3>::type, typename Functor<T4>::type, typename Functor<T5>::type, typename Functor<T6>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type, typename Functor<T3>::type, typename Functor<T4>::type, typename Functor<T5>::type, typename Functor<T6>::type, typename Functor<T7>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type, typename Functor<T3>::type, typename Functor<T4>::type, typename Functor<T5>::type, typename Functor<T6>::type, typename Functor<T7>::type, typename Functor<T8>::type> type;
};

template <template <class D> class Functor,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
struct Tuple_converter<Functor,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >{
  typedef CGAL::cpp0x::tuple<typename Functor<T1>::type, typename Functor<T2>::type, typename Functor<T3>::type, typename Functor<T4>::type, typename Functor<T5>::type, typename Functor<T6>::type, typename Functor<T7>::type, typename Functor<T8>::type, typename Functor<T9>::type> type;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//Apply a functor to each element of a tuple
template<class Functor,class Tuple,int pos=CGAL::internal::My_length<Tuple>::value-1>
struct Apply_functor_to_each_tuple_element
{
  static void run(Tuple& t){
    Functor() ( CGAL::cpp0x::get<pos>(t) );
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

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
//if the number of elements in T is lower than size,
//append void
template <int size,class Tuple,int remaining = size - CGAL::internal::My_length<Tuple>::value>
struct Fill_tag_false;

template <int size, class ... T,int remaining>
struct Fill_tag_false<size,CGAL::cpp0x::tuple<T...>,remaining>
{
  typedef typename Fill_tag_false<size,CGAL::cpp0x::tuple<T...,Tag_false>,remaining - 1>::type type;
};

template <int size,class ... T>
struct Fill_tag_false<size,CGAL::cpp0x::tuple<T...>,0 >
{
  typedef CGAL::cpp0x::tuple<T...> type;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//if the number of elements in T is lower than size,
//append Type
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class Type,int size, class Tuple, 
	  int remaining = size - CGAL::internal::My_length<Tuple>::value>
struct Fill_type;

template <class Type, int size, class ... T, int remaining>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T...>,remaining>
{
  typedef typename Fill_type<Type,size,
			     CGAL::cpp0x::tuple<T...,Void>,remaining - 1>::type type;
};

template < class Type,int size, class ... T >
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T...>,0 >
{
  typedef CGAL::cpp0x::tuple<T...> type;
};

template < int size, class Tuple >
struct Fill_disabled
{
  typedef typename Fill_type<Void,size,Tuple>::type type;
};
#else
template <class Type,int size, class Tuple>
struct Fill_type{
  typedef void type;
   BOOST_STATIC_ASSERT(size < 10);
};

template <class Type,int size>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<> >{
  typedef typename boost::mpl::if_c<size==0,CGAL::cpp0x::tuple<>,typename Fill_type<Type,size,CGAL::cpp0x::tuple<Type> >::type >::type type;
};

template <class Type,int size, class T1>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1> >{
  typedef typename boost::mpl::if_c<size==1,CGAL::cpp0x::tuple<T1>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2> >{
  typedef typename boost::mpl::if_c<size==2,CGAL::cpp0x::tuple<T1,T2>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3> >{
  typedef typename boost::mpl::if_c<size==3,CGAL::cpp0x::tuple<T1,T2,T3>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3,class T4>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4> >{
  typedef typename boost::mpl::if_c<size==4,CGAL::cpp0x::tuple<T1,T2,T3,T4>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3,class T4,class T5>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5> >{
  typedef typename boost::mpl::if_c<size==5,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3,class T4,class T5,class T6>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >{
  typedef typename boost::mpl::if_c<size==6,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3,class T4,class T5,class T6,class T7>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >{
  typedef typename boost::mpl::if_c<size==7,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >{
  typedef typename boost::mpl::if_c<size==8,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,Type> >::type  >::type type;
};

template <class Type,int size, class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
struct Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >{
  typedef typename boost::mpl::if_c<size==9,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9>, typename Fill_type<Type,size,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,Type> >::type  >::type type;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES



// To scan a given tuple, and keep only type different from Type
// to build the tuple Attribute_type.
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class Type,class Res, class Tuple=CGAL::cpp0x::tuple<> >
struct Keep_type_different_of;

template < class Type,class ... Res >
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<>, CGAL::cpp0x::tuple<Res...> >
{
  typedef CGAL::cpp0x::tuple<Res...> type;
};

template < class Type,class ... T, class ... Res >
struct Keep_type_different_of<Type,
		    CGAL::cpp0x::tuple<Type,T ...>, CGAL::cpp0x::tuple<Res...> >
{  
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T ...>,
    CGAL::cpp0x::tuple<Res...> >::type type;
};

template < class Type, class Other, class ... T, class ... Res >
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<Other,T...>, 
				       CGAL::cpp0x::tuple<Res...> >
{
  typedef typename  Keep_type_different_of<Type, CGAL::cpp0x::tuple<T...>,
    CGAL::cpp0x::tuple<Res...,Other> >::type type;
};
#else
template <class Type,class Tuple>
struct Append_to_tuple;
  
template <class Type>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<> >{  typedef CGAL::cpp0x::tuple<Type> type; };

template <class Type,class T1>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1> >{  typedef CGAL::cpp0x::tuple<T1,Type> type; };

template <class Type,class T1,class T2>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2> >{  typedef CGAL::cpp0x::tuple<T1,T2,Type> type; };

template <class Type,class T1,class T2,class T3>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,Type> type; };

template <class Type,class T1,class T2,class T3,class T4>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,Type> type; };

template <class Type,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
struct Append_to_tuple<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> >{  typedef CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,Type> type; };

template <class Type,class Tuple,class Result=CGAL::cpp0x::tuple<> >
struct Keep_type_different_of;

template <class Type,class Result>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<>,Result> {
  typedef Result type;
};  
  
template <class Type,class T1,class Result>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type type;
};

template <class Type,class Result,class T1,class T2>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3,T4>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,class T5>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,class T5,class T6>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,class T5,class T6,class T7>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8>,Result> {
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8>,New_result>::type type;
};

template <class Type,class Result,class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
struct Keep_type_different_of<Type,CGAL::cpp0x::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9>,Result>{
  typedef typename boost::mpl::if_< typename boost::is_same<T1,Type>::type,Result,typename Append_to_tuple<T1,Result>::type  >::type New_result;
  typedef typename Keep_type_different_of<Type,CGAL::cpp0x::tuple<T2,T3,T4,T5,T6,T7,T8,T9>,New_result>::type type;
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES


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
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class Functor,int n>
struct Foreach_static{
  template <class  ... T>
  static void run(const T& ... t){
    Functor:: template run<n>(t...);
    Foreach_static<Functor,n-1>::run(t...);
  }
};

template <class Functor>
struct Foreach_static<Functor,0>{
  template <class  ... T>
  static void run(const T& ... t)
  {
    Functor:: template run<0>( t... );
  }
};
#else
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
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5);
  }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6,t7);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6,t7,t8);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8,const T9& t9)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8,t9);
    Foreach_static<Functor,n-1>::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }
};


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
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5);
  }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6,t7);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6,t7,t8);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8,const T9& t9)
  {
    Functor:: template run<0>(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }  
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//Helper function that is calling
//Functor if TAG is different from Void
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class Functor,int n,class Type>
struct Conditionnal_run{
  template <class ... T>
  static void run(const T& ... t){
    Functor:: template run<n>(t...);
  }
};

template <class Functor,int n>
struct Conditionnal_run<Functor,n,Void>
{
  template <class ... T>
  static void run(T...){}
};
#else
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
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5);
  }
  
  template<class T1,class T2,class T3,class T4,class T5,class T6>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8);
  }

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8,const T9& t9)
  {
    Functor:: template run<n>(t1,t2,t3,t4,t5,t6,t7,t8,t9);
  }  
};

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
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&,const T7&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&,const T7&,const T8&){}

  template<class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
  static void run(const T1&,const T2&,const T3&,const T4&,const T5&,const T6&,const T7&,const T8&,const T9&){}
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

//Same as Foreach_static excepted that Functor
//is called for case k only if the k'th type in the tuple
//is different from Void. Note that to the converse of Foreach_static
//Functor are called from n =0 to k
#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
template <class Functor,class T,int n=0>
struct Foreach_static_restricted;

template <class Functor,class Head, class ... Items,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<Head,Items...>,n>
{
  template <class  ... T>
  static void run(const T& ... t){
    Conditionnal_run<Functor,n,Head>::run(t...);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<Items...>,n+1>::run(t...);
  }
};

template <class Functor,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n>{
  template <class  ... T>
  static void run(const T& ... ){}
};
#else
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
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<>,n+1>::run(t1,t2,t3,t4);
  }
};

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
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2>,n+1>::run(t1,t2,t3,t4);
  }
};

template <class Functor,class TT1,class TT2,class TT3,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3>,n+1>::run(t1,t2,t3,t4);
  }
};

template <class Functor,class TT1,class TT2,class TT3,class TT4,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4>,n+1>::run(t1,t2,t3,t4);
  }
};

template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5>,n+1>::run(t1,t2,t3,t4);
  }
};

template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,class TT6,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5,TT6>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6>,n+1>::run(t1,t2,t3,t4);
  }
};


template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,class TT6,class TT7,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7>,n+1>::run(t1,t2,t3,t4);
  }
};

template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,class TT6,class TT7,class TT8,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8>,n+1>::run(t1,t2,t3,t4);
  }
};

template <class Functor,class TT1,class TT2,class TT3,class TT4,class TT5,class TT6,class TT7,class TT8,class TT9,int n>
struct Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT1,TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n>
{
  static void run(){
    Conditionnal_run<Functor,n,TT1>::run();
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run();
  }

  template<class T1>
  static void run(const T1& t1){
    Conditionnal_run<Functor,n,TT1>::run(t1);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run(t1);
  }

  template<class T1,class T2>
  static void run(const T1& t1,const T2& t2){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run(t1,t2);
  }

  template<class T1,class T2,class T3>
  static void run(const T1& t1,const T2& t2,const T3& t3){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run(t1,t2,t3);
  }

  template<class T1,class T2,class T3,class T4>
  static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4){
    Conditionnal_run<Functor,n,TT1>::run(t1,t2,t3,t4);
    Foreach_static_restricted<Functor,CGAL::cpp0x::tuple<TT2,TT3,TT4,TT5,TT6,TT7,TT8,TT9>,n+1>::run(t1,t2,t3,t4);
  }
};


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
};
#endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

template <class CMap>
struct Combinatorial_map_helper
{
  
  // defines as type Compact_container<T>
  template <class T>
  struct Add_compact_container{
    typedef typename CMap::Alloc::template rebind<T>::other Attr_allocator;
    typedef CGAL::Compact_container<T,Attr_allocator> type;
  };
  
  // defines as type Compact_container<T>::iterator
  template <class T>
  struct Add_compact_container_iterator{
    typedef typename CMap::Alloc::template rebind<T>::other Attr_allocator;
    typedef typename CGAL::Compact_container<T,Attr_allocator>::iterator type;
  };
  
  // defines as type Compact_container<T>::const_iterator
  template <class T>
  struct Add_compact_container_const_iterator{
    typedef typename CMap::Alloc::template rebind<T>::other Attr_allocator;
    typedef typename CGAL::Compact_container<T,Attr_allocator>::const_iterator type;
  };

  // All the attributes (with CGAL::Void)
  typedef typename internal::Convert_tuple_with_void<typename CMap::Attributes>::type Attributes;
  // typedef typename CMap::Attributes Attributes;

  // Enabled attributes (without CGAL::Void)
  typedef typename Keep_type_different_of<Void,
    Attributes>::type Enabled_attributes;

  // Number of all attributes
  /*  Does not compile on windows !!
      static const unsigned int number_of_attributes = 
      CGAL::internal::My_length<Attributes>::value; */

  // Number of enabled attributes
  static const unsigned int nb_attribs = 
    Number_of_different_type_in_tuple<Void,Enabled_attributes>::value;
      
  // Given a dimension of the cell, return the index of corresponding attribute
  template <int d>
  struct Dimension_index
  { static const int value=
        Nb_type_different_in_tuple_up_to_k<Void,d,Attributes>::value; };  

  // All these type contains as many entries than the number of enabled attributes
  typedef typename Tuple_converter< Add_compact_container,
				    Enabled_attributes >::type Attribute_containers; 

  typedef typename Tuple_converter< Add_compact_container_iterator,
				    Enabled_attributes >::type  Attribute_iterators; 
  typedef typename Tuple_converter< Add_compact_container_const_iterator,
				    Enabled_attributes >::type  Attribute_const_iterators; 

  typedef Attribute_containers      Attribute_ranges; 
  typedef Attribute_iterators       Attribute_handles;
  typedef Attribute_const_iterators Attribute_const_handles;

  //Helper class allowing to retrieve the type of the attribute of dimension d
  template<int d, int in_tuple=(d<CGAL::internal::My_length<Attributes>::value)>
  struct Attribute_type
  { typedef typename CGAL::cpp0x::tuple_element<d,Attributes>::type type; };

  template<int d>
  struct Attribute_type<d,0>
  { typedef Void type; };

  // Helper class allowing to retreive the d-cell-handle attribute
  template<int d, class Type=typename CGAL::cpp0x::tuple_element<d,Attributes>::type>
  struct Attribute_handle
  {
    typedef typename CGAL::cpp0x::tuple_element<Dimension_index<d>::value,Attribute_handles>::type 
      type;
  };
  
  template<int d>
  struct Attribute_handle<d,Void>
  { typedef Void type; };

  // Helper class allowing to retreive the d-cell-const handle attribute
  template<int d, class Type=typename CGAL::cpp0x::tuple_element<d,Attributes>::type  >
  struct Attribute_const_handle
  {
    typedef typename CGAL::cpp0x::tuple_element<Dimension_index<d>::value,
						Attribute_const_handles>::type 
    type;
  };
  
  template<int d>
  struct Attribute_const_handle<d,Void>
  { typedef Void type; };

  // Helper class allowing to retreive the d-cell-iterator attribute
  template<int d, class Type=typename CGAL::cpp0x::tuple_element<d,Attributes>::type>
  struct Attribute_iterator
  {
    typedef typename CGAL::cpp0x::tuple_element<Dimension_index<d>::value,
						Attribute_iterators>::type 
      type;
  };
  
  template<int d>
  struct Attribute_iterator<d,Void>
  { typedef Void type; };

  // Helper class allowing to retreive the d-cell-const handle attribute
  template<int d, class Type=typename CGAL::cpp0x::tuple_element<d,Attributes>::type  >
  struct Attribute_const_iterator
  {
    typedef typename CGAL::cpp0x::tuple_element<Dimension_index<d>::value,
						Attribute_const_iterators>::type 
      type;
  };
  
  template<int d>
  struct Attribute_const_iterator<d,Void>
  { typedef Void type; };

  // Helper class allowing to retreive the d-cell-attribute range
  template<int d, class Type=
	   typename CGAL::cpp0x::tuple_element<d,Attributes>::type>
  struct Attribute_range
  {
    typedef typename CGAL::cpp0x::tuple_element<Dimension_index<d>::value,
						Attribute_ranges>::type type;
  };
  
  template<int d>
  struct Attribute_range<d,Void>
  { typedef Void type; };

  // Helper class allowing to retreive the d-cell-attribute range
  template<int d, class Type=typename CGAL::cpp0x::tuple_element<d,Attributes>::type  >
  struct Attribute_const_range
  {
    typedef const typename CGAL::cpp0x::tuple_element<Dimension_index<d>::value,
						      Attribute_ranges >::type 
    type;
  };
  
  template<int d>
  struct Attribute_const_range<d,Void>
  { typedef Void type; };

  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  // To iterate onto each dimension of the map
  template <class Functor>
  struct Foreach_dimension
  {
    template <class ...Ts>
    static void run(const Ts& ... t)
    { Foreach_static<Functor,CMap::dimension >::run(t...); }
  };
  #endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

  // To iterate onto each enabled attributes
  #ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <class Functor>
  struct Foreach_enabled_attributes
  {
    template <class ...Ts>
    static void run(const Ts& ... t)
    { Foreach_static_restricted<Functor,Attributes >::run(t...); }
  };
  #else
  template <class Functor>
  struct Foreach_enabled_attributes
  {
    static void run() {Foreach_static_restricted<Functor,Attributes >::run();}
    
    template <class T1>
    static void run(const T1& t1) {Foreach_static_restricted<Functor,Attributes >::run(t1);}

    template <class T1,class T2>
    static void run(const T1& t1,const T2& t2) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2);}

    template <class T1,class T2,class T3>
    static void run(const T1& t1,const T2& t2,const T3& t3) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3);}

    template <class T1,class T2,class T3,class T4>
    static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4);}

    template <class T1,class T2,class T3,class T4,class T5>
    static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5);}            

    template <class T1,class T2,class T3,class T4,class T5,class T6>
    static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,t6);}

    template <class T1,class T2,class T3,class T4,class T5,class T6,class T7>
    static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,t6,t7);}

    template <class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8>
    static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,t6,t7,t8);}

    template <class T1,class T2,class T3,class T4,class T5,class T6,class T7,class T8,class T9>
    static void run(const T1& t1,const T2& t2,const T3& t3,const T4& t4,const T5& t5,const T6& t6,const T7& t7,const T8& t8,const T9& t9) {Foreach_static_restricted<Functor,Attributes >::run(t1,t2,t3,t4,t5,t6,t7,t8,t9);}
  };
  #endif //CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
};


} //namespace CGAL

} //namespace internal

#endif //CGAL_INTERNAL_COMBINATORIAL_MAP_UTILITY_H
