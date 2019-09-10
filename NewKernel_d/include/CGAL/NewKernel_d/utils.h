// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_MARCUTILS
#define CGAL_MARCUTILS

#include <CGAL/config.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro 'BOOST_PP_EXPAND_I'
                                // http://lists.boost.org/boost-users/2014/11/83291.php
#endif                          

#ifdef CGAL_CXX11
#include <type_traits>
#include <utility>
#define CGAL_FORWARDABLE(T) T&&
#define CGAL_FORWARD(T,t) std::forward<T>(t)
#define CGAL_MOVE(t) std::move(t)
#define CGAL_CONSTEXPR constexpr
#else
#define CGAL_FORWARDABLE(T) T const&
#define CGAL_FORWARD(T,t) t
#define CGAL_MOVE(t) t
#define CGAL_CONSTEXPR
#endif
#include <boost/utility/enable_if.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <CGAL/Rational_traits.h>
#include <CGAL/tuple.h>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/not.hpp>
#include <boost/type_traits.hpp>

#ifdef CGAL_CXX11
#define CGAL_BOOSTD std::
#else
#define CGAL_BOOSTD boost::
#endif

namespace CGAL {
namespace internal {
	BOOST_MPL_HAS_XXX_TRAIT_DEF(type)
}

template <class T, class No, bool=internal::has_type<T>::value /*false*/>
struct Has_type_different_from : boost::false_type {};
template <class T, class No>
struct Has_type_different_from <T, No, true>
: boost::mpl::not_<boost::is_same<typename T::type, No> > {};


	template <class T> struct Wrap_type { typedef T type; };

	// tell a function f(a,b,c) that its real argument is a(b,c)
	struct Eval_functor {};

	// forget the first argument. Useful to make something dependant
	// (and thus usable in SFINAE), although that's not a great design.
	template<class A,class B> struct Second_arg {
		typedef B type;
	};

	// like std::forward, except for basic types where it does a cast, to
	// avoid issues with narrowing conversions
#ifdef CGAL_CXX11
	template<class T,class U,class V> inline
		typename std::conditional<std::is_arithmetic<T>::value&&std::is_arithmetic<typename std::remove_reference<U>::type>::value,T,U&&>::type
	       	forward_safe(V&& u) { return std::forward<U>(u); }
#else
	template<class T,class U> inline U const& forward_safe(U const& u) {
		return u;
	}
#endif

#ifdef CGAL_CXX11
	template<class...> struct Constructible_from_each;
	template<class To,class From1,class...From> struct Constructible_from_each<To,From1,From...>{
		enum { value=std::is_convertible<From1,To>::value&&Constructible_from_each<To,From...>::value };
	};
	template<class To> struct Constructible_from_each<To>{
		enum { value=true };
	};
#else
// currently only used in C++0X code
#endif

	template<class T> struct Scale {
#ifndef CGAL_CXX11
		template<class> struct result;
		template<class FT> struct result<Scale(FT)> {
			typedef FT type;
		};
#endif
		T const& scale;
		Scale(T const& t):scale(t){}
		template<class FT>
#ifdef CGAL_CXX11
		auto operator()(FT&& x)const->decltype(scale*std::forward<FT>(x))
#else
		FT operator()(FT const& x)const
#endif
		{
			return scale*CGAL_FORWARD(FT,x);
		}
	};
	template<class NT,class T> struct Divide {
#if !defined(CGAL_CXX11) || !defined(BOOST_RESULT_OF_USE_DECLTYPE)
		// requires boost > 1.44
		// shouldn't be needed with C++0X
		//template<class> struct result;
		//template<class FT> struct result<Divide(FT)> {
		//	typedef FT type;
		//};
		typedef NT result_type;
#endif
		T const& scale;
		Divide(T const& t):scale(t){}
		template<class FT>
#ifdef CGAL_CXX11
		//FIXME: gcc complains for Gmpq
		//auto operator()(FT&& x)const->decltype(Rational_traits<NT>().make_rational(std::forward<FT>(x),scale))
		NT operator()(FT&& x)const
#else
		NT operator()(FT const& x)const
#endif
		{
			return Rational_traits<NT>().
				make_rational(CGAL_FORWARD(FT,x),scale);
		}
	};

	template <class NT> struct has_cheap_constructor : boost::is_arithmetic<NT>{};
	template <bool p> struct has_cheap_constructor<Interval_nt<p> > {
		        enum { value=true };
	};

	// like std::multiplies but allows mixing types
	// in C++11 in doesn't need to be a template
	template < class Ret >
	struct multiplies {
		template<class A,class B>
#ifdef CGAL_CXX11
		auto operator()(A&&a,B&&b)const->decltype(std::forward<A>(a)*std::forward<B>(b))
#else
		Ret operator()(A const& a, B const& b)const
#endif
		{
			return CGAL_FORWARD(A,a)*CGAL_FORWARD(B,b);
		}
	};
	template < class Ret >
	struct division {
		template<class A,class B>
#ifdef CGAL_CXX11
		auto operator()(A&&a,B&&b)const->decltype(std::forward<A>(a)/std::forward<B>(b))
#else
		Ret operator()(A const& a, B const& b)const
#endif
		{
			return CGAL_FORWARD(A,a)/CGAL_FORWARD(B,b);
		}
	};

#ifdef CGAL_CXX11
	using std::decay;
#else
	template<class T> struct decay : boost::remove_cv<typename boost::decay<T>::type> {};
#endif

	template<class T,class U> struct Type_copy_ref { typedef U type; };
	template<class T,class U> struct Type_copy_ref<T&,U> { typedef U& type; };
#ifdef CGAL_CXX11
	template<class T,class U> struct Type_copy_ref<T&&,U> { typedef U&& type; };
#endif
	template<class T,class U> struct Type_copy_cv { typedef U type; };
	template<class T,class U> struct Type_copy_cv<T const,U> { typedef U const type; };
	template<class T,class U> struct Type_copy_cv<T volatile,U> { typedef U volatile type; };
	template<class T,class U> struct Type_copy_cv<T const volatile,U> { typedef U const volatile type; };

	template<class T,class U> struct Type_copy_cvref :
		Type_copy_ref<T,typename Type_copy_cv<typename boost::remove_reference<T>::type,U>::type> {};

	struct Dereference_functor {
		template<class> struct result{};
		template<class It> struct result<Dereference_functor(It)> {
			typedef typename std::iterator_traits<It>::reference type;
		};
		template<class It> typename result<Dereference_functor(It)>::type
			operator()(It const&i)const{
				return *i;
			}
	};

#ifdef CGAL_CXX11
	template<int...> struct Indices{};
	template<class> struct Next_increasing_indices;
	template<int...I> struct Next_increasing_indices<Indices<I...> > {
		typedef Indices<I...,sizeof...(I)> type;
	};
	template<int N> struct N_increasing_indices {
		typedef typename Next_increasing_indices<typename N_increasing_indices<N-1>::type>::type type;
	};
	template<> struct N_increasing_indices<0> { typedef Indices<> type; };
	namespace internal {
	template<class F,class...U,int...I> inline typename std::result_of<F&&(U...)>::type
	do_call_on_tuple_elements(F&&f, std::tuple<U...>&&t, Indices<I...>&&) {
		return f(std::get<I>(std::move(t))...);
	}
	} // internal
	template<class/*result type, ignored*/,class F,class...U>
	inline typename std::result_of<F&&(U...)>::type
	call_on_tuple_elements(F&&f, std::tuple<U...>&&t) {
		return internal::do_call_on_tuple_elements(std::forward<F>(f),std::move(t),
				typename N_increasing_indices<sizeof...(U)>::type());
	}
#else
#define CGAL_VAR(Z,N,_) cpp0x::get<N>(t)
#define CGAL_CODE(Z,N,_) template<class Res, class F BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N,class U)> \
	inline Res call_on_tuple_elements(F const&f, \
			cpp0x::tuple<BOOST_PP_ENUM_PARAMS(N,U)> const&t) { \
		return f(BOOST_PP_ENUM(N,CGAL_VAR,)); \
	}
	template<class Res, class F>
	inline Res call_on_tuple_elements(F const&f, cpp0x::tuple<>) {
		return f();
	}
BOOST_PP_REPEAT_FROM_TO(1, 8, CGAL_CODE, _ )
#undef CGAL_CODE
#undef CGAL_VAR
#endif

	template<class A> struct Factory {
	  typedef A result_type;
#ifdef CGAL_CXX11
	  template<class...U> result_type operator()(U&&...u)const{
	    return A(std::forward<U>(u)...);
	  }
#else
	  result_type operator()()const{
	    return A();
	  }
#define CGAL_CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class U)> \
	  result_type operator()(BOOST_PP_ENUM_BINARY_PARAMS(N,U,const&u))const{ \
	    return A(BOOST_PP_ENUM_PARAMS(N,u)); \
	  }
BOOST_PP_REPEAT_FROM_TO(1, 8, CGAL_CODE, _ )
#undef CGAL_CODE
#endif
	};
}

// TODO: make a Cartesian-only variant
// WARNING: do not use the Req* parameters too much, they can cause circular instanciations and are only useful for dispatching.
#define CGAL_STRIP_PAREN_(...) __VA_ARGS__
#define CGAL_STRIP_PAREN(...) CGAL_STRIP_PAREN_ __VA_ARGS__
// What to do with O? pass it down to other functors or drop it?
#define CGAL_KD_DEFAULT_FUNCTOR(Tg,Name,ReqTyp,ReqFun) \
    template <class K, class O> \
    struct Get_functor<K, Tg, O, \
      typename boost::mpl::if_c< \
        Provides_functor_i<K, Tg, O>::value \
        || !Provides_types<K, boost::mpl::vector<CGAL_STRIP_PAREN_ ReqTyp> >::value \
        || !Provides_functors<K, boost::mpl::vector<CGAL_STRIP_PAREN_ ReqFun> >::value \
      , int, void>::type> \
    { \
      typedef CGAL_STRIP_PAREN_ Name type; \
      typedef K Bound_kernel; \
    }

// Not used yet, may need some changes.
#define CGAL_KD_DEFAULT_TYPE(Tg,Name,ReqTyp,ReqFun) \
    template <class K> \
    struct Get_type<K, Tg, \
      typename boost::mpl::if_c< \
        Provides_type_i<K, Tg>::value \
        || !Provides_types<K, boost::mpl::vector<CGAL_STRIP_PAREN_ ReqTyp> >::value \
        || !Provides_functors<K, boost::mpl::vector<CGAL_STRIP_PAREN_ ReqFun> >::value \
      , int, void>::type> \
    { \
      typedef CGAL_STRIP_PAREN_ Name type; \
      typedef K Bound_kernel; \
    }

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif
