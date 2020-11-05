// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_MARCUTILS
#define CGAL_MARCUTILS

#include <CGAL/config.h>

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4003) // not enough actual parameters for macro 'BOOST_PP_EXPAND_I'
                                // https://lists.boost.org/boost-users/2014/11/83291.php
#endif

#include <type_traits>
#include <utility>
#include <boost/utility/enable_if.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <CGAL/Rational_traits.h>
#include <CGAL/tuple.h>
#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/not.hpp>
#include <boost/type_traits.hpp>

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
        template<class T,class U,class V> inline
                typename std::conditional<std::is_arithmetic<T>::value&&std::is_arithmetic<typename std::remove_reference<U>::type>::value,T,U&&>::type
                       forward_safe(V&& u) { return std::forward<U>(u); }

        template<class...> struct Constructible_from_each;
        template<class To,class From1,class...From> struct Constructible_from_each<To,From1,From...>{
                enum { value=std::is_convertible<From1,To>::value&&Constructible_from_each<To,From...>::value };
        };
        template<class To> struct Constructible_from_each<To>{
                enum { value=true };
        };

        template<class T> struct Scale {
                T const& scale;
                Scale(T const& t):scale(t){}
                template<class FT>
                decltype(auto) operator()(FT&& x)const
                {
                        return (scale*std::forward<FT>(x));
                }
        };
        template<class NT,class T> struct Divide {
                T const& scale;
                Divide(T const& t):scale(t){}
                template<class FT>
                //FIXME: gcc complains for Gmpq
                //decltype(auto) operator()(FT&& x)const
                NT operator()(FT&& x)const
                {
                        return Rational_traits<NT>().
                                make_rational(std::forward<FT>(x),scale);
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
                decltype(auto) operator()(A&&a,B&&b)const
                {
                        return std::forward<A>(a)*std::forward<B>(b);
                }
        };
        template < class Ret >
        struct division {
                template<class A,class B>
                decltype(auto) operator()(A&&a,B&&b)const
                {
                        return std::forward<A>(a)/std::forward<B>(b);
                }
        };

        using std::decay;

        template<class T,class U> struct Type_copy_ref { typedef U type; };
        template<class T,class U> struct Type_copy_ref<T&,U> { typedef U& type; };
        template<class T,class U> struct Type_copy_ref<T&&,U> { typedef U&& type; };
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
                template<class It> decltype(auto)
                        operator()(It const&i)const{
                                return *i;
                        }
        };

        namespace internal {
        template<class F,class...U,std::size_t...I> inline decltype(auto)
        do_call_on_tuple_elements(F&&f, std::tuple<U...>&&t, std::index_sequence<I...>&&) {
                return f(std::get<I>(std::move(t))...);
        }
        } // internal
        template<class F,class...U>
        inline decltype(auto)
        call_on_tuple_elements(F&&f, std::tuple<U...>&&t) {
                return internal::do_call_on_tuple_elements(std::forward<F>(f),std::move(t),
                                std::make_index_sequence<sizeof...(U)>());
        }
        template<class A> struct Factory {
          typedef A result_type;
          template<class...U> result_type operator()(U&&...u)const{
            return A(std::forward<U>(u)...);
          }
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
