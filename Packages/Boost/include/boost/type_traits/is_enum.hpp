
//  (C) Copyright Dave Abrahams, Steve Cleary, Beman Dawes, Howard
//  Hinnant & John Maddock 2000.  
//  Use, modification and distribution are subject to the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt).
//
//  See http://www.boost.org/libs/type_traits for most recent version including documentation.


#ifndef BOOST_TT_IS_ENUM_HPP_INCLUDED
#define BOOST_TT_IS_ENUM_HPP_INCLUDED

#include "boost/type_traits/add_reference.hpp"
#include "boost/type_traits/is_arithmetic.hpp"
#include "boost/type_traits/is_reference.hpp"
#include "boost/type_traits/is_convertible.hpp"
#ifdef __GNUC__
#include <boost/type_traits/is_function.hpp>
#endif
#include "boost/type_traits/config.hpp"
#if defined(BOOST_TT_HAS_CONFORMING_IS_CLASS_IMPLEMENTATION) 
#  include <boost/type_traits/is_class.hpp>
#  include <boost/type_traits/is_union.hpp>
#endif


// should be the last #include
#include "boost/type_traits/detail/bool_trait_def.hpp"

namespace boost {

#if !(defined(__BORLANDC__) && (__BORLANDC__ <= 0x551))

namespace detail {

#if defined(BOOST_TT_HAS_CONFORMING_IS_CLASS_IMPLEMENTATION) 

template <typename T>
struct is_class_or_union
{
   BOOST_STATIC_CONSTANT(bool, value =
      (::boost::type_traits::ice_or<
           ::boost::is_class<T>::value
         , ::boost::is_union<T>::value
      >::value));
};

#else

template <typename T>
struct is_class_or_union
{
# if BOOST_WORKAROUND(BOOST_MSVC, == 1200) || BOOST_WORKAROUND(__BORLANDC__, <= 0x570)// we simply can't detect it this way.
    BOOST_STATIC_CONSTANT(bool, value = false);
# else
    template <class U> static ::boost::type_traits::yes_type is_class_or_union_tester(void(U::*)(void));

#  if BOOST_WORKAROUND(BOOST_MSVC, == 1300)                 \
    || BOOST_WORKAROUND(__MWERKS__, <= 0x3000) // no SFINAE
    static ::boost::type_traits::no_type is_class_or_union_tester(...);
    BOOST_STATIC_CONSTANT(
        bool, value = sizeof(is_class_or_union_tester(0)) == sizeof(::boost::type_traits::yes_type));
#  else
    template <class U>
    static ::boost::type_traits::no_type is_class_or_union_tester(...);
    BOOST_STATIC_CONSTANT(
        bool, value = sizeof(is_class_or_union_tester<T>(0)) == sizeof(::boost::type_traits::yes_type));
#  endif
# endif
};
#endif

struct int_convertible
{
    int_convertible(int);
};

// Don't evaluate convertibility to int_convertible unless the type
// is non-arithmetic. This suppresses warnings with GCC.
template <bool is_typename_arithmetic_or_reference = true>
struct is_enum_helper
{
    template <typename T> struct type
    {
        BOOST_STATIC_CONSTANT(bool, value = false);
    };
};

template <>
struct is_enum_helper<false>
{
    template <typename T> struct type
        : ::boost::is_convertible<T,::boost::detail::int_convertible>
    {
    };
};

template <typename T> struct is_enum_impl
{
   typedef ::boost::add_reference<T> ar_t;
   typedef typename ar_t::type r_type;

#if defined(__GNUC__)

#ifdef BOOST_TT_HAS_CONFORMING_IS_CLASS_IMPLEMENTATION
   BOOST_STATIC_CONSTANT(bool, selector =
      (::boost::type_traits::ice_or<
           ::boost::is_arithmetic<T>::value
         , ::boost::is_reference<T>::value
         , ::boost::is_function<T>::value
         , is_class_or_union<T>::value
      >::value));
#else
   BOOST_STATIC_CONSTANT(bool, selector =
      (::boost::type_traits::ice_or<
           ::boost::is_arithmetic<T>::value
         , ::boost::is_reference<T>::value
         , ::boost::is_function<T>::value
      >::value));
#endif // BOOST_TT_HAS_CONFORMING_IS_CLASS_IMPLEMENTATION

#else
   BOOST_STATIC_CONSTANT(bool, selector =
      (::boost::type_traits::ice_or<
           ::boost::is_arithmetic<T>::value
         , ::boost::is_reference<T>::value
         , is_class_or_union<T>::value
       // However, not doing this on non-conforming compilers prevents
       // a dependency recursion.
      >::value));
#endif

#if BOOST_WORKAROUND(__BORLANDC__, < 0x600)
    typedef ::boost::detail::is_enum_helper<
          ::boost::detail::is_enum_impl<T>::selector
        > se_t;
#else
    typedef ::boost::detail::is_enum_helper<selector> se_t;
#endif

    typedef typename se_t::template type<r_type> helper;
    BOOST_STATIC_CONSTANT(bool, value = helper::value);
};

// these help on compilers with no partial specialization support:
BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(is_enum,void,false)
#ifndef BOOST_NO_CV_VOID_SPECIALIZATIONS
BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(is_enum,void const,false)
BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(is_enum,void volatile,false)
BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(is_enum,void const volatile,false)
#endif

} // namespace detail

BOOST_TT_AUX_BOOL_TRAIT_DEF1(is_enum,T,::boost::detail::is_enum_impl<T>::value)

#else // __BORLANDC__
//
// buggy is_convertible prevents working
// implementation of is_enum:
BOOST_TT_AUX_BOOL_TRAIT_DEF1(is_enum,T,false)

#endif

} // namespace boost

#include "boost/type_traits/detail/bool_trait_undef.hpp"

#endif // BOOST_TT_IS_ENUM_HPP_INCLUDED
