#ifndef BOOST_TT_IS_ABSTRACT_CLASS_HPP
#define BOOST_TT_IS_ABSTRACT_CLASS_HPP

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// is_abstract_class.hpp:
//
//  (C) Copyright 2002 Rani Sharoni (rani_sharoni@hotmail.com) and Robert Ramey
//  Use, modification and distribution is subject to the Boost Software
//  License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//  
//  See http://www.boost.org for updates, documentation, and revision history.
//

// Compile type discovery whether given type is abstract class or not.
//
//   Requires DR 337 to be supported by compiler
//   (http://anubis.dkuug.dk/jtc1/sc22/wg21/docs/cwg_active.html#337).
//
//
// Believed (Jan 2004) to work on:
//  - GCC 3.4
//  - VC++ 7.1
//  - compilers with new EDG frontend (Intel C++ 7, Comeau 4.3.2)
//
// Doesn't work on:
//  - VC++6, VC++7.0 and less
//  - GCC 3.3.X and less
//  - Borland C++ 6 and less
//      
//
// History:
//  - Originally written by Rani Sharoni, see
//    http://groups.google.com/groups?selm=df893da6.0207110613.75b2fe90%40posting.google.com
//    At this time supported by EDG (Intel C++ 7, Comeau 4.3.2) and VC7.1.
//  - Adapted and added into Boost.Serialization library by Robert Ramey 
//    (starting with submission #10).
//  - Jan 2004: GCC 3.4 fixed to suport DR337 (Giovanni Bajo).
//  - Jan 2004: modified to be part of Boost.TypeTraits (Pavel Vozenilek).
//

#include <boost/type_traits/detail/yes_no_type.hpp>
#include <boost/type_traits/is_class.hpp>
#include "boost/type_traits/detail/ice_and.hpp"
// should be the last #include
#include "boost/type_traits/detail/bool_trait_def.hpp"


namespace boost {
namespace detail{

template<class T>
struct is_abstract_imp
{
   // Deduction fails if T is void, function type, 
   // reference type (14.8.2/2)or an abstract class type 
   // according to review status issue #337
   //
   template<class U>
   static type_traits::no_type check_sig(U (*)[1]);
   template<class U>
   static type_traits::yes_type check_sig(...);

   // GCC2 won't even parse this template if we embed the computation
   // of s1 in the computation of value.
#ifdef __GNUC__
   BOOST_STATIC_CONSTANT(unsigned, s1 = sizeof(is_abstract_imp<T>::template check_sig<T>(0)));
#else
   BOOST_STATIC_CONSTANT(unsigned, s1 = sizeof(check_sig<T>(0)));
#endif
    
   BOOST_STATIC_CONSTANT(bool, value = 
      (::boost::type_traits::ice_and<
         ::boost::is_class<T>::value,
         (s1 == sizeof(type_traits::yes_type))
      >::value));
};

}

BOOST_TT_AUX_BOOL_TRAIT_DEF1(is_abstract,T,::boost::detail::is_abstract_imp<T>::value)

} // namespace boost

#include "boost/type_traits/detail/bool_trait_undef.hpp"

#endif //BOOST_TT_IS_ABSTRACT_CLASS_HPP
