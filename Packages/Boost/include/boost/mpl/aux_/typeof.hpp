//-----------------------------------------------------------------------------
// boost mpl/aux_/typeof.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2002
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_TYPEOF_HPP_INCLUDED
#define BOOST_MPL_AUX_TYPEOF_HPP_INCLUDED
#include <boost/detail/workaround.hpp>

#if defined(__BORLANDC__) 
#   define BOOST_MPL_AUX_TYPEOF(T,x) typename T::value_type
#elif BOOST_WORKAROUND(__MWERKS__, <= 0x2407) || BOOST_WORKAROUND(BOOST_MSVC, <= 1300) || BOOST_WORKAROUND(__EDG_VERSION__, <= 243)
#   define BOOST_MPL_AUX_TYPEOF(T,x) long
#elif defined(__GCC__) && !BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3003))
#   define BOOST_MPL_AUX_TYPEOF(T,x) __typeof__(x)
#else
#   include "boost/config.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"

namespace boost {
namespace mpl {

// the implementation below is based on "A Portable typeof Operator" article
// by Bill Gibbons, C++ User Journal, November 2000

namespace aux {
template< BOOST_MPL_AUX_NTTP_DECL(long, N) > struct typeof_answer { typedef char type[N]; };
template< BOOST_MPL_AUX_NTTP_DECL(long, S) > struct typeof_c;
}

#define BOOST_MPL_AUX_REGISTER_TYPE(index, T) \
namespace boost { namespace mpl { namespace aux { \
template<> struct typeof_c<index> { typedef T type; }; \
typeof_answer<index>::type& type_index(T const&); \
}}} \
/**/

#define BOOST_MPL_AUX_TYPEOF(T,x) \
typename boost::mpl::aux::typeof_c< \
    sizeof(::boost::mpl::aux::type_index(x)) \
    >::type \
/**/

} // namespace mpl
} // namespace boost

BOOST_MPL_AUX_REGISTER_TYPE(1, bool)
BOOST_MPL_AUX_REGISTER_TYPE(2, signed char)
BOOST_MPL_AUX_REGISTER_TYPE(3, unsigned char)
BOOST_MPL_AUX_REGISTER_TYPE(4, char)
#if !defined(BOOST_NO_INTRINSIC_WCHAR_T)
BOOST_MPL_AUX_REGISTER_TYPE(5, wchar_t)
#endif
BOOST_MPL_AUX_REGISTER_TYPE(6, short)
BOOST_MPL_AUX_REGISTER_TYPE(7, unsigned short)
BOOST_MPL_AUX_REGISTER_TYPE(8, int)
BOOST_MPL_AUX_REGISTER_TYPE(9, unsigned int)
BOOST_MPL_AUX_REGISTER_TYPE(10, long)
BOOST_MPL_AUX_REGISTER_TYPE(11, unsigned long)
//BOOST_MPL_AUX_REGISTER_TYPE(12, float)
//BOOST_MPL_AUX_REGISTER_TYPE(13, double)
//BOOST_MPL_AUX_REGISTER_TYPE(14, long double)

#endif // __GCC__

#endif // BOOST_MPL_AUX_TYPEOF_HPP_INCLUDED
