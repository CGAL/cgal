/*
 *
 * Copyright (c) 1998-2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         instances.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Defines those template instances that are placed in the
  *                library rather than in the users object files.
  */

//
// note no include guard, we may include this multiple times:
//
#ifndef BOOST_REGEX_NO_EXTERNAL_TEMPLATES

namespace boost{

//
// this header can be included multiple times, each time with
// a different character type, BOOST_REGEX_CHAR_T must be defined
// first:
//
#ifndef BOOST_REGEX_CHAR_T
#  error "BOOST_REGEX_CHAR_T not defined"
#endif

//
// what follows is compiler specific:
//

#if  defined(__BORLANDC__) && (__BORLANDC__ < 0x600)

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

#  ifndef BOOST_REGEX_INSTANTIATE
#     pragma option push -Jgx
#  endif

template class BOOST_REGEX_DECL reg_expression< BOOST_REGEX_CHAR_T >;

#  ifndef BOOST_REGEX_INSTANTIATE
#     pragma option pop
#  endif

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#elif defined(BOOST_MSVC) || defined(__GNUC__)

#  ifndef BOOST_REGEX_INSTANTIATE
#     define template extern template
#  endif

#  ifdef BOOST_MSVC
#     pragma warning(push)
#     pragma warning(disable : 4251 4231 4660)
#  endif

template class BOOST_REGEX_DECL reg_expression< BOOST_REGEX_CHAR_T >;

#  ifdef BOOST_MSVC
#     pragma warning(pop)
#  endif

#  ifdef template
#     undef template
#  endif

#endif

} // namespace boost

#endif // BOOST_REGEX_NO_EXTERNAL_TEMPLATES




