/*
 *
 * Copyright (c) 2002
 * Dr John Maddock
 *
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */

 /*
  *   LOCATION:    see http://www.boost.org for most recent version.
  *   FILE         char_regex_traits.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Declares deprecated traits classes char_regex_traits<>.
  */


#ifndef BOOST_REGEX_V4_CHAR_REGEX_TRAITS_HPP
#define BOOST_REGEX_V4_CHAR_REGEX_TRAITS_HPP

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost{

namespace deprecated{
//
// class char_regex_traits_i
// provides case insensitive traits classes (deprecated):
template <class charT>
class char_regex_traits_i : public regex_traits<charT> {};

template<>
class char_regex_traits_i<char> : public regex_traits<char>
{
public:
   typedef char char_type;
   typedef unsigned char uchar_type;
   typedef unsigned int size_type;
   typedef regex_traits<char> base_type;

   char BOOST_REGEX_CALL translate(char c, bool)const
   {
      return static_cast<const regex_traits<char>*>(this)->translate(c, true);
   }
};

#ifndef BOOST_NO_WREGEX
template<>
class char_regex_traits_i<wchar_t> : public regex_traits<wchar_t>
{
public:
   typedef wchar_t char_type;
   typedef unsigned short uchar_type;
   typedef unsigned int size_type;
   typedef regex_traits<wchar_t> base_type;

   wchar_t BOOST_REGEX_CALL translate(wchar_t c, bool)const
   {
      return static_cast<const regex_traits<wchar_t>*>(this)->translate(c, true);
   }
   boost::uint_fast32_t BOOST_REGEX_CALL lookup_classname(const wchar_t* first, const wchar_t* last)const
   {
      boost::uint_fast32_t result = static_cast<const regex_traits<wchar_t>*>(this)->lookup_classname(first, last);
      if((result & base_type::char_class_upper) == base_type::char_class_upper)
         result |= base_type::char_class_alpha;
      return result;
   }
};
#endif
} // namespace deprecated
} // namespace boost

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#endif // include

