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
  *   FILE         regbase.cpp
  *   VERSION      see <boost/version.hpp>
  *   DESCRIPTION: Declares class regbase.
  */

#ifndef BOOST_REGEX_V4_REGBASE_HPP
#define BOOST_REGEX_V4_REGBASE_HPP

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost{
//
// class regbase
// handles error codes and flags
//
class BOOST_REGEX_DECL regbase
{
public:
   enum flag_type_
   {
      escape_in_lists = 1,                     // '\' special inside [...]
      char_classes = escape_in_lists << 1,     // [[:CLASS:]] allowed
      intervals = char_classes << 1,           // {x,y} allowed
      limited_ops = intervals << 1,            // all of + ? and | are normal characters
      newline_alt = limited_ops << 1,          // \n is the same as |
      bk_plus_qm = newline_alt << 1,           // uses \+ and \?
      bk_braces = bk_plus_qm << 1,             // uses \{ and \}
      bk_parens = bk_braces << 1,              // uses \( and \)
      bk_refs = bk_parens << 1,                // \d allowed
      bk_vbar = bk_refs << 1,                  // uses \|

      use_except = bk_vbar << 1,               // exception on error
      failbit = use_except << 1,               // error flag
      literal = failbit << 1,                  // all characters are literals
      icase = literal << 1,                    // characters are matched regardless of case
      nocollate = 0,                           // don't use locale specific collation (deprecated)
      collate = icase << 1,                    // use locale specific collation
      perlex = collate << 1,                 // perl extensions
      nosubs = perlex << 1,                    // don't mark sub-expressions
      optimize = 0,                            // not really supported

      basic = char_classes | intervals | limited_ops | bk_braces | bk_parens | bk_refs | collate,
      extended = char_classes | intervals | bk_refs | collate,
      normal = perlex | escape_in_lists | char_classes | intervals | bk_refs | nocollate,
      emacs = bk_braces | bk_parens | bk_refs | bk_vbar,
      awk = extended | escape_in_lists,
      grep = basic | newline_alt,
      egrep = extended | newline_alt,
      sed = basic,
      perl = normal,
      ECMAScript = normal,
      JavaScript = normal,
      JScript = normal
   };
   typedef unsigned int flag_type;

   enum restart_info
   {
      restart_any = 0,
      restart_word = 1,
      restart_line = 2,
      restart_buf = 3,
      restart_continue = 4,
      restart_lit = 5,
      restart_fixed_lit = 6, 
      restart_count = 7
   };

   flag_type BOOST_REGEX_CALL flags()const
   {
      return _flags;
   }

   regbase();
   regbase(const regbase& b);
   void swap(regbase& that)
   { std::swap(_flags, that._flags); }
protected:
   flag_type _flags;
};

//
// provide std lib proposal compatible constants:
//
namespace regex_constants{

   enum flag_type_
   {
      escape_in_lists = ::boost::regbase::escape_in_lists,
      char_classes = ::boost::regbase::char_classes,
      intervals = ::boost::regbase::intervals,
      limited_ops = ::boost::regbase::limited_ops,
      newline_alt = ::boost::regbase::newline_alt,
      bk_plus_qm = ::boost::regbase::bk_plus_qm,
      bk_braces = ::boost::regbase::bk_braces,
      bk_parens = ::boost::regbase::bk_parens,
      bk_refs = ::boost::regbase::bk_refs,
      bk_vbar = ::boost::regbase::bk_vbar,

      use_except = ::boost::regbase::use_except,
      failbit = ::boost::regbase::failbit,
      literal = ::boost::regbase::literal,
      icase = ::boost::regbase::icase,
      nocollate = ::boost::regbase::nocollate,
      collate = ::boost::regbase::collate,
      perlex = ::boost::regbase::perlex,
      nosubs = ::boost::regbase::nosubs,
      optimize = ::boost::regbase::optimize,

      basic = ::boost::regbase::basic,
      extended = ::boost::regbase::extended,
      normal = ::boost::regbase::normal,
      emacs = ::boost::regbase::emacs,
      awk = ::boost::regbase::awk,
      grep = ::boost::regbase::grep,
      egrep = ::boost::regbase::egrep,
      sed = basic,
      perl = normal,
      ECMAScript = normal,
      JavaScript = normal,
      JScript = normal
   };
   typedef ::boost::regbase::flag_type syntax_option_type;

} // namespace regex_constants

} // namespace boost

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#endif

