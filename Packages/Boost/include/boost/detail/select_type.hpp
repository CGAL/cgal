// (C) Copyright David Abrahams 2001. Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
//
// See http://www.boost.org for most recent version including documentation.

// Revision History
// 09 Feb 01  Applied John Maddock's Borland patch Moving <true>
//            specialization to unspecialized template (David Abrahams)
// 06 Feb 01  Created (David Abrahams)

#ifndef SELECT_TYPE_DWA20010206_HPP
# define SELECT_TYPE_DWA20010206_HPP

namespace boost { namespace detail {

  // Template class if_true -- select among 2 types based on a bool constant expression
  // Usage:
  //   typename if_true<(bool_const_expression)>::template then<true_type, false_type>::type

  // HP aCC cannot deal with missing names for template value parameters
  template <bool b> struct if_true
  {
      template <class T, class F>
      struct then { typedef T type; };
  };

  template <>
  struct if_true<false>
  {
      template <class T, class F>
      struct then { typedef F type; };
  };
}}
#endif // SELECT_TYPE_DWA20010206_HPP
