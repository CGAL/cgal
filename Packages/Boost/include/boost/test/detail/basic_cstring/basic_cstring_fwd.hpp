//  (C) Copyright Gennadiy Rozental 2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : basic_cstring class wrap C string and provide std_string like 
//                interface
// ***************************************************************************

#ifndef BASIC_CSTRING_FWD_HPP_071894GER
#define BASIC_CSTRING_FWD_HPP_071894GER

#include <boost/detail/workaround.hpp>

namespace boost {

namespace unit_test {

template<class CharT> class         basic_cstring;
typedef basic_cstring<char const>   const_string;
#if BOOST_WORKAROUND(__DECCXX_VER, BOOST_TESTED_AT(60590041))
typedef const_string                literal_string;
#else
typedef const_string const          literal_string;
#endif

typedef char const* const           c_literal_string;

} // namespace unit_test

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:22  spion
//  Initial revision
//
//  Revision 1.5  2004/08/18 05:28:57  rogeeff
//  another tru64cxx65 workaround
//
//  Revision 1.4  2004/07/19 12:28:17  rogeeff
//  guard rename
//
//  Revision 1.3  2004/06/05 11:02:15  rogeeff
//  std::traits usage reworked
//
//  Revision 1.2  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.1  2004/05/11 11:00:55  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
// ***************************************************************************

#endif // BASIC_CSTRING_FWD_HPP_071894GER

