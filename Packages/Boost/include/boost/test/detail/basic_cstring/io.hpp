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
//  Description : basic_cstring i/o implementation
// ***************************************************************************

#ifndef  BASIC_CSTRING_IO_HPP_071894GER
#define  BASIC_CSTRING_IO_HPP_071894GER

// Boost.Test
#include <boost/test/detail/basic_cstring/basic_cstring.hpp>

// STL
#include <iosfwd>
#include <string>

//____________________________________________________________________________//

namespace boost {

namespace unit_test {

#ifdef BOOST_CLASSIC_IOSTREAMS

template<typename CharT>
inline std::ostream&
operator<<( std::ostream& os, basic_cstring<CharT> const& str )
{
    typedef typename ut_detail::bcs_base_char<CharT>::type char_type;
    char_type const* const beg = reinterpret_cast<char_type const* const>( str.begin() );
    char_type const* const end = reinterpret_cast<char_type const* const>( str.end() );
    os << std::basic_string<char_type>( beg, end - beg );

    return os;
}

#else

template<typename CharT1, typename Tr,typename CharT2>
inline std::basic_ostream<CharT1,Tr>&
operator<<( std::basic_ostream<CharT1,Tr>& os, basic_cstring<CharT2> const& str )
{
    CharT1 const* const beg = reinterpret_cast<CharT1 const*>( str.begin() ); //!!
    CharT1 const* const end = reinterpret_cast<CharT1 const*>( str.end() );
    os << std::basic_string<CharT1,Tr>( beg, end - beg );

    return os;
}

#endif

//____________________________________________________________________________//


} // namespace unit_test

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:23  spion
//  Initial revision
//
//  Revision 1.9  2004/09/19 09:22:13  rogeeff
//  ios fix for classic iostreams
//
//  Revision 1.8  2004/07/21 16:28:09  dgregor
//  io.hpp: Eliminate useless "const" qualifiers on an rvalue
//
//  Revision 1.7  2004/07/19 12:27:05  rogeeff
//  guard rename
//
//  Revision 1.6  2004/06/30 07:52:56  rogeeff
//  typo fix
//
//  Revision 1.5  2004/06/29 04:31:49  rogeeff
//  gcc 2.95 fix
//
//  Revision 1.4  2004/06/05 11:02:15  rogeeff
//  std::traits usage reworked
//
//  Revision 1.3  2004/05/27 06:24:44  rogeeff
//  workaround for gcc 2.95 io
//
//  Revision 1.2  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.1  2004/05/11 11:00:55  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
// ***************************************************************************

#endif // BASIC_CSTRING_IO_HPP_071894GER
