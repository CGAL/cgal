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
//  Description : 
// ***************************************************************************

#ifndef BOOST_IFSTREAM_LINE_ITERATOR_HPP_071894GER
#define BOOST_IFSTREAM_LINE_ITERATOR_HPP_071894GER

// Boost
#include <boost/test/detail/iterator/istream_line_iterator.hpp>

// STL
#include <fstream>

namespace boost {

namespace unit_test {

namespace ut_detail {

// ************************************************************************** //
// **************                ifstream_holder               ************** //
// ************************************************************************** //

template<typename CharT>
class ifstream_holder {
public:
    // Constructor
    explicit    ifstream_holder( basic_cstring<CharT const> file_name )
    {
        if( file_name.is_empty() )
            return;

        m_stream.open( file_name.begin(), std::ios::in );
    }

    bool is_valid()
    {
        return m_stream.is_open();
    }

protected:
#ifdef BOOST_CLASSIC_IOSTREAMS
    typedef std::ifstream                                       stream_t;
#else
    typedef std::basic_ifstream<CharT,std::char_traits<CharT> > stream_t;
#endif

    // Data members
    stream_t    m_stream;
};

} // namespace ut_detail

// ************************************************************************** //
// **************         basic_ifstream_line_iterator         ************** //
// ************************************************************************** //

template<typename CharT>
class basic_ifstream_line_iterator : ut_detail::ifstream_holder<CharT>, public basic_istream_line_iterator<CharT>
{
public:
    basic_ifstream_line_iterator( basic_cstring<CharT const> file_name, CharT delimeter )
    : ut_detail::ifstream_holder<CharT>( file_name ), basic_istream_line_iterator<CharT>( this->m_stream, delimeter ) {}

    explicit basic_ifstream_line_iterator( basic_cstring<CharT const> file_name = basic_cstring<CharT const>() )
    : ut_detail::ifstream_holder<CharT>( file_name ), basic_istream_line_iterator<CharT>( this->m_stream ) {}
};

typedef basic_ifstream_line_iterator<char>      ifstream_line_iterator;
typedef basic_ifstream_line_iterator<wchar_t>   wifstream_line_iterator;

} // namespace unit_test

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:23  spion
//  Initial revision
//
//  Revision 1.8  2004/09/19 09:22:13  rogeeff
//  ios fix for classic iostreams
//
//  Revision 1.7  2004/07/19 12:29:57  rogeeff
//  guard rename
//  mingw port
//
//  Revision 1.6  2004/06/29 04:32:25  rogeeff
//  no message
//
//  Revision 1.5  2004/06/07 07:33:49  rogeeff
//  detail namespace renamed
//
//  Revision 1.4  2004/06/05 11:03:12  rogeeff
//  input_iterator_adaptor simplified
//  token_iterator added
//
//  Revision 1.3  2004/05/27 07:01:49  rogeeff
//  portability workarounds
//
//  Revision 1.2  2004/05/25 10:29:09  rogeeff
//  use standard getline
//  eliminate initialize
//  proper handle \n in wide case
//
//  Revision 1.1  2004/05/21 06:30:10  rogeeff
//  ifstream_line_iterator added
//
// ***************************************************************************

#endif // BOOST_IFSTREAM_LINE_ITERATOR_HPP_071894GER

