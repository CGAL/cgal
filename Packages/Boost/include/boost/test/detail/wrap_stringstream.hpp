//  (C) Copyright Gennadiy Rozental 2002-2003.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : wraps strstream and stringstream (depends with one is present )
//                to prodive the unified interface
// ***************************************************************************

#ifndef BOOST_WRAP_STRINGSTREAM_HPP_071894GER
#define BOOST_WRAP_STRINGSTREAM_HPP_071894GER

// STL
#ifdef BOOST_NO_STRINGSTREAM
#include <strstream>        // for std::ostrstream
#else
#include <sstream>          // for std::ostringstream
#endif // BOOST_NO_STRINGSTREAM

#ifdef BOOST_MSVC
# pragma warning(push)
# pragma warning(disable: 4511) // copy constructor could not be generated
# pragma warning(disable: 4512) // assignment operator could not be generated
#endif

namespace boost {

// ************************************************************************** //
// **************            basic_wrap_stringstream           ************** //
// ************************************************************************** //

template<typename CharT>
class basic_wrap_stringstream {
#ifdef BOOST_CLASSIC_IOSTREAMS
    typedef std::ostringstream               wrapped_stream;
#else 
#ifdef BOOST_NO_STRINGSTREAM
    typedef std::basic_ostrstream<CharT>     wrapped_stream;
#else
    typedef std::basic_ostringstream<CharT>  wrapped_stream;
#endif // BOOST_NO_STRINGSTREAM

#endif // WORKAROUND

public:
    // Access methods
    basic_wrap_stringstream&        ref();
    wrapped_stream&                 stream();
    std::basic_string<CharT> const& str();

private:
    // Data members
    wrapped_stream                  m_stream;
    std::basic_string<CharT>        m_str;
};

//____________________________________________________________________________//

template <typename CharT, typename T>
inline basic_wrap_stringstream<CharT>&
operator<<( basic_wrap_stringstream<CharT>& targ, T const& t )
{
    targ.stream() << t;
    return targ;
}

//____________________________________________________________________________//

template <typename CharT>
inline typename basic_wrap_stringstream<CharT>::wrapped_stream&
basic_wrap_stringstream<CharT>::stream()
{
    return m_stream;
}

//____________________________________________________________________________//

template <typename CharT>
inline basic_wrap_stringstream<CharT>&
basic_wrap_stringstream<CharT>::ref()
{ 
    return *this;
}

//____________________________________________________________________________//

template <typename CharT>
inline std::basic_string<CharT> const&
basic_wrap_stringstream<CharT>::str()
{

#ifdef BOOST_NO_STRINGSTREAM
    m_str.assign( m_stream.str(), m_stream.pcount() );
    m_stream.freeze( false );
#else
    m_str = m_stream.str();
#endif

    return m_str;
}

//____________________________________________________________________________//

template <typename CharT>
inline basic_wrap_stringstream<CharT>&
operator<<( basic_wrap_stringstream<CharT>& targ, basic_wrap_stringstream<CharT>& src )
{
    targ << src.str();
    return targ;
}

//____________________________________________________________________________//

#if !defined(BOOST_NO_STD_LOCALE) && BOOST_WORKAROUND(BOOST_MSVC, >= 1310)

template <typename CharT>
inline basic_wrap_stringstream<CharT>&
operator<<( basic_wrap_stringstream<CharT>& targ, std::ios_base& (*man)(std::ios_base&) )
{
    targ.stream() << man;
    return targ;
}

//____________________________________________________________________________//

template<typename CharT,typename Elem,typename Tr>
inline basic_wrap_stringstream<CharT>&
operator<<( basic_wrap_stringstream<CharT>& targ, std::basic_ostream<Elem,Tr>& (*man)(std::basic_ostream<Elem, Tr>&) )
{
    targ.stream() << man;
    return targ;
}

//____________________________________________________________________________//

template<typename CharT,typename Elem,typename Tr>
inline basic_wrap_stringstream<CharT>&
operator<<( basic_wrap_stringstream<CharT>& targ, std::basic_ios<Elem, Tr>& (*man)(std::basic_ios<Elem, Tr>&) )
{
    targ.stream() << man;
    return targ;
}

//____________________________________________________________________________//

#endif

// ************************************************************************** //
// **************               wrap_stringstream              ************** //
// ************************************************************************** //

typedef basic_wrap_stringstream<char>       wrap_stringstream;
typedef basic_wrap_stringstream<wchar_t>    wrap_wstringstream;

}  // namespace boost

#ifdef BOOST_MSVC
# pragma warning(default: 4511) // copy constructor could not be generated
# pragma warning(default: 4512) // assignment operator could not be generated
# pragma warning(pop)
#endif

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:21  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.14  2004/09/19 09:22:12  rogeeff
//  ios fix for classic iostreams
//
//  Revision 1.13  2004/07/19 12:24:32  rogeeff
//  guard rename
//
//  Revision 1.12  2004/05/27 06:23:22  rogeeff
//  workaround for gcc 2.95 io
//  workaround for msvc < 7.1 for manipulator usage
//
//  Revision 1.11  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.10  2004/05/11 11:00:53  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.9  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif  // BOOST_WRAP_STRINGSTREAM_HPP_071894GER
