//  (C) Copyright Gennadiy Rozental 2002-2003.
//  Use, modification, and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : wraps strstream and stringstream (depends with one is present )
//                to prodive the unified interface
// ***************************************************************************

#ifndef BOOST_WRAP_STRINGSTREAM_HPP
#define BOOST_WRAP_STRINGSTREAM_HPP

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
// **************               wrap_stringstream              ************** //
// ************************************************************************** //

class wrap_stringstream {
#ifdef BOOST_NO_STRINGSTREAM
    typedef std::ostrstream     wrapped_stream;
#else
    typedef std::ostringstream  wrapped_stream;
#endif // BOOST_NO_STRINGSTREAM
public:

    // access methods
    wrap_stringstream&          ref();
    wrapped_stream&             stream();
    std::string const&          str();

private:
    // Data members
    wrapped_stream              m_stream;
    std::string                 m_str;
};

//____________________________________________________________________________//

template <class T>
inline wrap_stringstream&
operator<<( wrap_stringstream& targ, T const& t )
{
    targ.stream() << t;
    return targ;
}

//____________________________________________________________________________//

inline wrap_stringstream::wrapped_stream&
wrap_stringstream::stream()
{
    return m_stream;
}

//____________________________________________________________________________//

inline wrap_stringstream&
wrap_stringstream::ref()
{ 
    return *this;
}

//____________________________________________________________________________//

inline std::string const&
wrap_stringstream::str()
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

inline wrap_stringstream&
operator<<( wrap_stringstream& targ, wrap_stringstream& src )
{
    targ << src.str();
    return targ;
}

#ifndef BOOST_NO_STD_LOCALE
//____________________________________________________________________________//

inline wrap_stringstream&
operator<<( wrap_stringstream& targ, std::ios_base& (*man)(std::ios_base&) )
{
    targ.stream() << man;
    return targ;
}

//____________________________________________________________________________//

template<typename Elem,typename Tr>
inline wrap_stringstream&
operator<<( wrap_stringstream& targ, std::basic_ostream<Elem,Tr>& (*man)(std::basic_ostream<Elem, Tr>&) )
{
    targ.stream() << man;
    return targ;
}

//____________________________________________________________________________//

template<typename Elem,typename Tr>
inline wrap_stringstream&
operator<<( wrap_stringstream& targ, std::basic_ios<Elem, Tr>& (*man)(std::basic_ios<Elem, Tr>&) )
{
    targ.stream() << man;
    return targ;
}

#endif
//____________________________________________________________________________//

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
//  Revision 1.1  2004/05/23 10:51:39  spion
//  Initial revision
//
//  Revision 1.9  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif  // BOOST_WRAP_STRINGSTREAM_HPP
