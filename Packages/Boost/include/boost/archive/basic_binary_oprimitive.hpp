#ifndef BOOST_ARCHIVE_BASIC_BINARY_OPRIMITIVE_HPP
#define BOOST_ARCHIVE_BASIC_BINARY_OPRIMITIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_binary_oprimitive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

// archives stored as native binary - this should be the fastest way
// to archive the state of a group of obects.  It makes no attempt to
// convert to any canonical form.

// IN GENERAL, ARCHIVES CREATED WITH THIS CLASS WILL NOT BE READABLE
// ON PLATFORM APART FROM THE ONE THEY ARE CREATE ON

#include <iosfwd>
#include <cassert>
#include <locale>
#include <cstddef> // size_t

#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/cstdint.hpp>
#include <boost/limits.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/throw_exception.hpp>

#include <boost/archive/archive_exception.hpp>

namespace boost {
namespace archive {

/////////////////////////////////////////////////////////////////////////
// class basic_binary_oprimitive - binary output of prmitives

template<class Archive, class OStream>
class basic_binary_oprimitive
{
#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
    friend class save_access;
protected:
#else
public:
#endif
    // return a pointer to the most derived class
    Archive * This(){
        return static_cast<Archive *>(this);
    }
    // native binary streams are handled as bytes
    OStream &os;
    boost::scoped_ptr<std::locale> archive_locale;
    io::basic_ios_locale_saver<
        BOOST_DEDUCED_TYPENAME OStream::char_type, 
        BOOST_DEDUCED_TYPENAME OStream::traits_type
    > locale_saver;

    // default saving of primitives.
    template<class T>
    void save(const T & t)
    {
        save_binary(& t, sizeof(T));
    }

    void save(const char * t);
    void save(const wchar_t * t);
    void save(const std::string &s);
    #ifndef BOOST_NO_STD_WSTRING
    void save(const std::wstring &ws);
    #endif

    void init();
    basic_binary_oprimitive(OStream & os, bool no_codecvt);
    ~basic_binary_oprimitive();
public:
    void save_binary(const void *address, std::size_t count);
};

template<class Archive, class OStream>
inline void basic_binary_oprimitive<Archive, OStream>::save_binary(
    const void *address, 
    std::size_t count
){
    assert(
        static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()) >= count
    );
    // note: if the following assertions fail
    // a likely cause is that the output stream is set to "text"
    // mode where by cr characters recieve special treatment.
    // be sure that the output stream is opened with ios::binary
    if(os.fail())
        boost::throw_exception(archive_exception(archive_exception::stream_error));
    // figure number of elements to output - round up
    count = ( count + sizeof(BOOST_DEDUCED_TYPENAME OStream::char_type) - 1) 
        / sizeof(BOOST_DEDUCED_TYPENAME OStream::char_type);
    os.write(
        static_cast<const BOOST_DEDUCED_TYPENAME OStream::char_type *>(address), 
        count
    );
    assert(os.good());
}

} //namespace boost 
} //namespace archive 

#endif // BOOST_ARCHIVE_BASIC_BINARY_OPRIMITIVE_HPP
