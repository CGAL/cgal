#ifndef BOOST_ARCHIVE_BINARY_IPRIMITIVE_HPP
#define BOOST_ARCHIVE_BINARY_IPRIMITIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_binary_iprimitive.hpp
//
// archives stored as native binary - this should be the fastest way
// to archive the state of a group of obects.  It makes no attempt to
// convert to any canonical form.

// IN GENERAL, ARCHIVES CREATED WITH THIS CLASS WILL NOT BE READABLE
// ON PLATFORM APART FROM THE ONE THEY ARE CREATED ON

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <iosfwd>
#include <cassert>

#include <cstddef> // std::size_t
#include <cstring>

#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::memcpy; 
    using ::strcpy;
    using ::size_t;
} // namespace std
#endif

#include <boost/throw_exception.hpp>
#include <boost/limits.hpp>
#include <boost/cstdint.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/scoped_ptr.hpp>

#include <boost/archive/archive_exception.hpp>
#include <boost/archive/codecvt_null.hpp>

namespace boost { 
namespace archive {

/////////////////////////////////////////////////////////////////////////////
// class binary_iarchive - read serialized objects from a input binary stream
template<class Archive, class IStream>
class basic_binary_iprimitive
{
#ifndef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
    friend class load_access;
protected:
#else
public:
#endif
    // return a pointer to the most derived class
    Archive * This(){
        return static_cast<Archive *>(this);
    }
    // native streams are always handled as bytes
    IStream &is;
    boost::scoped_ptr<std::locale> archive_locale;
//    boost::scoped_ptr<
//        codecvt_null<BOOST_DEDUCED_TYPENAME IStream::char_type> 
//    > archive_codecvt;
    io::basic_ios_locale_saver<
        BOOST_DEDUCED_TYPENAME IStream::char_type, BOOST_DEDUCED_TYPENAME IStream::traits_type
    > locale_saver;

    // main template for serilization of primitive types
    template<class T>
    void load(T & t){
        load_binary(& t, sizeof(T));
    }

    void load(char * t);
    void load(wchar_t * t);
    void load(std::string &s);
    #ifndef BOOST_NO_STD_WSTRING
    void load(std::wstring &ws);
    #endif

    void init();
    basic_binary_iprimitive(IStream  &is_, bool no_codecvt);
    ~basic_binary_iprimitive();
public:
    void load_binary(void *address, std::size_t count);
};

template<class Archive, class IStream>
inline void basic_binary_iprimitive<Archive, IStream>::load_binary(
    void *address, 
    std::size_t count
){
    assert(
        static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max()) >= count
    );
    if(is.fail())
        boost::throw_exception(archive_exception(archive_exception::stream_error));
    // note: an optimizer should eliminate the following for char files
    std::size_t s = count / sizeof(BOOST_DEDUCED_TYPENAME IStream::char_type);
    is.read(
        static_cast<BOOST_DEDUCED_TYPENAME IStream::char_type *>(address), 
        s
    );
    // note: an optimizer should eliminate the following for char files
    s = count % sizeof(BOOST_DEDUCED_TYPENAME IStream::char_type);
    if(0 < s){
        if(is.fail())
            boost::throw_exception(archive_exception(archive_exception::stream_error));
        BOOST_DEDUCED_TYPENAME IStream::char_type t;
        is.read(& t, 1);
        std::memcpy(address, &t, s);
    }
}

} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_BINARY_IPRIMITIVE_HPP
