/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_binary_iprimitive.ipp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <cassert>
#include <cstddef> // size_t

#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/detail/workaround.hpp> // fixup for RogueWave

#include <boost/throw_exception.hpp>
#include <boost/scoped_ptr.hpp>

#include <boost/archive/archive_exception.hpp>
#include <boost/archive/codecvt_null.hpp>
#include <boost/archive/add_facet.hpp>

namespace boost {
namespace archive {

//////////////////////////////////////////////////////////////////////
// implementation of basic_binary_iprimitive

template<class Archive, class IStream>
void basic_binary_iprimitive<Archive, IStream>::init()
{
    // Detect  attempts to pass native binary archives across
    // incompatible platforms. This is not fool proof but its
    // better than nothing.
    unsigned char size;
    this->This()->load(size);
    if(sizeof(int) != size)
        boost::throw_exception(
            archive_exception(archive_exception::incompatible_native_format)
        );
    this->This()->load(size);
    if(sizeof(long) != size)
        boost::throw_exception(
            archive_exception(archive_exception::incompatible_native_format)
        );
    this->This()->load(size);
    if(sizeof(float) != size)
        boost::throw_exception(
            archive_exception(archive_exception::incompatible_native_format)
        );
    this->This()->load(size);
    if(sizeof(double) != size)
        boost::throw_exception(
            archive_exception(archive_exception::incompatible_native_format)
        );

    // for checking endian
    int i;
    this->This()->load(i);
    if(1 != i)
        boost::throw_exception(
            archive_exception(archive_exception::incompatible_native_format)
        );
}

template<class Archive, class IStream>
void basic_binary_iprimitive<Archive, IStream>::load(wchar_t * ws)
{
    std::size_t l;
    this->This()->load(l);
    load_binary(ws, l);
    ws[l / sizeof(wchar_t)] = L'\0';
}

template<class Archive, class IStream>
void basic_binary_iprimitive<Archive, IStream>::load(std::string & s)
{
    std::size_t l;
    this->This()->load(l);
    // borland de-allocator fixup
    #if BOOST_WORKAROUND(_RWSTD_VER, BOOST_TESTED_AT(20101))
    if(NULL != s.data())
    #endif
        s.resize(l);
    // note breaking a rule here - could be a problem on some platform
    load_binary(const_cast<char *>(s.data()), l);
}

#ifndef BOOST_NO_CWCHAR
template<class Archive, class IStream>
void basic_binary_iprimitive<Archive, IStream>::load(char * s)
{
    std::size_t l;
    this->This()->load(l);
    load_binary(s, l);
    s[l] = '\0';
}
#endif

#ifndef BOOST_NO_STD_WSTRING
template<class Archive, class IStream>
void basic_binary_iprimitive<Archive, IStream>::load(std::wstring & ws)
{
    std::size_t l;
    this->This()->load(l);
    // borland de-allocator fixup
    #if BOOST_WORKAROUND(_RWSTD_VER, BOOST_TESTED_AT(20101))
    if(NULL != ws.data())
    #endif
        ws.resize(l);
    // note breaking a rule here - is could be a problem on some platform
    load_binary(const_cast<wchar_t *>(ws.data()), l * sizeof(wchar_t) / sizeof(char));
}
#endif

template<class Archive, class IStream>
basic_binary_iprimitive<Archive, IStream>::basic_binary_iprimitive(
    IStream &is_, 
    bool no_codecvt
) :
    is(is_),
    archive_locale(NULL),
    locale_saver(is)
{
    if(! no_codecvt){
        archive_locale.reset(
            boost::archive::add_facet(
                std::locale::classic(),
                new codecvt_null<BOOST_DEDUCED_TYPENAME IStream::char_type>
            )
        );
        is.imbue(* archive_locale);
    }
}

// scoped_ptr requires that archive_locale be a complete type at time of
// destruction so define destructor here rather than in the header
template<class Archive, class IStream>
basic_binary_iprimitive<Archive, IStream>::~basic_binary_iprimitive(){
}

} // namespace archive
} // namespace boost
