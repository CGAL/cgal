#ifndef BOOST_ARCHIVE_TEXT_IARCHIVE_HPP
#define BOOST_ARCHIVE_TEXT_IARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// text_iarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <istream>

#include <boost/config.hpp>
#include <boost/archive/basic_text_iarchive.hpp>
#include <boost/archive/basic_text_iprimitive.hpp>

namespace boost { 
namespace archive {

template<class Archive>
class text_iarchive_impl : 
    public basic_text_iprimitive<std::istream>,
    public basic_text_iarchive<Archive>
{
protected:
#ifdef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
public:
#else
    friend class detail::interface_iarchive<Archive>;
    friend class basic_text_iarchive<Archive>;
    friend class load_access;
protected:
#endif
    template<class T>
    void load(T & t){
        basic_text_iprimitive<std::istream>::load(t);
    }
    void load(char * t);
    #ifndef BOOST_NO_INTRINSIC_WCHAR_T
    void load(wchar_t * t);
    #endif
    void load(std::string &s);
    #ifndef BOOST_NO_STD_WSTRING
    void load(std::wstring &ws);
    #endif
    // note: the following should not needed - but one compiler (vc 7.1)
    // fails to compile one test (test_shared_ptr) without it !!!
    // make this protected so it can be called from a derived archive
    template<class T>
    void load_override(T & t, BOOST_PFTO int){
        basic_text_iarchive<Archive>::load_override(t, 0);
    }
    text_iarchive_impl(std::istream & is, unsigned int flags = 0) :
        basic_text_iprimitive<std::istream>(
            is, 
            0 != (flags & no_codecvt)
        ),
        basic_text_iarchive<Archive>()
    {
        if(0 == (flags & no_header))
            basic_text_iarchive<Archive>::init();
    }
};

// do not derive from this class.  If you want to extend this functionality
// via inhertance, derived from text_iarchive_impl instead.  This will
// preserve correct static polymorphism.
class text_iarchive : 
    public text_iarchive_impl<text_iarchive>
{
public:
    text_iarchive(std::istream & is, unsigned int flags = 0) :
        text_iarchive_impl<text_iarchive>(is, flags)
    {
    }
};

} // namespace archive
} // namespace boost

// required by smart_cast for compilers not implementing 
// partial template specialization
BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(boost::archive::text_iarchive)

#endif // BOOST_ARCHIVE_TEXT_IARCHIVE_HPP
