#ifndef BOOST_ARCHIVE_TEXT_OARCHIVE_HPP
#define BOOST_ARCHIVE_TEXT_OARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// text_oarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <ostream>

#include <cstddef> // std::size_t
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/archive/basic_text_oprimitive.hpp>
#include <boost/archive/basic_text_oarchive.hpp>

namespace boost { 
namespace archive {

template<class Archive>
class text_oarchive_impl : 
     /* protected ? */ public basic_text_oprimitive<std::ostream>,
     public basic_text_oarchive<Archive>
{
#ifdef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
public:
#else
    friend class detail::interface_oarchive<Archive>;
    friend class basic_text_oarchive<Archive>;
    friend class save_access;
protected:
#endif
    template<class T>
    void save(const T & t){
        this->newtoken();
        basic_text_oprimitive<std::ostream>::save(t);
    }
    void save(const char * t);
    #ifndef BOOST_NO_INTRINSIC_WCHAR_T
    void save(const wchar_t * t);
    #endif
    void save(const std::string &s);
    #ifndef BOOST_NO_STD_WSTRING
    void save(const std::wstring &ws);
    #endif
protected:
    text_oarchive_impl(std::ostream & os, unsigned int flags = 0) :
        basic_text_oprimitive<std::ostream>(
            os, 
            0 != (flags & no_codecvt)
        ),
        basic_text_oarchive<Archive>()
    {
        if(0 == (flags & no_header))
            basic_text_oarchive<Archive>::init();
    }
public:
    void save_binary(const void *address, std::size_t count){
        put('\n');
        this->end_preamble();
        #if ! defined(__MWERKS__)
        this->basic_text_oprimitive<std::ostream>::save_binary(
        #else
        this->basic_text_oprimitive::save_binary(
        #endif
            address, 
            count
        );
        this->delimiter = this->eol;
    }
};

// do not derive from this class.  If you want to extend this functionality
// via inhertance, derived from text_oarchive_impl instead.  This will
// preserve correct static polymorphism.
class text_oarchive : 
    public text_oarchive_impl<text_oarchive>
{
public:
    text_oarchive(std::ostream & os, unsigned int flags = 0) :
        text_oarchive_impl<text_oarchive>(os, flags)
    {
    }
};

} // namespace archive
} // namespace boost

// required by smart_cast for compilers not implementing 
// partial template specialization
BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(boost::archive::text_oarchive)

#endif // BOOST_ARCHIVE_TEXT_OARCHIVE_HPP
