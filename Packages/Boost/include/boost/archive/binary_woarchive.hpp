#ifndef BOOST_ARCHIVE_BINARY_WOARCHIVE_HPP
#define BOOST_ARCHIVE_BINARY_WOARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// binary_woarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
#ifdef BOOST_NO_STD_WSTREAMBUF
#error "wide char i/o not supported on this platform"
#else

#include <ostream>
#include <boost/archive/detail/auto_link_warchive.hpp>
#include <boost/archive/basic_binary_oprimitive.hpp>
#include <boost/archive/basic_binary_oarchive.hpp>

#include <boost/archive/detail/abi_prefix.hpp> // must be the last header

namespace boost { 
namespace archive {

template<class Archive>
class binary_woarchive_impl : 
    public basic_binary_oprimitive<Archive, std::wostream>,
    public basic_binary_oarchive<Archive>
{
#ifdef BOOST_NO_MEMBER_TEMPLATE_FRIENDS
public:
#else
    friend class detail::interface_oarchive<Archive>;
    friend class basic_binary_oarchive<Archive>;
    friend class save_access;
protected:
#endif
    void init(){
        basic_binary_oarchive<Archive>::init();
        basic_binary_oprimitive<Archive, std::wostream>::init();
    }
    binary_woarchive_impl(std::wostream & os, unsigned int flags) :
        basic_binary_oprimitive<Archive, std::wostream>(
            os, 
            0 != (flags & no_codecvt)
        ),
        basic_binary_oarchive<Archive>(flags)
    {
       if(0 == (flags & no_header))
           init();
    }
};

// do not derive from this class.  If you want to extend this functionality
// via inhertance, derived from binary_oarchive_impl instead.  This will
// preserve correct static polymorphism.
class binary_woarchive : 
    public binary_woarchive_impl<binary_woarchive>
{
public:
    binary_woarchive(std::wostream & os, unsigned int flags = 0) :
        binary_woarchive_impl<binary_woarchive>(os, flags)
    {}
    ~binary_woarchive(){}
};

} // namespace archive
} // namespace boost

// required by smart_cast for compilers not implementing 
// partial template specialization
BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(boost::archive::binary_woarchive)

#include <boost/archive/detail/abi_suffix.hpp> // pops abi_suffix.hpp pragmas

#endif // BOOST_NO_STD_WSTREAMBUF
#endif // BOOST_ARCHIVE_BINARY_WOARCHIVE_HPP
