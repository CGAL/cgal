#ifndef BOOST_SERIALIZATION_EXTENDED_TYPE_INFO_HPP
#define BOOST_SERIALIZATION_EXTENDED_TYPE_INFO_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// extended_type_info.hpp: interface for portable version of type_info

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

// for now, extended type info is part of the serialization libraries
// this could change in the future.
#include <boost/noncopyable.hpp>
#include <boost/serialization/config.hpp>

#include <boost/config/abi_prefix.hpp> // must be the last header
#ifdef BOOST_MSVC
#  pragma warning(push)
#  pragma warning(disable : 4251 4231 4660 4275)
#endif

#define BOOST_SERIALIZATION_MAX_KEY_SIZE 128

namespace boost { 
namespace serialization {

class BOOST_SERIALIZATION_DECL(BOOST_PP_EMPTY()) extended_type_info : 
    private boost::noncopyable 
{
private:
    virtual bool
    less_than(const extended_type_info &rhs) const = 0;
    // used to uniquely identify the type of class derived from this one
    // so that different derivations of this class can be simultaneously
    // included in implementation of sets and maps.
    const char * type_info_key;
    int type_info_key_cmp(const extended_type_info & rhs) const;
protected:
    const char * key;
    extended_type_info(const char * type_info_key_);
    virtual ~extended_type_info() = 0;
public:
    void self_register();
    void key_register(const char *key);
    bool operator<(const extended_type_info &rhs) const;
    bool operator==(const extended_type_info &rhs) const {
        return this == & rhs;
    }
    bool operator!=(const extended_type_info &rhs) const {
        return this != & rhs;
    }
    const char * get_key() const {
        return key;
    }
    static const extended_type_info * find(const char *key);
    static const extended_type_info * find(const extended_type_info * t);
};

} // namespace serialization 
} // namespace boost

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
#include <boost/config/abi_suffix.hpp> // pops abi_suffix.hpp pragmas

#endif // BOOST_SERIALIZATION_EXTENDED_TYPE_INFO_HPP

