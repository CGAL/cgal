#ifndef BOOST_EXTENDED_TYPE_INFO_HPP
#define BOOST_EXTENDED_TYPE_INFO_HPP

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

#include <boost/config.hpp>
#include <boost/noncopyable.hpp>

#define BOOST_SERIALIZATION_MAX_KEY_SIZE 128

namespace boost { 
namespace serialization {

class extended_type_info : public boost::noncopyable 
{
private:
    friend bool
    operator<(const extended_type_info &lhs, const extended_type_info &rhs);
    friend bool
    operator==(const extended_type_info &lhs, const extended_type_info &rhs);
    friend bool
    operator!=(const extended_type_info &lhs, const extended_type_info &rhs);
    virtual bool
    less_than(const extended_type_info &rhs) const = 0;
    virtual bool
    equal_to(const extended_type_info &rhs) const = 0;
    virtual bool
    not_equal_to(const extended_type_info &rhs) const = 0;
protected:
    void self_register();
    extended_type_info(const char * type_info_key_) :
        type_info_key(type_info_key_),
        key(NULL)
    {}
    virtual ~extended_type_info(){};
    int type_info_key_cmp(const extended_type_info & rhs) const;
public:
    // used to uniquely identify the type of class derived from this one
    // so that different derivations of this class can be simultaneously
    // included in implementation of sets and maps.
    const char * type_info_key;
    // text string used as key for exporting serialized classes
    const char * key;
    static const extended_type_info * find(const char *key);
    static const extended_type_info * find(const extended_type_info * t);
    void key_register(const char *key);
};

bool
operator<(const extended_type_info &lhs, const extended_type_info &rhs);
bool
operator==(const extended_type_info &lhs, const extended_type_info &rhs);
bool
operator!=(const extended_type_info &lhs, const extended_type_info &rhs);

} // namespace serialization 
} // namespace boost

#endif // BOOST_EXTENDED_TYPE_INFO_HPP

