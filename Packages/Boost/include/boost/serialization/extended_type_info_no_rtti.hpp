#ifndef BOOST_EXTENDED_TYPE_INFO_NO_RTTI_HPP
#define BOOST_EXTENDED_TYPE_INFO_NO_RTTI_HPP
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// extended_type_info_no_rtti.hpp: implementation for version that depends
// on runtime typing (rtti - typeid) but uses a user specified string
// as the portable class identifier.

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
#include <cassert>

#include <boost/config.hpp>
#include <cstring>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ using ::strcmp; }
#endif
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_const.hpp>

#include <boost/serialization/extended_type_info.hpp>
#include <boost/mpl/bool.hpp>

namespace boost {
namespace serialization {

///////////////////////////////////////////////////////////////////////
// define a special type_info that doesn't depend on rtti which is not
// available in all situations.

// common base class to share type_info_key.  This is used to 
// identify the method used to keep track of the extended type
class extended_type_info_no_rtti_base : public extended_type_info
{
protected:
    virtual bool
    less_than(const boost::serialization::extended_type_info &rhs) const {
        return std::strcmp(key, rhs.key) < 0;
    }
    virtual bool
    equal_to(const boost::serialization::extended_type_info &rhs) const{
        return std::strcmp(key, rhs.key) == 0;
    }
    virtual bool
    not_equal_to(const boost::serialization::extended_type_info &rhs) const {
        return std::strcmp(key, rhs.key) != 0;
    }
public:
    static const char * type_info_key;
    struct is_polymorphic
    {
        typedef boost::mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value = is_polymorphic::type::value);
    };
    extended_type_info_no_rtti_base() :
        boost::serialization::extended_type_info(type_info_key)
    {}
};

template<class T>
class extended_type_info_no_rtti : public extended_type_info_no_rtti_base
{
public:
    // Note: this version of extended_type_info
    // relies on the key used for exporting data.  
    // So we have to have the key when the instance is created and
    // can't wait for it be exported as in other cases.
    static const char * type_key;

    extended_type_info_no_rtti(){
        key_register(type_key);
        self_register();    // add type to type table
   }
    static const boost::serialization::extended_type_info *
    get_derived_extended_type_info(const T & t){
        // find the type that corresponds to the most derived type.
        // this implementation doesn't depend on typeid() but assumes
        // that the specified type has a function of the following signature.
        // A common implemention of such a function is to define as a virtual
        // function. 
        const char * derived_key = t.get_key();
        assert(NULL != derived_key);
        return boost::serialization::extended_type_info::find(derived_key);
    }

    static boost::serialization::extended_type_info *
    get_instance(){
        static extended_type_info_no_rtti<const T> instance;
        return & instance;
    }
};

} // namespace serialization
} // namespace boost

///////////////////////////////////////////////////////////////////////////////
// If no other implementation has been designated as default, 
// use this one.  To use this implementation as the default, specify it
// before any of the other headers.

#ifndef BOOST_SERIALIZATION_DEFAULT_TYPE_INFO
#define BOOST_SERIALIZATION_DEFAULT_TYPE_INFO(T) \
    extended_type_info_no_rtti<const T>
/**/
#endif

#endif // BOOST_EXTENDED_TYPE_INFO_NO_RTTI_HPP
