#ifndef BOOST_SERIALIZATION_EXTENDED_TYPE_INFO_TYPEID_HPP
#define BOOST_SERIALIZATION_EXTENDED_TYPE_INFO_TYPEID_HPP
/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

// extended_type_info_typeid.hpp: implementation for version that depends
// on runtime typing (rtti - typeid) but uses a user specified string
// as the portable class identifier.

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <typeinfo>

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

//#include <boost/static_warning.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_polymorphic.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <boost/serialization/extended_type_info.hpp>

namespace boost {
namespace serialization {

namespace detail {

class extended_type_info_typeid_0 : public extended_type_info
{
private:
    static const char * type_info_key;
    virtual bool
    less_than(const extended_type_info &rhs) const
    {
        return 0 != get_type().before(
            static_cast<const extended_type_info_typeid_0 &>(rhs).get_type()
        );
    }
    virtual bool
    equal_to(const extended_type_info &rhs) const
    {
        return 0 != get_type().operator==(
            static_cast<const extended_type_info_typeid_0 &>(rhs).get_type()
        );
    }
    virtual bool
    not_equal_to(const extended_type_info &rhs) const
    {
        return 0 != get_type().operator!=(
            static_cast<const extended_type_info_typeid_0 &>(rhs).get_type()
        );
    }
protected:
    extended_type_info_typeid_0() :
        extended_type_info(type_info_key)
    {}
public:
    virtual const std::type_info & get_type() const = 0;
};

// this derivation is used for creating search arguments
class extended_type_info_typeid_arg : public extended_type_info_typeid_0
{
private:
    const std::type_info & ti;
    virtual const std::type_info &get_type() const
    {
        return ti;
    }
public:
    extended_type_info_typeid_arg(const std::type_info & ti_)
        : ti(ti_)
    { 
        // note absense of self register and key as this is used only as
        // search argument given a type_info reference and is not to 
        // be added to the map.
    }
};

} // namespace detail

///////////////////////////////////////////////////////////////////////////////
template<class T>
class extended_type_info_typeid : public detail::extended_type_info_typeid_0
{
private:
    virtual const std::type_info & get_type() const {
        return typeid(T);
    }
    extended_type_info_typeid() :
        detail::extended_type_info_typeid_0()
    {
        self_register();    // add type to type table
    }
public:
    struct is_polymorphic
    {
        typedef BOOST_DEDUCED_TYPENAME boost::is_polymorphic<T>::type type;
        BOOST_STATIC_CONSTANT(bool, value = is_polymorphic::type::value);
    };
    static const extended_type_info *
    get_derived_extended_type_info(const T & t){
        // note: this implementation - based on usage of typeid (rtti)
        // only works if the class has at least one virtual function.
//      BOOST_STATIC_WARNING(
//          static_cast<bool>(is_polymorphic::value)
//      );
        detail::extended_type_info_typeid_arg etia(typeid(t));
        return extended_type_info::find(& etia);
    }
    static extended_type_info *
    get_instance(){
        static extended_type_info_typeid<T> instance;
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
    extended_type_info_typeid<const T>
/**/
#endif

#endif // BOOST_SERIALIZATION_EXTENDED_TYPE_INFO_TYPEID_HPP
