#ifndef BOOST_SERIALIZATION_BASE_OBJECT_HPP
#define BOOST_SERIALIZATION_BASE_OBJECT_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// base_object.hpp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/type_traits/is_pointer.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/static_assert.hpp>

#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/void_cast.hpp>

namespace boost {
    namespace archive{
        class basic_oarchive;
        class basic_iarchive;
    } // namespace archive
namespace serialization {

namespace detail {
    // only register void casts if the types are polymorphic
    template<class Base, class Derived>
    struct base_register
    {
        struct nothing {
            static void invoke(){}
        };
        template<class B, class D>
        struct reg{
            static void invoke(){
                void_cast_register<const D, const B>(
                    static_cast<const D *>(NULL),
                    static_cast<const B *>(NULL)
                );
            }
        };
        static void invoke(){
            typedef BOOST_DEDUCED_TYPENAME mpl::eval_if<
                BOOST_DEDUCED_TYPENAME type_info_implementation<Base>::type::is_polymorphic,
                mpl::identity<reg<Base, Derived> >,
                mpl::identity<nothing>
            >::type typex;
            typex::invoke();
        }
    };
    
    // get the base type for a given derived type
    // preserving the const-ness
    template<class Base, class Derived>
    struct base_cast
    {
        typedef BOOST_DEDUCED_TYPENAME
            mpl::if_<
                is_const<Derived>,
                const Base,
                Base
            >::type type;
        BOOST_STATIC_ASSERT(
            is_const<type>::value && is_const<Derived>::value
            || ! is_const<type>::value && ! is_const<Derived>::value
        );
    };
} // namespace detail

// BORLAND
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x551))
template<class Base, class Derived>
const Base & 
base_object(const Derived & d)
{
    BOOST_STATIC_ASSERT(! is_pointer<Derived>::value);
    detail::base_register<Base, Derived>::invoke();
    return static_cast<const Base &>(d);
}
#else
template<class Base, class Derived>
BOOST_DEDUCED_TYPENAME detail::base_cast<Base, Derived>::type & 
base_object(Derived &d)
{
    BOOST_STATIC_ASSERT(( is_base_and_derived<Base,Derived>::value));
    BOOST_STATIC_ASSERT(! is_pointer<Derived>::value);
    detail::base_register<Base, Derived>::invoke();
    typedef BOOST_DEDUCED_TYPENAME detail::base_cast<Base, Derived>::type type;
    return static_cast<type &>(d);
}
#endif

} // namespace serialization
} // namespace boost

#endif // BOOST_SERIALIZATION_BASE_OBJECT_HPP
