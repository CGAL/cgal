#ifndef BOOST_SERIALIZATION_TRACKING_HPP
#define BOOST_SERIALIZATION_TRACKING_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// tracking.hpp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/greater.hpp>

#include <boost/type_traits/is_base_and_derived.hpp>
#include <boost/serialization/traits.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/tracking_enum.hpp>

namespace boost {
namespace serialization {

// default tracking level
template<class T>
struct tracking_level {
    template<class U>
    struct traits_class_tracking {
        typedef BOOST_DEDUCED_TYPENAME U::tracking type;
    };
    
    typedef mpl::integral_c_tag tag;
    typedef
        BOOST_DEDUCED_TYPENAME mpl::eval_if<
            is_base_and_derived<basic_traits, T>,
            traits_class_tracking<T>,
        //else
        BOOST_DEDUCED_TYPENAME mpl::eval_if<
            // for primitives
            BOOST_DEDUCED_TYPENAME mpl::equal_to<
                implementation_level<T>,
                mpl::int_<primitive_type> 
            >,
            // is never
            mpl::int_<track_never>,
            // otherwise its selective
            mpl::int_<track_selectivly>
            >
        >::type type;
    BOOST_STATIC_CONSTANT(int, value = tracking_level::type::value);
};


template<class T, enum tracking_type L>
inline bool operator>=(tracking_level<T> t, enum tracking_type l)
{
    return t.value >= (int)l;
}

} // namespace serialization
} // namespace boost

// specify the current tracking behavior for the class
#define BOOST_CLASS_TRACKING(T, E)           \
namespace boost {                            \
namespace serialization {                    \
template<>                                   \
struct tracking_level<T >                    \
{                                            \
    typedef mpl::integral_c_tag tag;         \
    typedef mpl::int_< E> type;              \
    BOOST_STATIC_CONSTANT(                   \
        int,                                 \
        value = tracking_level::type::value  \
    );                                       \
    /* tracking for a class  */              \
    BOOST_STATIC_ASSERT((                    \
        mpl::greater<                        \
            /* that is a prmitive */         \
            implementation_level<T >,        \
            mpl::int_<primitive_type>        \
        >::value                             \
    ));                                      \
};                                           \
}}

#endif // BOOST_SERIALIZATION_TRACKING_HPP
