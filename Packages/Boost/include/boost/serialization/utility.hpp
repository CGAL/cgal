#ifndef  BOOST_SERIALIZATION_UTILITY_HPP
#define BOOST_SERIALIZATION_UTILITY_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// serialization/utility.hpp:
// serialization for stl utility templates

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <utility>
#include <boost/config.hpp>

#include <boost/serialization/nvp.hpp>

// function specializations must be defined in the appropriate
// namespace - boost::serialization
#if defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)
#define STD _STLP_STD
#else
#define STD std
#endif

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
namespace boost { namespace serialization {
#else
namespace STD {
#endif

// pair
template<class Archive, class F, class S>
inline void serialize(
    Archive & ar,
    STD::pair<F, S> & p,
    const unsigned int /* file_version */
){
    ar & boost::serialization::make_nvp("first", p.first);
    ar & boost::serialization::make_nvp("second", p.second);
}

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
}} // namespace boost::serialization
#else
} // namespace std
#endif
#undef STD

#endif // BOOST_SERIALIZATION_UTILITY_HPP
