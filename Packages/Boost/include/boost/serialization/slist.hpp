#ifndef BOOST_SERIALIZATION_SLIST_HPP
#define BOOST_SERIALIZATION_SLIST_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// slist.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
#ifdef BOOST_HAS_SLIST

#include <slist>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/nvp.hpp>

// function specializations must be defined in the appropriate
// namespace - boost::serialization
#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
    namespace boost { namespace serialization {
#else
    #if defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)
    #define STD _STLP_STD
    #else
    #define STD BOOST_STD_EXTENSION_NAMESPACE
    #endif
    namespace STD {
#endif

template<class Archive, class U, class Allocator>
inline void save(
    Archive & ar,
    const STD::slist<U, Allocator> &t,
    const unsigned int file_version
){
    boost::serialization::stl::save_collection<
        Archive,
        STD::slist<U, Allocator> 
    >(ar, t);
}

template<class Archive, class U, class Allocator>
inline void load(
    Archive & ar,
    STD::slist<U, Allocator> &t,
    const unsigned int file_version
){
    // retrieve number of elements
    t.clear();
    // retrieve number of elements
    unsigned int count;
    ar >> BOOST_SERIALIZATION_NVP(count);
    if(0 == count)
        return;

    boost::serialization::stl::stack_construct<Archive, U> u(ar);
    ar >> boost::serialization::make_nvp("item", u.reference());
    t.push_front(u.reference());
    BOOST_DEDUCED_TYPENAME BOOST_STD_EXTENSION_NAMESPACE::slist<U, Allocator>::iterator last;
    last = t.begin();
    while(--count > 0){
        boost::serialization::stl::stack_construct<Archive, U> u(ar);
        ar >> boost::serialization::make_nvp("item", u.reference());
        last = t.insert_after(last, u.reference());
    }
}

// split non-intrusive serialization function member into separate
// non intrusive save/load member functions
template<class Archive, class U, class Allocator>
inline void serialize(
    Archive & ar,
    STD::slist<U, Allocator> &t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
}} // namespace boost::serialization
#else
} // BOOST_STD_EXTENSION_NAMESPACE
#endif

#include <boost/serialization/collection_traits.hpp>

BOOST_SERIALIZATION_COLLECTION_TRAITS(STD::slist)
#undef STD

#endif  // BOOST_HAS_SLIST
#endif  // BOOST_SERIALIZATION_SLIST_HPP
