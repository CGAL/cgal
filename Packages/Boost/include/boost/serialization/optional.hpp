/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8

// (C) Copyright 2002-4 Pavel Vozenilek . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Provides non-intrusive serialization for boost::optional.

#ifndef BOOST_SERIALIZATION_OPTIONAL_HPP_
#define BOOST_SERIALIZATION_OPTIONAL_HPP_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#include <boost/config.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/assert.hpp>

// function specializations must be defined in the appropriate
// namespace - boost::serialization
namespace boost { 
#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
namespace serialization {
#endif
    
    template<class Archive, class T>
    void save(
        Archive & ar, 
        const boost::optional<T> & t, 
        const unsigned int /*version*/
    ){
        bool tflag = t;
        ar << boost::serialization::make_nvp("initialized", tflag);
        if (tflag)
           ar << boost::serialization::make_nvp("value", *t);
    }

    template<class Archive, class T>
    void load(
        Archive & ar, 
        boost::optional<T> & t, 
        const unsigned int /*version*/
    ){
        bool tflag;
        ar >> boost::serialization::make_nvp("initialized", tflag);
        if (tflag){
            T aux;
            ar >> boost::serialization::make_nvp("value", aux);
            t.reset(aux);
        }
        else {
            t.reset();
        }
    }

    template<class Archive, class T>
    void serialize(
        Archive & ar, 
        boost::optional<T> & t, 
        const unsigned int version
    ){
        boost::serialization::split_free(ar, t, version);
    }

    // the following would be slightly more efficient.  But it
    // would mean that archives created with programs that support
    // TPS wouldn't be readable by programs that don't support TPS.
    // Hence we decline to support this otherwise convenient optimization.
    //#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
    #if 0

    template <class T>
    struct implementation_level<optional<T> >
    {
        typedef mpl::integral_c_tag tag;
        typedef mpl::int_<boost::serialization::object_serializable> type;
        BOOST_STATIC_CONSTANT(
            int , 
            value = boost::serialization::implementation_level::type::value
        );
    };

    template<class T>
    struct tracking_level<optional<T> >
    {
        typedef mpl::integral_c_tag tag;
        typedef mpl::int_<boost::serialization::track_never> type;
        BOOST_STATIC_CONSTANT(
            int , 
            value = boost::serialization::tracking_level::type::value
        );
    };

    #endif

#ifdef BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP
} // serialization
#endif
} // namespace boost
#endif
