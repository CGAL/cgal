#ifndef BOOST_ARCHIVE_KNOWN_ARCHIVE_TYPES_HPP
#define BOOST_ARCHIVE_KNOWN_ARCHIVE_TYPES_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// known_archive_types.hpp: set traits of classes to be serialized

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

// list of archive type shipped with the serialization system

#include <boost/mpl/list.hpp>
#include <boost/mpl/pop_front.hpp>

#include <boost/archive/detail/known_archive_types_fwd.hpp>

namespace boost { 
namespace archive {
namespace detail {

class null_type;

template<>
struct known_archive_types<false> {
    typedef 
        mpl::pop_front<
            mpl::list<
                null_type
                #if defined(BOOST_ARCHIVE_TEXT_OARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::text_oarchive
                #endif
                #if defined(BOOST_ARCHIVE_TEXT_IARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::text_iarchive
                #endif
                #if defined(BOOST_ARCHIVE_TEXT_WOARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::text_woarchive
                #endif
                #if defined(BOOST_ARCHIVE_TEXT_WIARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::text_wiarchive
                #endif
                #if defined(BOOST_ARCHIVE_BINARY_OARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::binary_oarchive
                #endif
                #if defined(BOOST_ARCHIVE_BINARY_IARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::binary_iarchive
                #endif
                #if defined(BOOST_ARCHIVE_BINARY_WOARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::binary_woarchive
                #endif
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                #if defined(BOOST_ARCHIVE_BINARY_WIARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::binary_wiarchive
                #endif
                #if defined(BOOST_ARCHIVE_XML_OARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::xml_oarchive
                #endif
                #if defined(BOOST_ARCHIVE_XML_IARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::xml_iarchive
                #endif
                #if defined(BOOST_ARCHIVE_XML_WOARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::xml_woarchive
                #endif
                #if defined(BOOST_ARCHIVE_XML_WIARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::xml_wiarchive
                #endif
                #if defined(BOOST_ARCHIVE_POLYMORPHIC_OARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::polymorphic_oarchive
                #endif
                #if defined(BOOST_ARCHIVE_POLYMORPHIC_IARCHIVE_HPP)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , boost::archive::polymorphic_iarchive
                #endif
                #if defined(BOOST_ARCHIVE_CUSTOM_IARCHIVE_TYPES)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , BOOST_ARCHIVE_CUSTOM_IARCHIVE_TYPES
                #endif
                #if defined(BOOST_ARCHIVE_CUSTOM_OARCHIVE_TYPES)
                    #if ! defined(BOOST_ARCHIVE_EXPORT)
                        #define BOOST_ARCHIVE_EXPORT
                    #endif
                    , BOOST_ARCHIVE_CUSTOM_OARCHIVE_TYPES
                #endif
            >::type
        >::type type;
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_KNOWN_ARCHIVE_TYPES_HPP
