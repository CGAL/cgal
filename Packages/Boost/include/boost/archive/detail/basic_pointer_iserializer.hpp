#ifndef BOOST_ARCHIVE_BASIC_ARCHIVE_POINTER_ISERIALIZER_HPP
#define BOOST_ARCHIVE_BASIC_ARCHIVE_POINTER_ISERIALIZER_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_pointer_oserializer.hpp: extenstion of type_info required for 
// serialization.

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/archive/detail/basic_serializer.hpp>

namespace boost {

// forward declarations
namespace archive {
namespace detail {

class basic_iarchive;
class basic_iserializer;

class basic_pointer_iserializer : public basic_serializer {
protected:
    explicit basic_pointer_iserializer(
        const boost::serialization::extended_type_info & type_
    ) :
        basic_serializer(type_)
    {}
    virtual ~basic_pointer_iserializer(){};
public:
    virtual const basic_iserializer & get_basic_serializer() const = 0;
    virtual void load_object_ptr(
        basic_iarchive & ar, 
        void * & x,
        const unsigned int file_version
    ) const = 0;
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_BASIC_ARCHIVE_POINTER_ISERIALIZER_HPP
