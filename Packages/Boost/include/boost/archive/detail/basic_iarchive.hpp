#ifndef BOOST_ARCHIVE_BASIC_IARCHIVE_HPP
#define BOOST_ARCHIVE_BASIC_IARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_iarchive.hpp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>
// can't use this - much as I'd like to as borland doesn't support it
// #include <boost/scoped_ptr.hpp>

#include <boost/serialization/tracking.hpp>

#include <boost/archive/basic_archive.hpp>

namespace boost {

namespace serialization {
    class extended_type_info;
} // namespace serialization

namespace archive {
namespace detail {

class basic_iarchive_impl;
class basic_iserializer;
class basic_pointer_iserializer;
//////////////////////////////////////////////////////////////////////
// class basic_iarchive - read serialized objects from a input stream
class basic_iarchive 
{
    friend class basic_iarchive_impl;
    // hide implementation of this class to minimize header conclusion
    // in client code. I couldn't used scoped pointer with borland
    // boost::scoped_ptr<basic_iarchive_impl> pimpl;
    basic_iarchive_impl * pimpl;

    virtual void vload(version_type &t) =  0;
    virtual void vload(object_id_type &t) =  0;
    virtual void vload(class_id_type &t) =  0;
    virtual void vload(class_id_optional_type &t) = 0;
    virtual void vload(class_name_type &t) = 0;
    virtual void vload(tracking_type &t) = 0;
    version_type archive_library_version;
public: // note: not part of the public API.
    void
    next_object_pointer(void *t);
protected:
    void init(unsigned int archive_library_version_){
        archive_library_version = archive_library_version_;
    }
    basic_iarchive();
    virtual ~basic_iarchive();
public:
    unsigned int library_version() const{
        return archive_library_version;
    }
    void load_object(
        void *t, 
        const basic_iserializer & bis
    );
    const basic_pointer_iserializer * load_pointer(
        void * & t, 
        const basic_pointer_iserializer * bpis_ptr,
        const basic_pointer_iserializer * (*finder)(
            const boost::serialization::extended_type_info & type
        )
    );
    void register_basic_serializer(const basic_iserializer & bis);
    void delete_created_pointers();
};

} // namespace detail
} // namespace archive
} // namespace boost

// required by smart_cast for compilers not implementing 
// partial template specialization
BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(
    boost::archive::detail::basic_iarchive
)

#endif //BOOST_ARCHIVE_BASIC_IARCHIVE_HPP
