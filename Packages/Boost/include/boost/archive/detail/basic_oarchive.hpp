#ifndef BOOST_ARCHIVE_BASIC_OARCHIVE_HPP
#define BOOST_ARCHIVE_BASIC_OARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_oarchive.hpp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/config.hpp>

#include <boost/detail/workaround.hpp>

// can't use this - much as I'd like to as borland doesn't support it
// #include <boost/scoped_ptr.hpp>

#include <boost/archive/basic_archive.hpp>
#include <boost/serialization/tracking.hpp>

namespace boost {
namespace archive {
namespace detail {

class basic_oarchive_impl;
class basic_oserializer;
class basic_pointer_oserializer;
//////////////////////////////////////////////////////////////////////
// class basic_oarchive - write serialized objects to an output stream
class basic_oarchive
{
    friend class basic_oarchive_impl;
    // hide implementation of this class to minimize header conclusion
    // in client code. note: borland can't use scoped_ptr
    //boost::scoped_ptr<basic_oarchive_impl> pimpl;
    basic_oarchive_impl * pimpl;

    // overload these to bracket object attributes. Used to implement
    // xml archives
    virtual void vsave(const version_type t) =  0;
    virtual void vsave(const object_id_type t) =  0;
    virtual void vsave(const object_reference_type t) =  0;
    virtual void vsave(const class_id_type t) =  0;
    virtual void vsave(const class_id_optional_type t) = 0;
    virtual void vsave(const class_id_reference_type t) =  0;
    virtual void vsave(const class_name_type & t) = 0;
    virtual void vsave(const tracking_type t) = 0;

protected:
    basic_oarchive();
    virtual ~basic_oarchive();

public:
    unsigned int library_version() const{
        return ARCHIVE_VERSION;
    }
    void save_object(
        const void *x, 
        const basic_oserializer & bos
    );
    void save_pointer(
        const void * t, 
        const basic_pointer_oserializer * bpos_ptr
    );
    void register_basic_serializer(const basic_oserializer & bos);
    void save_null_pointer(){
        vsave(null_pointer_tag);
    }
    void end_preamble(){} // default implementation does nothing
};

} // namespace detail
} // namespace archive
} // namespace boost


// required by smart_cast for compilers not implementing 
// partial template specialization
BOOST_BROKEN_COMPILER_TYPE_TRAITS_SPECIALIZATION(
    boost::archive::detail::basic_oarchive
)

#endif //BOOST_ARCHIVE_BASIC_OARCHIVE_HPP
