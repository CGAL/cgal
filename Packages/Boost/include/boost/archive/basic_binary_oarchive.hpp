#ifndef BOOST_ARCHIVE_BASIC_BINARY_OARCHIVE_HPP
#define BOOST_ARCHIVE_BASIC_BINARY_OARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_binary_oarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

// archives stored as native binary - this should be the fastest way
// to archive the state of a group of obects.  It makes no attempt to
// convert to any canonical form.

// IN GENERAL, ARCHIVES CREATED WITH THIS CLASS WILL NOT BE READABLE
// ON PLATFORM APART FROM THE ONE THEY ARE CREATE ON

#include <boost/config.hpp>
#include <boost/pfto.hpp>

#include <boost/detail/workaround.hpp>
#include <boost/archive/detail/interface_oarchive.hpp>
#include <boost/archive/detail/common_oarchive.hpp>

#include <boost/serialization/string.hpp>

namespace boost {
namespace archive {

//////////////////////////////////////////////////////////////////////
// class basic_binary_oarchive - write serialized objects to a binary output stream
// note: this archive has no pretensions to portability.  Archive format
// may vary across machine architectures and compilers.  About the only
// guarentee is that an archive created with this code will be readable
// by a program built with the same tools for the same machne.  This class
// does have the virtue of buiding the smalles archive in the minimum amount
// of time.  So under some circumstances it may be he right choice.

/////////////////////////////////////////////////////////////////////////
// class basic_text_iarchive - read serialized objects from a input text stream
template<class Archive>
class basic_binary_oarchive : public detail::common_oarchive<Archive>
{
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
public:
#elif defined(BOOST_MSVC)
    // for some inexplicable reason insertion of "class" generates compile erro
    // on msvc 7.1
    friend detail::interface_oarchive<Archive>;
protected:
#else
    friend class detail::interface_oarchive<Archive>;
protected:
#endif
    // any datatype not specifed below will be handled
    // by this function
    template<class T>
    void save_override(const T & t, BOOST_PFTO int)
    {
        archive::save(* this->This(), t);
    }
    // binary files don't include the optional information 
    void save_override(const class_id_optional_type & /* t */, int){}

    void save_override(const version_type & t, int){
        // upto 255 versions
        // note:t.t resolves borland ambguity
        unsigned  char x = t.t;
        * this->This() << x;
    }
    void save_override(const class_id_type & t, int){
        // upto 32K classes
        int_least16_t x = t.t;
        * this->This() << x;
    }
    void save_override(const class_id_reference_type & t, int){
        // upto 32K classes
        int_least16_t x = t.t;
        * this->This() << x;
    }
    void save_override(const object_id_type & t, int){
        // upto 2G objects
        uint_least32_t x = t.t;
        * this->This() << x;
    }
    void save_override(const object_reference_type & t, int){
        // upto 2G objects
        uint_least32_t x = t.t;
        * this->This() << x;
    }
    void save_override(const tracking_type & t, int){
        char x = t.t;
        * this->This() << x;
    }

    // explicitly convert to char * to avoid compile ambiguities
    void save_override(const class_name_type & t, int){
        * this->This() << std::string(static_cast<const char *>(t));
    }

    basic_binary_oarchive(unsigned int flags = 0){}
};

} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_BASIC_BINARY_OARCHIVE_HPP
