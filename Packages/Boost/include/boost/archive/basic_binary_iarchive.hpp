#ifndef BOOST_ARCHIVE_BASIC_BINARY_IARCHIVE_HPP
#define BOOST_ARCHIVE_BASIC_BINARY_IARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_binary_iarchive.hpp
//
// archives stored as native binary - this should be the fastest way
// to archive the state of a group of obects.  It makes no attempt to
// convert to any canonical form.

// IN GENERAL, ARCHIVES CREATED WITH THIS CLASS WILL NOT BE READABLE
// ON PLATFORM APART FROM THE ONE THEY ARE CREATED ON

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
#include <cstring>

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
    namespace std{ using ::memcpy; }
#endif

#include <boost/throw_exception.hpp>
#include <boost/pfto.hpp>

#include <boost/archive/detail/interface_iarchive.hpp>
#include <boost/archive/detail/common_iarchive.hpp>

#include <boost/serialization/string.hpp>

namespace boost { 
namespace archive {

/////////////////////////////////////////////////////////////////////////
// class basic_binary_iarchive - read serialized objects from a input binary stream
template<class Archive>
class basic_binary_iarchive : public detail::common_iarchive<Archive>
{
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1300)
public:
#elif defined(BOOST_MSVC)
    // for some inexplicable reason insertion of "class" generates compile erro
    // on msvc 7.1
    friend detail::interface_iarchive<Archive>;
protected:
#else
    friend class detail::interface_iarchive<Archive>;
protected:
#endif
    // intermediate level to support override of operators
    // fot templates in the absence of partial function 
    // template ordering
    template<class T>
    void load_override(T & t, BOOST_PFTO int)
    {
        archive::load(* this->This(), t);
    }
    // binary files don't include the optional information 
    void load_override(class_id_optional_type & /* t */, int){}

    // the following have been overridden to provide specific sizes
    // for these pseudo prmitive types.
    void load_override(version_type & t, int){ 
        // upto 255 versions
        unsigned char x;
        * this->This() >> x;
        t = version_type(x);
    }
    void load_override(class_id_type & t, int){
        // upto 32K classes
        int_least16_t x;
        * this->This() >> x;
        t = class_id_type(x);
    }
    void load_override(class_id_reference_type & t, int){
        // upto 32K classes
        int_least16_t x;
        * this->This() >> x;
        t = class_id_reference_type(x);
    }
    void load_override(object_id_type & t, int){
        // upto 2G objects
        uint_least32_t x;
        * this->This() >> x;
        t = object_id_type(x);
    }
    void load_override(object_reference_type & t, int){
        // upto 2G objects
        uint_least32_t x;
        * this->This() >> x;
        t = object_reference_type(x);
    }
    void load_override(tracking_type & t, int){
        char x;
        * this->This() >> x;
        t = (0 != x);
    }

    void load_override(class_name_type & t, int){
        std::string cn;
        cn.reserve(BOOST_SERIALIZATION_MAX_KEY_SIZE);
        load_override(cn, 0);
        if(cn.size() > (BOOST_SERIALIZATION_MAX_KEY_SIZE - 1))
            boost::throw_exception(
                archive_exception(archive_exception::invalid_class_name)
           );
        std::memcpy(t, cn.data(), cn.size());
        // .t is a borland tweak
        t.t[cn.size()] = '\0';
    }

    basic_binary_iarchive() {}
};

} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_BASIC_BINARY_IARCHIVE_HPP
