#ifndef BOOST_ARCHIVE_DETAIL_COMMON_IARCHIVE_HPP
#define BOOST_ARCHIVE_DETAIL_COMMON_IARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// common_iarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/throw_exception.hpp>

#include <boost/config.hpp>

#include <boost/archive/archive_exception.hpp>
#include <boost/archive/detail/basic_iarchive.hpp>
#include <boost/archive/detail/interface_iarchive.hpp>

namespace boost {
namespace archive {
namespace detail {

// note: referred to as Curiously Recurring Template Patter (CRTP)
template<class Archive>
class common_iarchive : 
    public basic_iarchive,
    public interface_iarchive<Archive>
{
private:
    virtual void vload(version_type & t){
        * this->This() >> t; 
    }
    virtual void vload(object_id_type & t){
        * this->This() >> t;
    }
    virtual void vload(class_id_type & t){
        * this->This() >> t;
    }
    virtual void vload(class_id_optional_type & t){
        * this->This() >> t;
    }
    virtual void vload(tracking_type & t){
        * this->This() >> t;
    }
    virtual void vload(class_name_type &s){
        * this->This() >> s;
    }
protected:
    // default implementations of functions which emit start/end tags for
    // archive types that require them.
    void load_start(const char *name){}
    void load_end(const char *name){}
    // default archive initialization
    void init(){
        // read signature in an archive version independent manner
        std::string file_signature;
        * this->This() >> file_signature;
        if(file_signature != ARCHIVE_SIGNATURE)
            boost::throw_exception(
                archive_exception(archive_exception::invalid_signature)
           );

        // make sure the version of the reading archive library can
        // support the format of the archive being read
        version_type input_library_version;
        * this->This() >> input_library_version;
        basic_iarchive::init(input_library_version);
        const unsigned int current_library_version = ARCHIVE_VERSION;
        // extra little .t is to get around borland quirk
        if(current_library_version < input_library_version.t)
            boost::throw_exception(
                archive_exception(archive_exception::unsupported_version)
            );
    }
    common_iarchive() : 
        basic_iarchive(),
        interface_iarchive<Archive>()
    {}
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_DETAIL_COMMON_IARCHIVE_HPP

