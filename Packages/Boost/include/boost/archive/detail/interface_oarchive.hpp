#ifndef BOOST_ARCHIVE_DETAIL_INTERFACE_OARCHIVE_HPP
#define BOOST_ARCHIVE_DETAIL_INTERFACE_OARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// interface_oarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
#include <string>
#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/cstdint.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/static_warning.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/level.hpp>
#include <boost/archive/detail/oserializer.hpp>

namespace boost { 
namespace archive {
namespace detail {

class basic_oserializer;
class basic_pointer_oserializer;

template<class Archive>
class interface_oarchive 
{
protected:
    interface_oarchive(){};
public:
    /////////////////////////////////////////////////////////
    // archive public interface

    struct is_loading {
        typedef mpl::bool_<false> type;
        BOOST_STATIC_CONSTANT(bool, value=false);
    };
    struct is_saving {
        typedef mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value=true);
    };

    // return a pointer to the most derived class
    Archive * This(){
        return static_cast<Archive *>(this);
    }

    template<class T>
    const basic_pointer_oserializer * register_type(T * t = NULL){
        const basic_pointer_oserializer & bpos =
            instantiate_pointer_oserializer(
                static_cast<Archive *>(NULL),
                static_cast<T *>(NULL)
            );
        this->This()->register_basic_serializer(bpos.get_basic_serializer());
        return & bpos;
    }

    // default processing - invoke serialization library
    template<class T>
    void save_override(T & t, /*BOOST_PFTO*/ int){
        archive::save(* this->This(), t);
    }

    // note: we presume that older compilers will never create a const
    // argument from a non-const by copyiing
    template<class T>
    Archive & operator<<(const T & t){
        this->This()->save_override(t, 0);
        return * this->This();
    }

    // the & operator 
    template<class T>
    Archive & operator&(const T & t){
        this->This()->save_override(t, 0);
        return * this->This();
    }

    // define operators for non-const arguments.  Don't depend one the const
    // ones below because the compiler MAY make a temporary copy to
    // create the const parameter (Though I havn't seen this happen). 
    #ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
        // the << operator
        template<class T>
        Archive & operator<<(T & t){
            // if trap here, we're saving a tracted non-const
            // value - this could be a stack variable with the same
            // address for multiple items. This would be the source of very
            // subtle errors and should be double checked
            // BOOST_STATIC_WARNING(
            //     serialization::tracking_level == serialization::track_never
            // );
            return *this << const_cast<const T &>(t);
        }
        // the & operator
        template<class T>
        Archive & operator&(T & t){
            return *this << const_cast<const T &>(t);
        }
    #endif
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_DETAIL_INTERFACE_IARCHIVE_HPP
