#ifndef BOOST_ARCHIVE_DETAIL_INTERFACE_IARCHIVE_HPP
#define BOOST_ARCHIVE_DETAIL_INTERFACE_IARCHIVE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// interface_iarchive.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.
#include <string>
#include <boost/config.hpp>
#include <boost/cstdint.hpp>
#include <boost/mpl/bool.hpp>

#include <boost/serialization/nvp.hpp>
#include <boost/archive/detail/iserializer.hpp>

namespace boost { 
namespace archive {
namespace detail {

class basic_iserializer;
class basic_pointer_iserializer;

template<class Archive>
class interface_iarchive 
{
protected:
    interface_iarchive(){};
public:
    /////////////////////////////////////////////////////////
    // archive public interface
    struct is_loading {
        typedef mpl::bool_<true> type;
        BOOST_STATIC_CONSTANT(bool, value=true);
    };
    struct is_saving {
        typedef mpl::bool_<false> type;
        BOOST_STATIC_CONSTANT(bool, value=false);
    };

    // return a pointer to the most derived class
    Archive * This(){
        return static_cast<Archive *>(this);
    }

    template<class T>
    const basic_pointer_iserializer * register_type(T * t = NULL){
        const basic_pointer_iserializer & bpis =
            archive::detail::instantiate_pointer_iserializer(
                static_cast<Archive *>(NULL),
                static_cast<T *>(NULL)
            );
        this->This()->register_basic_serializer(bpis.get_basic_serializer());
        return & bpis;
    }

    // default processing - invoke serialization library
    template<class T>
    void load_override(T & t, /*BOOST_PFTO*/ int){
        archive::load(* this->This(), t);
    }

    // define operators for non-const arguments.  Don't depend one the const
    // ones below because the compiler MAY make a temporary copy to
    // create the const parameter (Though I havn't seen this happen). 
    // the >> operator
    template<class T>
    Archive & operator>>(T & t){
        // if this assertion trips. It means we're trying to load a
        // const object with a compiler that doesn't have correct
        // funtion template ordering.  On other compilers, this is
        // handled below.
        BOOST_STATIC_ASSERT(! boost::is_const<T>::value);
        this->This()->load_override(t, 0);
        return * this->This();
    }

    // the & operator 
    template<class T>
    Archive & operator&(T & t){
        // see above
        BOOST_STATIC_ASSERT(! boost::is_const<T>::value);
        this->This()->load_override(t, 0);
        return * this->This();
    }

    // define the following pair in order to permit passing of const and non_const
    // temporary objects. These are needed to properly implement serialization
    // wrappers.

    #ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    // the >> operator
    template<class T>
    Archive & operator>>(const T & t){
        // this should only be used for wrappers.  Check that here
        This()->load_override(const_cast<T &>(t), 0);
        return * this->This();
    }
    // the & operator 
    template<class T>
    Archive & operator&(const T & t){
        return * this >> t;
    }
    #endif
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_DETAIL_INTERFACE_IARCHIVE_HPP
