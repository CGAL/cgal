#ifndef BOOST_SERIALIZATION_NVP_HPP
#define BOOST_SERIALIZATION_NVP_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// nvp.hpp: interface for serialization system.

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <utility>

#include <boost/config.hpp>
#include <boost/serialization/traits.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/tracking.hpp>

#include <boost/serialization/base_object.hpp>

namespace boost {
namespace serialization {

template<class T>
struct nvp : 
    public std::pair<const char *, T *>,
    public traits<nvp<T>, object_serializable, track_never>
{
    explicit nvp(const char * name, T & t) :
        // note: rudundant cast works around borland issue
        std::pair<const char *, T *>(name, (T*)(& t))
    {}
    nvp(const nvp & rhs) : 
        // note: rudundant cast works around borland issue
        std::pair<const char *, T *>(rhs.first, (T*)rhs.second)
    {}

    const char * name() const {
        return this->first;
    }
    T & value() const {
        return *(this->second);
    }
    // default treatment for name-value pairs. The name is
    // just discarded and only the value is serialized.  Note the unusual
    // fact that his is "const".  This is because wrappers themselves are
    // in fact "const" - even though the things they wrap may not be.
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */) const
    {
        ar & value();
    }
};

template<class T>
inline nvp<T> make_nvp(const char * name, T & t){
    return nvp<T>(name, t);
}

// to maintain efficiency and portability, we want to assign
// specific serialization traits to all instances of this wrappers.
// we can't strait forward method below as it depends upon
// Partial Template Specialization and doing so would mean that wrappers
// wouldn't be treated the same on different platforms.  This would
// break archive portability. Leave this here as reminder not to use it !!!
#if 0 // #ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

template <class T>
struct implementation_level<nvp<T> >
{
    typedef mpl::integral_c_tag tag;
    typedef mpl::int_<object_serializable> type;
    BOOST_STATIC_CONSTANT(int, value = implementation_level::type::value);
};

// nvp objects are generally created on the stack and are never tracked
template<class T>
struct tracking_level<nvp<T> >
{
    typedef mpl::integral_c_tag tag;
    typedef mpl::int_<track_never> type;
    BOOST_STATIC_CONSTANT(int, value = tracking_level::type::value);
};

#endif

} // seralization
} // boost

#include <boost/preprocessor/stringize.hpp>

#define BOOST_SERIALIZATION_NVP(name)                              \
    boost::serialization::make_nvp(BOOST_PP_STRINGIZE(name), name)
/**/

#define BOOST_SERIALIZATION_BASE_OBJECT_NVP(name)                  \
    boost::serialization::make_nvp(                                \
        BOOST_PP_STRINGIZE(name),                                  \
        boost::serialization::base_object<name >(*this)            \
    )
/**/

#endif // BOOST_SERIALIZATION_NVP_HPP
