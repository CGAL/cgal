#ifndef BOOST_SERIALIZATION_BINARY_OBJECT_HPP
#define BOOST_SERIALIZATION_BINARY_OBJECT_HPP

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

#include <cassert>

#include <cstddef> // std::size_t
#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/preprocessor/stringize.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/nvp.hpp>

namespace boost {
namespace serialization {

struct binary_object {
    /* const */ void * const m_t;
    const std::size_t m_size;
    template<class Archive>
    void save(Archive & ar, const unsigned int /* file_version */) const {
        ar.save_binary(m_t, m_size);
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int /* file_version */){
        ar.load_binary(const_cast<void *>(m_t), m_size);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
        binary_object(/* const */ void * const t, std::size_t size) :
        m_t(t),
        m_size(size)
    {}
    binary_object(const binary_object & rhs) :
        m_t(rhs.m_t),
        m_size(rhs.m_size)
    {}
};

// just a little helper to support the convention that all serialization
// wrappers follow the naming convention make_xxxxx
inline 
binary_object make_binary_object(/* const */ void * t, std::size_t size){
    return binary_object(t, size);
}

// make special version of nvp which for binary object wrapper rather than
// a pointer to them.  This permits a better composition of nvp(binary_object)
// than would otherwise be possible.

template<>
struct nvp<binary_object> : 
    public std::pair<const char *, binary_object>,
    public traits<nvp<binary_object>, object_serializable, track_never>
{
    explicit nvp(const char * name, binary_object & t) : 
        std::pair<const char *, binary_object>(name, t)
    {}
    nvp(const nvp<binary_object> & rhs) : 
        std::pair<const char *, binary_object>(rhs.first, rhs.second)
    {}

    const char * name() const {
        return this->first;
    }
    binary_object & value() {
        return this->second;
    }
    const binary_object & value() const {
        return this->second;
    }
    // default treatment for name-value pairs. The name is
    // just discarded and only the value is serialized.  Note the unusual
    // fact that his is "const".  This is because wrappers themselves are
    // in fact "const" - even though the things they wrap may not be.
    // note: removed "const" because it confuses MSVC 6.0
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */) /*const*/
    {
        ar & value();
    }
};

inline nvp<binary_object> make_nvp(const char * name, binary_object t){
    return nvp<binary_object>(name, t);
}

} // namespace serialization
} // boost

// don't need versioning info for this type
BOOST_CLASS_IMPLEMENTATION(
    binary_object, 
    boost::serialization::object_serializable
)

// don't track binary objects - usually they will be created on the stack
// and tracking algorithm (which uses the object address) might get
// confused.  note that these address will likely be members of some
// other structure which itself is tracked, so as a practical matter
// suppressing tracking shouldn't cause any redundancy.

BOOST_CLASS_TRACKING(binary_object, boost::serialization::track_never) 

#endif // BOOST_SERIALIZATION_BINARY_OBJECT_HPP
