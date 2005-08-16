#ifndef  BOOST_SERIALIZATION_COLLECTIONS_LOAD_IMP_HPP
#define BOOST_SERIALIZATION_COLLECTIONS_LOAD_IMP_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

#if defined(_MSC_VER) && (_MSC_VER <= 1020)
#  pragma warning (disable : 4786) // too long name, harmless warning
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// collections_load_imp.hpp: serialization for loading stl collections

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

// helper function templates for serialization of collections

#include <cassert>
#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

#include <boost/aligned_storage.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/serialization.hpp>

namespace boost{
namespace serialization {

// reserve space on stack for an object of type T without actually
// construction such an object
template<typename T > 
struct stack_allocate
{
    T * address() {
        return static_cast<T*>(storage_.address()); 
    }
    T & reference() {
        return * address();
    }
private:
    typedef BOOST_DEDUCED_TYPENAME boost::aligned_storage<
        sizeof(T), 
        #if BOOST_WORKAROUND(__BORLANDC__,BOOST_TESTED_AT(0x560))
            8
        #else
            boost::alignment_of<T>::value
        #endif
    > type;
    type storage_;
};

// construct element on the stack
template<class Archive, class T>
struct stack_construct : public stack_allocate<T>
{
    stack_construct(Archive & ar){
        // note borland emits a no-op without the explicit namespace
        boost::serialization::load_construct_data_adl(ar, this->address(), 0U);
    }
    ~stack_construct(){
        this->address()->~T(); // undo load_construct_data above
    }
};

namespace stl {

//////////////////////////////////////////////////////////////////////
// implementation of serialization for STL containers
//

// sequential container input
template<class Archive, class Container>
struct archive_input_seq
{
    inline void operator()(Archive &ar, Container &s)
    {
        typedef BOOST_DEDUCED_TYPENAME Container::value_type type;
        stack_construct<Archive, type> t(ar);
        // borland fails silently w/o full namespace
        ar >> boost::serialization::make_nvp("item", t.reference());
        s.push_back(t.reference());
        ar.reset_object_address(& s.back() , & t.reference());
    }
};

// map input
template<class Archive, class Container>
struct archive_input_map
{
    inline void operator()(Archive &ar, Container &s)
    {
        #if BOOST_WORKAROUND(BOOST_DINKUMWARE_STDLIB, == 1)
            typedef BOOST_DEDUCED_TYPENAME std::pair<
                BOOST_DEDUCED_TYPENAME Container::key_type,
                BOOST_DEDUCED_TYPENAME Container::referent_type
            > type;
        #else
            typedef BOOST_DEDUCED_TYPENAME std::pair<
                BOOST_DEDUCED_TYPENAME Container::key_type,
                BOOST_DEDUCED_TYPENAME Container::mapped_type
            > type;
        #endif
        stack_construct<Archive, type> t(ar);
        // borland fails silently w/o full namespace
        ar >> boost::serialization::make_nvp("item", t.reference());
        std::pair<BOOST_DEDUCED_TYPENAME Container::const_iterator, bool> result = 
            s.insert(t.reference());
        assert(result.second); // make sure we inserted a new element
        ar.reset_object_address(& (* result.first), & t);
    }
};

// set input
template<class Archive, class Container>
struct archive_input_set
{
    inline void operator()(Archive &ar, Container &s)
    {   
        typedef BOOST_DEDUCED_TYPENAME Container::value_type type;
        stack_construct<Archive, type> t(ar);
        // borland fails silently w/o full namespace
        ar >> boost::serialization::make_nvp("item", t.reference());
        std::pair<BOOST_DEDUCED_TYPENAME Container::const_iterator, bool> result = 
            s.insert(t.reference());
        assert(result.second); // make sure we inserted a new element
        ar.reset_object_address(& (* result.first), & t);
    }
};

// multimap input
template<class Archive, class Container>
struct archive_input_multimap
{
    inline void operator()(Archive &ar, Container &s)
    {
        #if BOOST_WORKAROUND(BOOST_DINKUMWARE_STDLIB, == 1)
            typedef BOOST_DEDUCED_TYPENAME std::pair<
                BOOST_DEDUCED_TYPENAME Container::key_type,
                BOOST_DEDUCED_TYPENAME Container::referent_type
            > type;
        #else
            typedef BOOST_DEDUCED_TYPENAME std::pair<
                BOOST_DEDUCED_TYPENAME Container::key_type,
                BOOST_DEDUCED_TYPENAME Container::mapped_type
            > type;
        #endif
        stack_construct<Archive, type> t(ar);
        // borland fails silently w/o full namespace
        ar >> boost::serialization::make_nvp("item", t.reference());
        BOOST_DEDUCED_TYPENAME Container::const_iterator result 
            = s.insert(t.reference());
        ar.reset_object_address(& (* result), & t);
    }
};

// multiset input
template<class Archive, class Container>
struct archive_input_multiset
{
    inline void operator()(Archive &ar, Container &s)
    {   
        typedef BOOST_DEDUCED_TYPENAME Container::value_type type;
        stack_construct<Archive, type> t(ar);
        // borland fails silently w/o full namespace
        ar >> boost::serialization::make_nvp("item", t.reference());
        BOOST_DEDUCED_TYPENAME Container::const_iterator result 
            = s.insert(t.reference());
        ar.reset_object_address(& (* result), & t);
    }
};

template<class Container>
class reserve_imp
{
public:
    void operator()(Container &s, unsigned int count) const {
        s.reserve(count);
    }
};

template<class Container>
class no_reserve_imp
{
public:
    void operator()(Container & /* s */, unsigned int /* count */) const{}
};

template<class Archive, class Container, class InputFunction, class R>
inline void rebuild_collection(Archive & ar, Container &s)
{
    s.clear();
    // retrieve number of elements
    unsigned int count;
    ar >> BOOST_SERIALIZATION_NVP(count);
    R rx;
    rx(s, count);
    InputFunction ifunc;
    while(count-- > 0){
        ifunc(ar, s);
    }
}

template<class Archive, class Container>
inline void copy_collection(Archive & ar, Container &s)
{
    // retrieve number of elements
    unsigned int count;
    ar >> BOOST_SERIALIZATION_NVP(count);
    assert(count == s.size());
    BOOST_DEDUCED_TYPENAME Container::iterator it = s.begin();
    while(count-- > 0){
        ar >> boost::serialization::make_nvp("item", *it++);
    }
}

template<class Archive, class Container, class InputFunction, class R>
inline void load_collection(Archive & ar, Container &s){
//    if(0 != (ar.get_flags() & boost::archive::no_object_creation))
//        copy_collection<Archive, Container>(ar, s);
//    else
        rebuild_collection<Archive, Container, InputFunction, R>(ar, s);
}


} // namespace stl 
} // namespace serialization
} // namespace boost

#endif //BOOST_SERIALIZATION_COLLECTIONS_LOAD_IMP_HPP
