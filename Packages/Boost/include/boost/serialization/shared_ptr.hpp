#ifndef BOOST_SERIALIZATION_SHARED_PTR_HPP
#define BOOST_SERIALIZATION_SHARED_PTR_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// shared_ptr.hpp: serialization for boost shared pointer

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

// note: totally unadvised hack to gain access to private variables
// in shared_ptr and shared_count. Unfortunately its the only way to
// do this without changing shared_ptr and shared_count
// the best we can do is to detect a conflict here
#include <boost/config.hpp>

#ifdef BOOST_SHARED_PTR_HPP_INCLUDED
#error "include <boost/serialization/shared_ptr> first"
#endif

#define private public
#include <boost/shared_ptr.hpp>
#undef private

#include <boost/serialization/is_abstract.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/serialization/void_cast.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/tracking.hpp>

// mark base class as an (uncreatable) base class
BOOST_IS_ABSTRACT(boost::detail::sp_counted_base)

// function specializations must be defined in the appropriate
// namespace - boost::serialization
namespace boost { 
#if defined(BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP)
namespace serialization {
#else
namespace detail {
#endif

/////////////////////////////////////////////////////////////
// sp_counted_base_impl serialization

template<class Archive, class P, class D>
inline void serialize(
    Archive & /* ar */,
    boost::detail::sp_counted_base_impl<P, D> & /* t */,
    const unsigned int /*file_version*/
){
    // register the relationship between each derived class
    // its polymorphic base
    boost::serialization::void_cast_register<
        boost::detail::sp_counted_base_impl<P, D>,
        boost::detail::sp_counted_base 
    >(
        static_cast<boost::detail::sp_counted_base_impl<P, D> *>(NULL),
        static_cast<boost::detail::sp_counted_base *>(NULL)
    );
}

template<class Archive, class P, class D>
inline void save_construct_data(
    Archive & ar,
    const boost::detail::sp_counted_base_impl<P, D> *t, 
    const unsigned int /* file_version */
){
    // variables used for construction
    ar << boost::serialization::make_nvp("ptr", t->ptr);
}

template<class Archive, class P, class D>
inline void load_construct_data(
    Archive & ar,
    boost::detail::sp_counted_base_impl<P, D> * t, 
    const unsigned int /* file_version */
){
    P ptr_;
    ar >> boost::serialization::make_nvp("ptr", ptr_);
    ::new(t)boost::detail::sp_counted_base_impl<P, D>(ptr_,  D()); // placement
    // compensate for that fact that a new shared count always is 
    // initialized with one. the add_ref_copy below will increment it
    // every time its serialized so without this adjustment
    // the use and weak counts will be off by one.
    t->use_count_ = 0;
}

/////////////////////////////////////////////////////////////
// shared_count serialization

template<class Archive>
inline void save(
    Archive & ar,
    const boost::detail::shared_count &t,
    const unsigned int /* file_version */
){
    ar << boost::serialization::make_nvp("pi", t.pi_);
}

template<class Archive>
inline void load(
    Archive & ar,
    boost::detail::shared_count &t,
    const unsigned int /* file_version */
){
    ar >> boost::serialization::make_nvp("pi", t.pi_);
    if(NULL != t.pi_)
        t.pi_->add_ref_copy();
}

template<class Archive>
inline void serialize(
    Archive & ar,
    boost::detail::shared_count &t,
    const unsigned int file_version
){
    boost::serialization::split_free(ar, t, file_version);
}

#if ! defined(BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP)
} // detail
#endif

/////////////////////////////////////////////////////////////
// implement serialization for shared_ptr<T>

template<class Archive, class T>
inline void serialize(
    Archive & ar,
    boost::shared_ptr<T> &t,
    const unsigned int /* file_version */
){
    // correct shared_ptr serialization depends upon object tracking
    // being used.
    BOOST_STATIC_ASSERT(
        boost::serialization::tracking_level<T>::value
        != boost::serialization::track_never
    );
    // only the raw pointer has to be saved
    // the ref count is maintained automatically as shared pointers are loaded
    ar.register_type(static_cast<
        boost::detail::sp_counted_base_impl<T *, boost::checked_deleter<T> > *
    >(NULL));
    ar & boost::serialization::make_nvp("px", t.px);
    ar & boost::serialization::make_nvp("pn", t.pn);
}

#if defined(BOOST_NO_ARGUMENT_DEPENDENT_LOOKUP)
} // serialization
#endif
} // namespace boost

// This macro is used to export GUIDS for shared pointers to allow
// the serialization system to export them properly. David Tonge
#define BOOST_SHARED_POINTER_EXPORT_GUID(T, K)    \
    typedef boost::detail::sp_counted_base_impl<  \
        T *,                                      \
        boost::checked_deleter< T >               \
    > __shared_ptr_ ## T;                         \
    BOOST_CLASS_EXPORT_GUID(__shared_ptr_ ## T, "__shared_ptr_" K)\
    BOOST_CLASS_EXPORT_GUID(T, K)                 \
    /**/

#define BOOST_SHARED_POINTER_EXPORT(T)            \
    BOOST_SHARED_POINTER_EXPORT_GUID(             \
        T,                                        \
        BOOST_PP_STRINGIZE(T)                     \
    )                                             \
    /**/

#endif // BOOST_SERIALIZATION_SHARED_PTR_HPP
