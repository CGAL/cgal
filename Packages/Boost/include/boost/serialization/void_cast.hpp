#ifndef  BOOST_SERIALIZATION_VOID_CAST_HPP
#define BOOST_SERIALIZATION_VOID_CAST_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// void_cast.hpp:   interface for run-time casting of void pointers.

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// gennadiy.rozental@tfn.com

//  See http://www.boost.org for updates, documentation, and revision history.

#include <cassert>

#include <boost/smart_cast.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>

#include <boost/type_traits/is_polymorphic.hpp>

#include <boost/serialization/force_include.hpp>
#include <boost/serialization/type_info_implementation.hpp>

namespace boost { 
namespace serialization { 

class extended_type_info;

// Given a void *, assume that it really points to an instance of one type
// and alter it so that it would point to an instance of a related type.
// Return the altered pointer. If there exists no sequence of casts that
// can transform from_type to to_type, return a NULL.  

void const *
void_upcast(
    extended_type_info const & derived_type,  
    extended_type_info const & base_type, 
    void const * t,
    bool top = true
);

inline void *
void_upcast(
    extended_type_info const & derived_type_,
    extended_type_info const & base_type_,
    void * t 
){
    return const_cast<void*>(void_upcast(
        derived_type_, 
        base_type_, 
        const_cast<void const *>(t)
    ));
}

void const *
void_downcast(
    extended_type_info const & derived_type,  
    extended_type_info const & base_type, 
    void const * t,
    bool top = true
);

inline void *
void_downcast(
    extended_type_info const & derived_type_,
    extended_type_info const & base_type_,
    void * t 
){
    return const_cast<void*>(void_downcast(
        derived_type_, 
        base_type_, 
        const_cast<void const *>(t)
    ));
}


namespace void_cast_detail {

// note: can't be abstract because an instance is used as a search argument
class void_caster
{
    friend struct void_caster_compare ;
    friend const void * boost::serialization::void_upcast(
        const extended_type_info & derived_type,
        const extended_type_info & base_type,
        const void * t,
        bool top
    );
    friend const void * boost::serialization::void_downcast(
        const extended_type_info & derived_type,
        const extended_type_info & base_type,
        const void * t,
        bool top
    );
    // each derived class must re-implement these;
    virtual void const * upcast(void const * t) const {
        assert(false);
        return NULL;
    }
    virtual void const * downcast(void const * t) const {
        assert(false);
        return NULL;
    }
    // Data members
    extended_type_info const & m_derived_type;
    extended_type_info const & m_base_type;
protected:
    void self_register();
public:
    // Constructor
    void_caster(
        extended_type_info const & derived_type_,
        extended_type_info const & base_type_ 
    ) :
        m_derived_type( derived_type_),
        m_base_type(base_type_)
    {}
    virtual ~void_caster(){};
};

class void_caster_derived : public void_caster
{
    std::ptrdiff_t difference;
    virtual void const*
    upcast( void const* t ) const{
        return static_cast<const char*> ( t ) + difference;
    }
    virtual void const*
    downcast( void const* t ) const{
        return static_cast<const char*> ( t ) - difference;
    }
public:
    void_caster_derived(
        extended_type_info const& derived_type_,
        extended_type_info const& base_type_,
        std::ptrdiff_t difference_
    ) :
        void_caster(derived_type_, base_type_),
        difference( difference_ )
    {
        self_register();
    }
};
    
template <class Derived, class Base>
class void_caster_primitive : public void_caster
{
    virtual void const* downcast( void const * t ) const {
        return boost::smart_cast<const Derived *>(
            boost::smart_cast<const Base *>(t)
        );
    }
    virtual void const* upcast(void const * t) const {
        return boost::smart_cast<const Base *>(
            boost::smart_cast<const Derived *>(t)
        );
    }

public:
    void_caster_primitive() :
        void_caster( 
            * type_info_implementation<Derived>::type::get_instance(), 
            * type_info_implementation<Base>::type::get_instance() 
        )
    {
        self_register();
    }
};

// this purpose of this class is to create to->from and from->to instances
// of void_caster_primitive for each related pair of types.  This has to be
// done a pre-execution time - hence the usage of static variable.

template<class Derived, class Base>
struct static_initializer
{
    static const static_initializer instance;
    static_initializer(){
        static void_caster_primitive<const Derived, const Base> instance1;
    }
    static const static_initializer & instantiate(){
        return instance;
    }
};

// This implementation is somewhat complicated by the following problem 
// with msvc 6.0

// when used with msvc 6.0, the conforming method emits error at link
// time LIN1179 - "invalid or corrupt file: duplicate comdat".  According
// to http://groups.google.com/groups?th=8a05c82c4ffee280 (look for P78)
//  A LNK1179 error occurs if these preconditions are met:
//  - The template class takes at least two arguments.
//  - The template is used at least two times with identical first
//    and different second arguments.
//  - The static member variable is of an object type with at least one
//    base class. (In another scenario it also occurred using a member
//    without a base class.)

// in the absence of multiple inheritance, a derived class will be only
// used once as the from element of a pair.  The effectively avoids the
// problem described abot.

// this problem will manifest itself in the serialization system when
// serializing an instance of a class with multiple base classes when
// using msvc 6.0.

// note: this instance is defined at compile time.  Compiler processing of
// void_cast_register instantiates this class for the base/derived pair so
// that all registrations occur before the program actually starts.  note
// that using new to register a new instance would not guarentee what
// the registration occurs soon enough
// just use correct C++ syntax
template<class Derived, class Base>
const static_initializer<Derived, Base> static_initializer<Derived, Base>::instance;

} // void_cast_detail 

// Register a base/derived pair.  This indicates that it is possible
// to upcast a void pointer from Derived to Base and downcast a
// void pointer from Base to Derived.  Note bogus arguments to workaround
// bug in msvc 6.0
template<class Derived, class Base>
inline BOOST_DLLEXPORT const void *  void_cast_register(
    const Derived * /* dnull = NULL */, 
    const Base * /* bnull = NULL */
) BOOST_USED;
template<class Derived, class Base>
inline BOOST_DLLEXPORT const void * void_cast_register(
    const Derived * /* dnull = NULL */, 
    const Base * /* bnull = NULL */
){
    return & void_cast_detail::static_initializer<Derived, Base>::instantiate();
}

} // namespace serialization
} // namespace boost

#endif // BOOST_SERIALIZATION_VOID_CAST_HPP
