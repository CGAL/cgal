//  (C) Copyright Gennadiy Rozental 2001-2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : simple facility that mimmic notion of read-only read-write 
//  properties in C++ classes. Original idea by Henrik Ravn.
// ***************************************************************************

#ifndef BOOST_TEST_CLASS_PROPERTIES_HPP_071894GER
#define BOOST_TEST_CLASS_PROPERTIES_HPP_071894GER

// Boost.Test
#include <boost/test/detail/unit_test_config.hpp>

// BOOST
#include <boost/preprocessor/repetition/repeat.hpp> 
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/call_traits.hpp>
#include <boost/type_traits/add_pointer.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/utility/addressof.hpp>

// STL
#include <iosfwd>

#if BOOST_WORKAROUND(__BORLANDC__, <= 0x570)           || \
    BOOST_WORKAROUND( __COMO__, <= 0x433 )             || \
    BOOST_WORKAROUND( __INTEL_COMPILER, <= 800 )       || \
    BOOST_WORKAROUND(__GNUC__, < 3)                    || \
    defined(__sgi) && _COMPILER_VERSION <= 730         || \
    BOOST_WORKAROUND(__IBMCPP__, BOOST_TESTED_AT(600)) || \
    defined(__DECCXX)
#define BOOST_TEST_NO_PROTECTED_USING
#endif

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************                 class_property               ************** //
// ************************************************************************** //

template<class PropertyType>
class class_property {
protected:
    typedef typename call_traits<PropertyType>::const_reference read_access_t;
    typedef typename call_traits<PropertyType>::param_type      write_param_t;
    typedef typename add_pointer<PropertyType const>::type      address_res_t;

public:
    // Constructor
                    class_property() : value( PropertyType() ) {}
    explicit        class_property( write_param_t init_value ) : value( init_value ) {}

    // Access methods
    operator        read_access_t() const   { return value; }
    read_access_t   get() const             { return value; }
    bool            operator!() const       { return !value; }
    address_res_t   operator&() const       { return &value; }

    // Data members
#ifndef BOOST_TEST_NO_PROTECTED_USING
protected:
#endif
    PropertyType        value;
};

//____________________________________________________________________________//

#ifdef BOOST_CLASSIC_IOSTREAMS

template<class PropertyType>
inline std::ostream&
operator<<( std::ostream& os, class_property<PropertyType> const& p )

#else

template<typename CharT1, typename Tr,class PropertyType>
inline std::basic_ostream<CharT1,Tr>&
operator<<( std::basic_ostream<CharT1,Tr>& os, class_property<PropertyType> const& p )

#endif
{
    return os << p.get();
}

//____________________________________________________________________________//

#define DEFINE_PROPERTY_FREE_BINARY_OPERATOR( op )                              \
template<class PropertyType>                                                    \
inline bool                                                                     \
operator op( PropertyType const& lhs, class_property<PropertyType> const& rhs ) \
{                                                                               \
    return lhs op rhs.get();                                                    \
}                                                                               \
template<class PropertyType>                                                    \
inline bool                                                                     \
operator op( class_property<PropertyType> const& lhs, PropertyType const& rhs ) \
{                                                                               \
    return lhs.get() op rhs;                                                    \
}                                                                               \
template<class PropertyType>                                                    \
inline bool                                                                     \
operator op( class_property<PropertyType> const& lhs,                           \
             class_property<PropertyType> const& rhs )                          \
{                                                                               \
    return lhs.get() op rhs.get();                                              \
}                                                                               \
/**/

DEFINE_PROPERTY_FREE_BINARY_OPERATOR( == )
DEFINE_PROPERTY_FREE_BINARY_OPERATOR( != )

#undef DEFINE_PROPERTY_FREE_BINARY_OPERATOR

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1200)

#define DEFINE_PROPERTY_LOGICAL_OPERATOR( op )                                  \
template<class PropertyType>                                                    \
inline bool                                                                     \
operator op( bool b, class_property<PropertyType> const& p )                    \
{                                                                               \
    return b op p.get();                                                        \
}                                                                               \
template<class PropertyType>                                                    \
inline bool                                                                     \
operator op( class_property<PropertyType> const& p, bool b )                    \
{                                                                               \
    return b op p.get();                                                        \
}                                                                               \
/**/

DEFINE_PROPERTY_LOGICAL_OPERATOR( && )
DEFINE_PROPERTY_LOGICAL_OPERATOR( || )

#endif

// ************************************************************************** //
// **************               readonly_property              ************** //
// ************************************************************************** //

template<class PropertyType>
class readonly_property : public class_property<PropertyType> {
    typedef class_property<PropertyType>    base;
    typedef typename base::address_res_t    arrow_res_t;
protected:
    typedef typename base::write_param_t    write_param_t;
public:
    // Constructor
                    readonly_property() {}
    explicit        readonly_property( write_param_t init_value ) : base( init_value ) {}

    // access methods
    arrow_res_t     operator->() const      { return boost::addressof( base::value ); }
};

//____________________________________________________________________________//

#if BOOST_WORKAROUND(__IBMCPP__, BOOST_TESTED_AT(600))

#define BOOST_READONLY_PROPERTY( property_type, friends ) boost::unit_test::readwrite_property<property_type >

#else

#define BOOST_READONLY_PROPERTY_DECLARE_FRIEND(r, data, elem) friend class elem;

#define BOOST_READONLY_PROPERTY( property_type, friends )                           \
class BOOST_JOIN( readonly_property, __LINE__ )                                     \
: public boost::unit_test::readonly_property<property_type > {                      \
    typedef boost::unit_test::readonly_property<property_type > base;               \
    BOOST_PP_SEQ_FOR_EACH( BOOST_READONLY_PROPERTY_DECLARE_FRIEND, ' ', friends )   \
    typedef base::write_param_t  write_param_t;                                     \
public:                                                                             \
                BOOST_JOIN( readonly_property, __LINE__ )() {}                      \
    explicit    BOOST_JOIN( readonly_property, __LINE__ )( write_param_t init_v  )  \
    : base( init_v ) {}                                                             \
}                                                                                   \
/**/

#endif

// ************************************************************************** //
// **************              readwrite_property              ************** //
// ************************************************************************** //

template<class PropertyType>
class readwrite_property : public class_property<PropertyType> {
    typedef class_property<PropertyType>                base;
    typedef typename add_pointer<PropertyType>::type    arrow_res_t;
    typedef typename base::address_res_t                const_arrow_res_t;
    typedef typename base::write_param_t                write_param_t;
public:
                    readwrite_property() : base() {}
    explicit        readwrite_property( write_param_t init_value ) : base( init_value ) {}

    // access methods
    void            set( write_param_t v )  { base::value = v; }
    arrow_res_t     operator->()            { return boost::addressof( base::value ); }
    const_arrow_res_t operator->() const    { return boost::addressof( base::value ); }

#ifndef BOOST_TEST_NO_PROTECTED_USING
    using           base::value;
#endif
};

//____________________________________________________________________________//

} // unit_test

} // namespace boost

#undef BOOST_TEST_NO_PROTECTED_USING

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:19  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.28  2004/10/05 07:38:09  rogeeff
//  typo fix
//
//  Revision 1.27  2004/10/05 02:18:18  rogeeff
//  gcc2.85 fix
//
//  Revision 1.26  2004/09/19 09:22:12  rogeeff
//  ios fix for classic iostreams
//
//  Revision 1.25  2004/08/10 04:08:30  rogeeff
//  first tru64cxx65 fix
//
//  Revision 1.24  2004/08/04 04:14:58  rogeeff
//  irix6 fix
//
//  Revision 1.23  2004/07/19 12:19:05  rogeeff
//  guard rename
//  no using workaroung reworked
//
//  Revision 1.22  2004/06/05 11:00:26  rogeeff
//  proper IBM VA port
//
//  Revision 1.21  2004/05/27 06:22:17  rogeeff
//  workaround for gcc 2.95 io
//  workaround for msvc logical properties operators
//
//  Revision 1.20  2004/05/25 10:16:22  rogeeff
//  upgrade workaround version for Intel
//
//  Revision 1.19  2004/05/23 08:58:12  rogeeff
//  add intel into workaround branch
//
//  Revision 1.18  2004/05/23 08:56:58  rogeeff
//  add intel into workaround branch
//
//  Revision 1.17  2004/05/21 06:19:11  rogeeff
//  hack for non-using version of readwrite properties
//  licence update
//
//  Revision 1.16  2004/05/18 13:39:32  dgregor
//  class_properties.hpp: Make the empty character constant into a single space (which isn't used), because the Sun compiler is very eager to spit out an error here.
//
//  Revision 1.15  2004/05/18 13:10:49  dgregor
//  class_properties.hpp: Borland C++ does not handle using declarations for data members properly; fixed the existing Borland workaround.
//
//  Revision 1.14  2004/05/11 11:00:53  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.13  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_TEST_CLASS_PROPERTIES_HPP_071894GER
