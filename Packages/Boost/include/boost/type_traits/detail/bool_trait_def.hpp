
// NO INCLUDE GUARDS, THE HEADER IS INTENDED FOR MULTIPLE INCLUSION

// Copyright Aleksey Gurtovoy 2002-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)

// $Source$
// $Date$
// $Revision$

#include <boost/type_traits/detail/template_arity_spec.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/aux_/lambda_support.hpp>
#include <boost/config.hpp>

#if defined(__SUNPRO_CC)
#   define BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    typedef BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::bool_< C > type; \
    enum { value = type::value }; \
    /**/
#   define BOOST_TT_AUX_BOOL_C_BASE(C)

#elif defined(BOOST_MSVC) && BOOST_MSVC <= 1200

#   define BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    typedef BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::bool_< C > base_; \
    using base_::value; \
    /**/

#endif

#ifndef BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL
#   define BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) /**/
#endif

#ifndef BOOST_TT_AUX_BOOL_C_BASE
#   define BOOST_TT_AUX_BOOL_C_BASE(C) : BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::bool_< C >
#endif 


#define BOOST_TT_AUX_BOOL_TRAIT_DEF1(trait,T,C) \
template< typename T > struct trait \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,trait,(T)) \
}; \
\
BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(1,trait) \
/**/


#define BOOST_TT_AUX_BOOL_TRAIT_DEF2(trait,T1,T2,C) \
template< typename T1, typename T2 > struct trait \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,trait,(T1,T2)) \
}; \
\
BOOST_TT_AUX_TEMPLATE_ARITY_SPEC(2,trait) \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_SPEC1(trait,sp,C) \
template<> struct trait< sp > \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(1,trait,(sp)) \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_SPEC2(trait,sp1,sp2,C) \
template<> struct trait< sp1,sp2 > \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(2,trait,(sp1,sp2)) \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC1(trait,sp,C) \
template<> struct trait##_impl< sp > \
{ \
    BOOST_STATIC_CONSTANT(bool, value = (C)); \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_IMPL_SPEC2(trait,sp1,sp2,C) \
template<> struct trait##_impl< sp1,sp2 > \
{ \
    BOOST_STATIC_CONSTANT(bool, value = (C)); \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_PARTIAL_SPEC1_1(param,trait,sp,C) \
template< param > struct trait< sp > \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_PARTIAL_SPEC1_2(param1,param2,trait,sp,C) \
template< param1, param2 > struct trait< sp > \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_PARTIAL_SPEC2_1(param,trait,sp1,sp2,C) \
template< param > struct trait< sp1,sp2 > \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(2,trait,(sp1,sp2)) \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_PARTIAL_SPEC2_2(param1,param2,trait,sp1,sp2,C) \
template< param1, param2 > struct trait< sp1,sp2 > \
    BOOST_TT_AUX_BOOL_C_BASE(C) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(C) \
}; \
/**/

#define BOOST_TT_AUX_BOOL_TRAIT_IMPL_PARTIAL_SPEC2_1(param,trait,sp1,sp2,C) \
template< param > struct trait##_impl< sp1,sp2 > \
{ \
    BOOST_STATIC_CONSTANT(bool, value = (C)); \
}; \
/**/

#ifndef BOOST_NO_CV_SPECIALIZATIONS
#   define BOOST_TT_AUX_BOOL_TRAIT_CV_SPEC1(trait,sp,value) \
    BOOST_TT_AUX_BOOL_TRAIT_SPEC1(trait,sp,value) \
    BOOST_TT_AUX_BOOL_TRAIT_SPEC1(trait,sp const,value) \
    BOOST_TT_AUX_BOOL_TRAIT_SPEC1(trait,sp volatile,value) \
    BOOST_TT_AUX_BOOL_TRAIT_SPEC1(trait,sp const volatile,value) \
    /**/
#else
#   define BOOST_TT_AUX_BOOL_TRAIT_CV_SPEC1(trait,sp,value) \
    BOOST_TT_AUX_BOOL_TRAIT_SPEC1(trait,sp,value) \
    /**/
#endif

#if 0  // there are true_type and false_type already in boost::
       // This also induces dependencies which may be undesirable
       // Let's wait until sometime not just before a release and clean
       // the whole ct_if mess up.
# ifndef BOOST_TT_INTEGRAL_CONSTANT
#  define BOOST_TT_INTEGRAL_CONSTANT
#  include <boost/mpl/integral_c.hpp>

//
// this is not a TR1 conforming integral_constant,
// but it is a first start:
//

namespace boost{

template <class T, T val>
struct integral_constant
: public BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::integral_c<T,val> {};


template<> struct integral_constant< bool, true > \
    BOOST_TT_AUX_BOOL_C_BASE(true) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(true) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(1,integral_constant,(bool)) \
};
template<> struct integral_constant< bool, false > \
    BOOST_TT_AUX_BOOL_C_BASE(false) \
{ \
    BOOST_TT_AUX_BOOL_TRAIT_VALUE_DECL(false) \
    BOOST_MPL_AUX_LAMBDA_SUPPORT_SPEC(1,integral_constant,(bool)) \
};

namespace pending {
typedef BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::true_ true_type;
typedef BOOST_MPL_AUX_ADL_BARRIER_NAMESPACE::false_ false_type;
}

}

# endif
#endif
