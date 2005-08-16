//  (C) Copyright John Maddock 2005. 
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_TYPE_TRAITS_INTEGRAL_CONSTANT_HPP
#define BOOST_TYPE_TRAITS_INTEGRAL_CONSTANT_HPP

#include <boost/config.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/integral_c.hpp>

namespace boost{

#if defined(BOOST_NO_DEPENDENT_TYPES_IN_TEMPLATE_VALUE_PARAMETERS) || defined(__BORLANDC__)
template <class T, int val>
#else
template <class T, T val>
#endif
struct integral_constant : public mpl::integral_c<T, val>
{
   //BOOST_STATIC_CONSTANT(T, value = val);
   //typedef T value_type;
   typedef integral_constant<T,val> type;

#if 0
   //
   // everything that follows now, is MPL-compatibility code:
   //
   typedef ::boost::mpl::integral_c_tag tag;

   // have to #ifdef here: some compilers don't like the 'val + 1' form (MSVC),
   // while some other don't like 'value + 1' (Borland), and some don't like
   // either
#if BOOST_WORKAROUND(__EDG_VERSION__, <= 243)
private:
   BOOST_STATIC_CONSTANT(T, next_value = BOOST_MPL_AUX_STATIC_CAST(T, (val + 1)));
   BOOST_STATIC_CONSTANT(T, prior_value = BOOST_MPL_AUX_STATIC_CAST(T, (val - 1)));
public:
   typedef integral_constant<T,next_value> next;
   typedef integral_constant<T,prior_value> prior;
#elif BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561)) \
   || BOOST_WORKAROUND(__IBMCPP__, BOOST_TESTED_AT(502)) \
   || BOOST_WORKAROUND(__HP_aCC, BOOST_TESTED_AT(53800))
   typedef integral_constant<T, ( BOOST_MPL_AUX_STATIC_CAST(T, (val + 1)) )> next;
   typedef integral_constant<T, ( BOOST_MPL_AUX_STATIC_CAST(T, (val - 1)) )> prior;
#else
   typedef integral_constant<T, ( BOOST_MPL_AUX_STATIC_CAST(T, (value + 1)) )> next;
   typedef integral_constant<T, ( BOOST_MPL_AUX_STATIC_CAST(T, (value - 1)) )> prior;
#endif

   // enables uniform function call syntax for families of overloaded 
   // functions that return objects of both arithmetic ('int', 'long',
   // 'double', etc.) and wrapped integral types (for an example, see 
   // "mpl/example/power.cpp")
   operator T() const { return static_cast<T>(this->value); } 
#endif
};

template<> struct integral_constant<bool,true> : public mpl::true_ 
{
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
   typedef mpl::true_ base_;
   using base_::value;
#endif
   typedef integral_constant<bool,true> type;
};
template<> struct integral_constant<bool,false> : public mpl::false_ 
{
#if BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
   typedef mpl::false_ base_;
   using base_::value;
#endif
   typedef integral_constant<bool,false> type;
};

typedef integral_constant<bool,true> true_type;
typedef integral_constant<bool,false> false_type;

}

#endif
