
#ifndef BOOST_MPL_AUX_YES_NO_HPP_INCLUDED
#define BOOST_MPL_AUX_YES_NO_HPP_INCLUDED

// + file: boost/mpl/aux_/yes_no.hpp
// + last modified: 05/nov/03

// Copyright Aleksey Gurtovoy 2000-03
//
// Use, modification and distribution are subject to the Boost Software 
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy 
// at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.


#include "boost/mpl/aux_/config/workaround.hpp"
#include "boost/mpl/aux_/config/msvc.hpp"

namespace boost { namespace mpl { namespace aux {

typedef char (&no_tag)[1];
typedef char (&yes_tag)[2];

template< bool C_ > struct yes_no_tag
{
    typedef no_tag type;
};

template<> struct yes_no_tag<true>
{
    typedef yes_tag type;
};


template< long n > struct weighted_tag
{
#if !BOOST_WORKAROUND(BOOST_MSVC, == 1200)
    typedef char (&type)[n];
#else
    char buf[n];
    typedef weighted_tag type;
#endif
};

#if    BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561)) \
    || BOOST_WORKAROUND(BOOST_MSVC, == 1300)
template<> struct weighted_tag<0>
{
    typedef char (&type)[1];
};
#endif

}}} // namespace boost::mpl::aux 

#endif // BOOST_MPL_AUX_YES_NO_HPP_INCLUDED
