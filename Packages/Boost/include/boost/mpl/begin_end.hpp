//-----------------------------------------------------------------------------
// boost mpl/begin_end.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_BEGIN_END_HPP_INCLUDED
#define BOOST_MPL_BEGIN_END_HPP_INCLUDED

#include "boost/mpl/begin_end_fwd.hpp"
#include "boost/mpl/aux_/begin_end_impl.hpp"
#include "boost/mpl/aux_/sequence_tag.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

namespace boost {
namespace mpl {

// agurt, 13/sep/02: switched from inheritance to typedef; MSVC is more
// happy this way (less ETI-related errors), and it doesn't affect 
// anything else
template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    >
struct begin
{
    typedef typename BOOST_MPL_AUX_SEQUENCE_TAG(Sequence) tag_;
    typedef typename begin_traits< tag_ >
        ::template algorithm< Sequence >::type type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,begin,(Sequence))
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    >
struct end
{
    typedef typename BOOST_MPL_AUX_SEQUENCE_TAG(Sequence) tag_;
    typedef typename end_traits< tag_ >
        ::template algorithm< Sequence >::type type;
    BOOST_MPL_AUX_LAMBDA_SUPPORT(1,end,(Sequence))
};

BOOST_MPL_AUX_VOID_SPEC(1, begin)
BOOST_MPL_AUX_VOID_SPEC(1, end)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_BEGIN_END_HPP_INCLUDED
