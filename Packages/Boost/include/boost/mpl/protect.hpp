//-----------------------------------------------------------------------------
// boost mpl/protect.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Peter Dimov, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_PROTECT_HPP_INCLUDED
#define BOOST_MPL_PROTECT_HPP_INCLUDED

#include "boost/mpl/aux_/arity.hpp"
#include "boost/mpl/aux_/config/dtp.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    >
struct protect : T
{
    typedef struct protect type;
};

#if defined(BOOST_BROKEN_DEFAULT_TEMPLATE_PARAMETERS_IN_NESTED_TEMPLATES)
namespace aux { 
template< BOOST_MPL_AUX_NTTP_DECL(int, N), typename T >
struct arity< protect<T>, N > 
    : arity<T,N>
{ 
};
} // namespace aux
#endif

BOOST_MPL_AUX_VOID_SPEC_MAIN(1, protect)
#if !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
BOOST_MPL_AUX_VOID_SPEC_TEMPLATE_ARITY(1, 1, protect)
#endif

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_PROTECT_HPP_INCLUDED
