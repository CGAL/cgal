//-----------------------------------------------------------------------------
// boost mpl/aux_/unwrap.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
//  Copyright (c) 2001, 2002 Peter Dimov and Multi Media Ltd.
//  Copyright (c) 2001 David Abrahams
//
//  Permission to copy, use, modify, sell and distribute this software
//  is granted provided this copyright notice appears in all copies.
//  This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.

#ifndef BOOST_MPL_AUX_UNWRAP_HPP_INCLUDED
#define BOOST_MPL_AUX_UNWRAP_HPP_INCLUDED

#include "boost/ref.hpp"

namespace boost { namespace mpl { namespace aux {

template< typename F >
inline
F& unwrap(F& f, long)
{
    return f;
}

template< typename F >
inline
F&
unwrap(reference_wrapper<F>& f, int)
{
    return f;
}

template< typename F >
inline
F&
unwrap(reference_wrapper<F> const& f, int)
{
    return f;
}

}}} // namespace boost::mpl::aux

#endif // BOOST_MPL_AUX_UNWRAP_HPP_INCLUDED
