//-----------------------------------------------------------------------------
// boost mpl/aux_/iter_fold_impl.hpp header file
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

#ifndef BOOST_MPL_AUX_ITER_FOLD_IMPL_HPP_INCLUDED
#define BOOST_MPL_AUX_ITER_FOLD_IMPL_HPP_INCLUDED

#include "boost/mpl/aux_/apply.hpp"
#include "boost/mpl/aux_/next.hpp"
#include "boost/mpl/aux_/config/eti.hpp"
#include "boost/config.hpp"

#if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) && \
    !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/if.hpp"
#   include "boost/type_traits/is_same.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) && \
    !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER iter_fold_impl.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   define BOOST_MPL_AUX_FOLD_IMPL_OP(iter) iter
#   define BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX iter_fold
#   include "boost/mpl/aux_/fold_impl_body.hpp"

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_ITER_FOLD_IMPL_HPP_INCLUDED
