
#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_PLACEHOLDERS_HPP_INCLUDED
#define BOOST_MPL_PLACEHOLDERS_HPP_INCLUDED

// + file: boost/mpl/placeholders.hpp
// + last modified: 25/feb/03

// Copyright (c) 2001-03
// Peter Dimov, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/arg.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
 && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER placeholders.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/aux_/config/nttp.hpp"
#   include "boost/mpl/limits/arity.hpp"
#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/cat.hpp"

namespace boost {
namespace mpl {

// watch out for GNU gettext users, who #define _(x)
#if !defined(_) || defined(BOOST_MPL_NO_UNNAMED_PLACEHOLDER_SUPPORT)
typedef arg<-1> _;

namespace placeholders {
using boost::mpl::_;
}
#endif

//: agurt, 17/mar/02: one more placeholder for the last 'apply#' 
//: specialization
#define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(1, BOOST_MPL_METAFUNCTION_MAX_ARITY + 1, "boost/mpl/placeholders.hpp"))
#include BOOST_PP_ITERATE()

} // namespace mpl
} // namespace boost 

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_PLACEHOLDERS_HPP_INCLUDED

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

typedef arg<i> BOOST_PP_CAT(_,i);

namespace placeholders {
using boost::mpl::BOOST_PP_CAT(_,i);
}

#undef i
#endif // BOOST_PP_IS_ITERATING
