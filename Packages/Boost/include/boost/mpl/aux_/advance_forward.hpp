
#if !defined(BOOST_PP_IS_ITERATING)

///// header body

#ifndef BOOST_MPL_AUX_ADVANCE_FORWARD_HPP_INCLUDED
#define BOOST_MPL_AUX_ADVANCE_FORWARD_HPP_INCLUDED

// + file: boost/mpl/aux_/advance_forward.hpp
// + last modified: 06/aug/03

// Copyright (c) 2000-03
// Aleksey Gurtovoy
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

#include "boost/mpl/aux_/apply.hpp"
#include "boost/mpl/aux_/next.hpp"
#include "boost/mpl/aux_/config/eti.hpp"

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if    !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) \
    && !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER advance_forward.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/limits/unrolling.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"

#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/cat.hpp"
#   include "boost/preprocessor/inc.hpp"

namespace boost {
namespace mpl {
namespace aux {

// forward declaration
template< BOOST_MPL_AUX_NTTP_DECL(long, N) > struct advance_forward;

#   define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_UNROLLING_LIMIT, "boost/mpl/aux_/advance_forward.hpp"))
#   include BOOST_PP_ITERATE()

// implementation for N that exceeds BOOST_MPL_UNROLLING_LIMIT
template< BOOST_MPL_AUX_NTTP_DECL(long, N) > 
struct advance_forward
{
    template< typename Iterator > struct apply
    {
        typedef typename BOOST_MPL_AUX_APPLY1(
              advance_forward<BOOST_MPL_UNROLLING_LIMIT>
            , Iterator
            )::type chunk_result_;

        typedef typename BOOST_MPL_AUX_APPLY1(
              advance_forward<(
                (N - BOOST_MPL_UNROLLING_LIMIT) < 0
                    ? 0
                    : N - BOOST_MPL_UNROLLING_LIMIT
                    )>
            , chunk_result_
            )::type type;
    };
};

} // namespace aux
} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_ADVANCE_FORWARD_HPP_INCLUDED

///// iteration, depth == 1

#elif BOOST_PP_ITERATION_DEPTH() == 1
#define i BOOST_PP_FRAME_ITERATION(1)

template<>
struct advance_forward< BOOST_PP_FRAME_ITERATION(1) >
{
    template< typename Iterator > struct apply
    {
        typedef Iterator iter0;

#if i > 0
#   define BOOST_PP_ITERATION_PARAMS_2 \
    (3,(1, i, "boost/mpl/aux_/advance_forward.hpp"))
#   include BOOST_PP_ITERATE()
#endif
        typedef BOOST_PP_CAT(iter,i) type;
    };

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
    //: ETI workaround
    template<> struct apply<int>
    {
        typedef int type;
    };
#endif
};

#undef i

///// iteration, depth == 2

#elif BOOST_PP_ITERATION_DEPTH() == 2

#   define AUX_ITER_0 BOOST_PP_CAT(iter,BOOST_PP_DEC(BOOST_PP_FRAME_ITERATION(2)))
#   define AUX_ITER_1 BOOST_PP_CAT(iter,BOOST_PP_FRAME_ITERATION(2))

        typedef typename BOOST_MPL_AUX_NEXT(AUX_ITER_0) AUX_ITER_1;
        
#   undef AUX_ITER_1
#   undef AUX_ITER_0

#endif // BOOST_PP_IS_ITERATING
