//-----------------------------------------------------------------------------
// boost mpl/aux_/iter_fold_if_impl.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy, David Abrahams
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_AUX_ITER_FOLD_IF_IMPL_HPP_INCLUDED
#define BOOST_MPL_AUX_ITER_FOLD_IF_IMPL_HPP_INCLUDED

#include "boost/mpl/aux_/apply.hpp"

#if !defined(BOOST_MPL_PREPROCESSING_MODE)
#   include "boost/mpl/identity.hpp"
#   include "boost/mpl/next.hpp"
#   include "boost/mpl/if.hpp"
#   include "boost/mpl/aux_/value_wknd.hpp"
#   include "boost/mpl/aux_/bool_value_wknd.hpp"
#endif

#include "boost/mpl/aux_/config/use_preprocessed.hpp"

#if !defined(BOOST_MPL_NO_PREPROCESSED_HEADERS) && \
    !defined(BOOST_MPL_PREPROCESSING_MODE)

#   define BOOST_MPL_PREPROCESSED_HEADER iter_fold_if_impl.hpp
#   include "boost/mpl/aux_/include_preprocessed.hpp"

#else

#   include "boost/mpl/limits/unrolling.hpp"
#   include "boost/preprocessor/arithmetic/sub.hpp"
#   include "boost/preprocessor/repeat.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/dec.hpp"
#   include "boost/preprocessor/cat.hpp"

namespace boost {
namespace mpl {
namespace aux {

template< typename Iterator, typename State >
struct iter_fold_if_null_step
{
    typedef State state;
    typedef Iterator iterator;
};

template< bool >
struct iter_fold_if_step_impl
{
    template<
          typename Iterator
        , typename State
        , typename StateOp
        , typename IteratorOp
        >
    struct result_
    {
        typedef typename BOOST_MPL_AUX_APPLY2(StateOp,State,Iterator)::type state;
        typedef typename IteratorOp::type iterator;
    };
};

template<>
struct iter_fold_if_step_impl<false>
{
    template<
          typename Iterator
        , typename State
        , typename StateOp
        , typename IteratorOp
        >
    struct result_
    {
        typedef State state;
        typedef Iterator iterator;
    };
};

// agurt, 25/jun/02: MSVC 6.5 workaround, had to get rid of inheritance 
// here and in 'iter_fold_if_backward_step', because sometimes it interfered 
// with the "early template instantiation bug" in _really_ ugly ways
template<
      typename Iterator
    , typename State
    , typename ForwardOp
    , typename Predicate
    >
struct iter_fold_if_forward_step
{
    typedef typename BOOST_MPL_AUX_APPLY2(Predicate,State,Iterator)::type not_last;
    typedef typename iter_fold_if_step_impl<
          BOOST_MPL_AUX_MSVC_VALUE_WKND(not_last)::value
        >::template result_< Iterator,State,ForwardOp,mpl::next<Iterator> > impl_;

    typedef typename impl_::state state;
    typedef typename impl_::iterator iterator;
};

template<
      typename Iterator
    , typename State
    , typename BackwardOp
    , typename Predicate
    >
struct iter_fold_if_backward_step
{
    typedef typename BOOST_MPL_AUX_APPLY2(Predicate,State,Iterator)::type not_last;
    typedef typename iter_fold_if_step_impl<
          BOOST_MPL_AUX_MSVC_VALUE_WKND(not_last)::value
        >::template result_< Iterator,State,BackwardOp,identity<Iterator> > impl_;

    typedef typename impl_::state state;
    typedef typename impl_::iterator iterator;
};


// local macros, #undef-ined at the end of the header

#   define AUX_ITER_FOLD_FORWARD_STEP(unused, i, unused2) \
    typedef iter_fold_if_forward_step< \
          typename BOOST_PP_CAT(forward_step,i)::iterator \
        , typename BOOST_PP_CAT(forward_step,i)::state \
        , ForwardOp \
        , ForwardPredicate \
        > BOOST_PP_CAT(forward_step, BOOST_PP_INC(i)); \
    /**/

#   define AUX_ITER_FOLD_BACKWARD_STEP_FUNC(i) \
    typedef iter_fold_if_backward_step< \
          typename BOOST_PP_CAT(forward_step,BOOST_PP_DEC(i))::iterator \
        , typename BOOST_PP_CAT(backward_step,i)::state \
        , BackwardOp \
        , BackwardPredicate \
        > BOOST_PP_CAT(backward_step,BOOST_PP_DEC(i)); \
    /**/

#   define AUX_ITER_FOLD_BACKWARD_STEP(unused, i, unused2) \
    AUX_ITER_FOLD_BACKWARD_STEP_FUNC( \
        BOOST_PP_SUB_D(1,BOOST_MPL_UNROLLING_LIMIT,i) \
        ) \
    /**/

#   define AUX_LAST_FORWARD_STEP \
    BOOST_PP_CAT(forward_step, BOOST_MPL_UNROLLING_LIMIT) \
    /**/

#   define AUX_LAST_BACKWARD_STEP \
    BOOST_PP_CAT(backward_step, BOOST_MPL_UNROLLING_LIMIT) \
    /**/

template<
      typename Iterator
    , typename State
    , typename ForwardOp
    , typename ForwardPredicate
    , typename BackwardOp
    , typename BackwardPredicate
    >
struct iter_fold_if_impl
{
 private:
    typedef iter_fold_if_null_step<Iterator,State> forward_step0;
    BOOST_PP_REPEAT_1(
          BOOST_MPL_UNROLLING_LIMIT
        , AUX_ITER_FOLD_FORWARD_STEP
        , unused
        )
    
    typedef typename if_<
          typename AUX_LAST_FORWARD_STEP::not_last
        , iter_fold_if_impl<
              typename AUX_LAST_FORWARD_STEP::iterator
            , typename AUX_LAST_FORWARD_STEP::state
            , ForwardOp
            , ForwardPredicate
            , BackwardOp
            , BackwardPredicate
            >
        , iter_fold_if_null_step<
              typename AUX_LAST_FORWARD_STEP::iterator
            , typename AUX_LAST_FORWARD_STEP::state
            >
        >::type AUX_LAST_BACKWARD_STEP;

    BOOST_PP_REPEAT_1(
          BOOST_MPL_UNROLLING_LIMIT
        , AUX_ITER_FOLD_BACKWARD_STEP
        , unused
        )

 public:
    typedef typename backward_step0::state state;
    typedef typename AUX_LAST_BACKWARD_STEP::iterator iterator;
};

#   undef AUX_LAST_BACKWARD_STEP
#   undef AUX_LAST_FORWARD_STEP
#   undef AUX_ITER_FOLD_BACKWARD_STEP
#   undef AUX_ITER_FOLD_BACKWARD_STEP_FUNC
#   undef AUX_ITER_FOLD_FORWARD_STEP

} // namespace aux
} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_USE_PREPROCESSED_HEADERS
#endif // BOOST_MPL_AUX_ITER_FOLD_IF_IMPL_HPP_INCLUDED
