//-----------------------------------------------------------------------------
// boost mpl/aux_/fold_backward_impl_body.hpp header file
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

// no include guards, the header is intended for multiple inclusion!

#if !defined(BOOST_PP_IS_ITERATING)

#   include "boost/mpl/aux_/apply.hpp"
#   include "boost/mpl/aux_/next.hpp"
#   include "boost/mpl/aux_/config/ctps.hpp"
#   include "boost/mpl/aux_/config/nttp.hpp"

#   include "boost/mpl/limits/unrolling.hpp"
#   include "boost/mpl/aux_/preprocessor/repeat.hpp"

#   include "boost/preprocessor/arithmetic/sub.hpp"
#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/dec.hpp"
#   include "boost/preprocessor/inc.hpp"
#   include "boost/preprocessor/cat.hpp"

// local macros, #undef-ined at the end of the header

#   define AUX_ITER_FOLD_FORWARD_STEP(unused, i, unused2) \
    typedef typename BOOST_MPL_AUX_APPLY2( \
          ForwardOp \
        , BOOST_PP_CAT(fwd_state,i) \
        , BOOST_MPL_AUX_FOLD_IMPL_OP(BOOST_PP_CAT(iter,i)) \
        )::type BOOST_PP_CAT(fwd_state,BOOST_PP_INC(i)); \
    typedef typename BOOST_MPL_AUX_NEXT(BOOST_PP_CAT(iter,i)) \
        BOOST_PP_CAT(iter,BOOST_PP_INC(i)); \
    /**/

#   define AUX_ITER_FOLD_BACKWARD_STEP_FUNC(i) \
    typedef typename BOOST_MPL_AUX_APPLY2( \
          BackwardOp \
        , BOOST_PP_CAT(bkwd_state,i) \
        , BOOST_MPL_AUX_FOLD_IMPL_OP(BOOST_PP_CAT(iter,BOOST_PP_DEC(i))) \
        )::type BOOST_PP_CAT(bkwd_state,BOOST_PP_DEC(i)); \
    /**/

#   define AUX_ITER_FOLD_BACKWARD_STEP(unused, i, j) \
    AUX_ITER_FOLD_BACKWARD_STEP_FUNC( \
        BOOST_PP_SUB_D(1,j,i) \
        ) \
    /**/

#   define AUX_FIRST_BACKWARD_STATE_TYPEDEF(i) \
    typedef typename nested_chunk::state BOOST_PP_CAT(bkwd_state,i);
    /**/

#   define AUX_FOLD_IMPL_NAME \
    BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_impl) \
    /**/

#   define AUX_FOLD_CHUNK_NAME \
    BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_chunk) \
    /**/

namespace boost {
namespace mpl {
namespace aux {

//: forward declaration
template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
    , typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME;

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
 && !defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION)

#   define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_UNROLLING_LIMIT, "boost/mpl/aux_/fold_backward_impl_body.hpp"))
#   include BOOST_PP_ITERATE()

// implementation for N that exceeds BOOST_MPL_UNROLLING_LIMIT
template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
    , typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME
{
    typedef First iter0;
    typedef State fwd_state0;

    BOOST_MPL_PP_REPEAT(
          BOOST_MPL_UNROLLING_LIMIT
        , AUX_ITER_FOLD_FORWARD_STEP
        , unused
        )

    typedef AUX_FOLD_IMPL_NAME<
          ( (N - BOOST_MPL_UNROLLING_LIMIT) < 0 ? 0 : N - BOOST_MPL_UNROLLING_LIMIT )
        , BOOST_PP_CAT(iter,BOOST_MPL_UNROLLING_LIMIT)
        , Last
        , BOOST_PP_CAT(fwd_state,BOOST_MPL_UNROLLING_LIMIT)
        , BackwardOp
        , ForwardOp
        > nested_chunk;
        
    AUX_FIRST_BACKWARD_STATE_TYPEDEF(BOOST_MPL_UNROLLING_LIMIT)

    BOOST_MPL_PP_REPEAT(
          BOOST_MPL_UNROLLING_LIMIT
        , AUX_ITER_FOLD_BACKWARD_STEP
        , BOOST_MPL_UNROLLING_LIMIT
        )

    typedef bkwd_state0 state;
    typedef typename nested_chunk::iterator iterator;
};

// fallback implementation for sequences of unknown size
template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME<-1,First,Last,State,BackwardOp,ForwardOp>
{
    typedef AUX_FOLD_IMPL_NAME<
          -1
        , typename BOOST_MPL_AUX_NEXT(First)
        , Last
        , typename BOOST_MPL_AUX_APPLY2(ForwardOp,State,BOOST_MPL_AUX_FOLD_IMPL_OP(First))::type
        , BackwardOp
        , ForwardOp
        > nested_step;

    typedef typename BOOST_MPL_AUX_APPLY2(
          BackwardOp
        , typename nested_step::state
        , BOOST_MPL_AUX_FOLD_IMPL_OP(First)
        )::type state;

    typedef typename nested_step::iterator iterator;
};

template<
      typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME<-1,Last,Last,State,BackwardOp,ForwardOp>
{
    typedef State state;
    typedef Last iterator;
};

#else // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

template< BOOST_MPL_AUX_NTTP_DECL(long, N) >
struct AUX_FOLD_CHUNK_NAME;

#   define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_UNROLLING_LIMIT, "boost/mpl/aux_/fold_backward_impl_body.hpp"))
#   include BOOST_PP_ITERATE()

// implementation for N that exceeds BOOST_MPL_UNROLLING_LIMIT
template< BOOST_MPL_AUX_NTTP_DECL(long, N) > 
struct AUX_FOLD_CHUNK_NAME
{
    template<
          typename First
        , typename Last
        , typename State
        , typename BackwardOp
        , typename ForwardOp
        > 
    struct result_
    {
        typedef First iter0;
        typedef State fwd_state0;

        BOOST_MPL_PP_REPEAT(
              BOOST_MPL_UNROLLING_LIMIT
            , AUX_ITER_FOLD_FORWARD_STEP
            , unused
            )

        typedef AUX_FOLD_IMPL_NAME<
              ( (N - BOOST_MPL_UNROLLING_LIMIT) < 0 ? 0 : N - BOOST_MPL_UNROLLING_LIMIT )
            , BOOST_PP_CAT(iter,BOOST_MPL_UNROLLING_LIMIT)
            , Last
            , BOOST_PP_CAT(fwd_state,BOOST_MPL_UNROLLING_LIMIT)
            , BackwardOp
            , ForwardOp
            > nested_chunk;
            
        AUX_FIRST_BACKWARD_STATE_TYPEDEF(BOOST_MPL_UNROLLING_LIMIT)

        BOOST_MPL_PP_REPEAT(
              BOOST_MPL_UNROLLING_LIMIT
            , AUX_ITER_FOLD_BACKWARD_STEP
            , BOOST_MPL_UNROLLING_LIMIT
            )

        typedef bkwd_state0 state;
        typedef typename nested_chunk::iterator iterator;
    };
};

// fallback implementation for sequences of unknown size
template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_step);

template<
      typename Last
    , typename State
    >
struct BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_null_step)
{
    typedef Last iterator;
    typedef State state;
};

template<> 
struct AUX_FOLD_CHUNK_NAME<-1>
{
    template<
          typename First
        , typename Last
        , typename State
        , typename BackwardOp
        , typename ForwardOp
        > 
    struct result_
    {
        typedef typename if_<
              typename is_same<First,Last>::type
            , BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_null_step)<Last,State>
            , BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_step)<First,Last,State,BackwardOp,ForwardOp>
            >::type res_;

        typedef typename res_::state state;
        typedef typename res_::iterator iterator;
    };

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
    //: ETI workaround
    template<> struct result_<int,int,int,int,int>
    {
        typedef int state;
        typedef int iterator;
    };
#endif
};

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_step)
{
    typedef AUX_FOLD_CHUNK_NAME<-1>::template result_<
          typename BOOST_MPL_AUX_NEXT(First)
        , Last
        , typename BOOST_MPL_AUX_APPLY2(ForwardOp,State,BOOST_MPL_AUX_FOLD_IMPL_OP(First))::type
        , BackwardOp
        , ForwardOp
        > nested_step;

    typedef typename BOOST_MPL_AUX_APPLY2(
          BackwardOp
        , typename nested_step::state
        , BOOST_MPL_AUX_FOLD_IMPL_OP(First)
        )::type state;

    typedef typename nested_step::iterator iterator;
};

template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
    , typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME
    : AUX_FOLD_CHUNK_NAME<N>
        ::template result_<First,Last,State,BackwardOp,ForwardOp>
{
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

} // namespace aux
} // namespace mpl
} // namespace boost

#   undef AUX_FIRST_BACKWARD_STATE_TYPEDEF
#   undef AUX_ITER_FOLD_BACKWARD_STEP
#   undef AUX_ITER_FOLD_BACKWARD_STEP_FUNC
#   undef AUX_ITER_FOLD_FORWARD_STEP

#undef BOOST_MPL_AUX_FOLD_IMPL_OP
#undef BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX

///// iteration

#else
#define i BOOST_PP_FRAME_ITERATION(1)

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION) \
 && !defined(BOOST_NO_NON_TYPE_TEMPLATE_PARTIAL_SPECIALIZATION)

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    >
struct AUX_FOLD_IMPL_NAME<i,First,Last,State,BackwardOp,ForwardOp>
{
    typedef First iter0;
    typedef State fwd_state0;

    BOOST_MPL_PP_REPEAT(
          i
        , AUX_ITER_FOLD_FORWARD_STEP
        , unused
        )

    typedef BOOST_PP_CAT(fwd_state,i) BOOST_PP_CAT(bkwd_state,i);

    BOOST_MPL_PP_REPEAT(
          i
        , AUX_ITER_FOLD_BACKWARD_STEP
        , i
        )

    typedef bkwd_state0 state;
    typedef BOOST_PP_CAT(iter,i) iterator;
};

#else

template<>
struct AUX_FOLD_CHUNK_NAME<BOOST_PP_FRAME_ITERATION(1)>
{
    template<
          typename First
        , typename Last
        , typename State
        , typename BackwardOp
        , typename ForwardOp
        >
    struct result_
    {
        typedef First iter0;
        typedef State fwd_state0;

        BOOST_MPL_PP_REPEAT(
              i
            , AUX_ITER_FOLD_FORWARD_STEP
            , unused
            )

        typedef BOOST_PP_CAT(fwd_state,i) BOOST_PP_CAT(bkwd_state,i);

        BOOST_MPL_PP_REPEAT(
              i
            , AUX_ITER_FOLD_BACKWARD_STEP
            , i
            )

        typedef bkwd_state0 state;
        typedef BOOST_PP_CAT(iter,i) iterator;
    };

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
    //: ETI workaround
    template<> struct result_<int,int,int,int,int>
    {
        typedef int state;
        typedef int iterator;
    };
#endif
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#undef i
#endif // BOOST_PP_IS_ITERATING
