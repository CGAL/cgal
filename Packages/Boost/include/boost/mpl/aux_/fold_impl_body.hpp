//-----------------------------------------------------------------------------
// boost mpl/aux_/fold_impl_body.hpp header file
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

#   include "boost/preprocessor/iterate.hpp"
#   include "boost/preprocessor/dec.hpp"
#   include "boost/preprocessor/cat.hpp"

// local macros, #undef-ined at the end of the header

#   define AUX_ITER_FOLD_STEP(unused, i, unused2) \
    typedef typename BOOST_MPL_AUX_APPLY2( \
          ForwardOp \
        , BOOST_PP_CAT(state,i) \
        , BOOST_MPL_AUX_FOLD_IMPL_OP(BOOST_PP_CAT(iter,i)) \
        )::type BOOST_PP_CAT(state,BOOST_PP_INC(i)); \
    typedef typename BOOST_MPL_AUX_NEXT(BOOST_PP_CAT(iter,i)) \
        BOOST_PP_CAT(iter,BOOST_PP_INC(i)); \
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
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME;

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#   if !defined(__BORLANDC__)

#   define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_UNROLLING_LIMIT, "boost/mpl/aux_/fold_impl_body.hpp"))
#   include BOOST_PP_ITERATE()

// implementation for N that exceeds BOOST_MPL_UNROLLING_LIMIT
template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME
{
    typedef AUX_FOLD_IMPL_NAME<
          BOOST_MPL_UNROLLING_LIMIT
        , First
        , Last
        , State
        , ForwardOp
        > chunk_;

    typedef AUX_FOLD_IMPL_NAME<
          ( (N - BOOST_MPL_UNROLLING_LIMIT) < 0 ? 0 : N - BOOST_MPL_UNROLLING_LIMIT )
        , typename chunk_::iterator
        , Last
        , typename chunk_::state
        , ForwardOp
        > res_;
        
    typedef typename res_::state state;
    typedef typename res_::iterator iterator;
};

// fallback implementation for sequences of unknown size
template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME<-1,First,Last,State,ForwardOp>
    : AUX_FOLD_IMPL_NAME<
          -1
        , typename BOOST_MPL_AUX_NEXT(First)
        , Last
        , typename BOOST_MPL_AUX_APPLY2(ForwardOp,State,BOOST_MPL_AUX_FOLD_IMPL_OP(First))::type
        , ForwardOp
        >
{
};

template<
      typename Last
    , typename State
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME<-1,Last,Last,State,ForwardOp>
{
    typedef State state;
    typedef Last iterator;
};

#   else // __BORLANDC__

// Borland have some serious problems with the unrolled version, so
// we always use a basic implementation
template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME
{
    typedef AUX_FOLD_IMPL_NAME<
          -1
        , typename First::next
        , Last
        , typename BOOST_MPL_AUX_APPLY2(ForwardOp,State,BOOST_MPL_AUX_FOLD_IMPL_OP(First))::type
        , ForwardOp
        > res_;

    typedef typename res_::state state;
    typedef typename res_::iterator iterator;
    typedef state type;
};

template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
     , typename Last
    , typename State
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME<N,Last,Last,State,ForwardOp >
{
    typedef State state;
    typedef Last iterator;
    typedef state type;
};

#   endif // __BORLANDC__
 
#else // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

template< BOOST_MPL_AUX_NTTP_DECL(long, N) >
struct AUX_FOLD_CHUNK_NAME;

#   define BOOST_PP_ITERATION_PARAMS_1 \
    (3,(0, BOOST_MPL_UNROLLING_LIMIT, "boost/mpl/aux_/fold_impl_body.hpp"))
#   include BOOST_PP_ITERATE()

// implementation for N that exceeds BOOST_MPL_UNROLLING_LIMIT
template< BOOST_MPL_AUX_NTTP_DECL(long, N) > 
struct AUX_FOLD_CHUNK_NAME
{
    template<
          typename First
        , typename Last
        , typename State
        , typename ForwardOp
        > 
    struct result_
    {
        typedef AUX_FOLD_IMPL_NAME<
              BOOST_MPL_UNROLLING_LIMIT
            , First
            , Last
            , State
            , ForwardOp
            > chunk_;

        typedef AUX_FOLD_IMPL_NAME<
              ( (N - BOOST_MPL_UNROLLING_LIMIT) < 0 ? 0 : N - BOOST_MPL_UNROLLING_LIMIT )
            , typename chunk_::iterator
            , Last
            , typename chunk_::state
            , ForwardOp
            > res_;

        typedef typename res_::state state;
        typedef typename res_::iterator iterator;
    };
};

// fallback implementation for sequences of unknown size
template<
      typename First
    , typename Last
    , typename State
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
        , typename ForwardOp
        > 
    struct result_
    {
        typedef typename if_<
              typename is_same<First,Last>::type
            , BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_null_step)<Last,State>
            , BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_step)<First,Last,State,ForwardOp>
            >::type res_;

        typedef typename res_::state state;
        typedef typename res_::iterator iterator;
    };

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
    //: ETI workaround
    template<> struct result_<int,int,int,int>
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
    , typename ForwardOp
    > 
struct BOOST_PP_CAT(BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX,_step)
{
    // can't inherit here - it breaks MSVC 7.0
    typedef AUX_FOLD_CHUNK_NAME<-1>::template result_<
          typename BOOST_MPL_AUX_NEXT(First)
        , Last
        , typename BOOST_MPL_AUX_APPLY2(ForwardOp,State,BOOST_MPL_AUX_FOLD_IMPL_OP(First))::type
        , ForwardOp
        > chunk_;

    typedef typename chunk_::state state;
    typedef typename chunk_::iterator iterator;
};

template<
      BOOST_MPL_AUX_NTTP_DECL(long, N)
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct AUX_FOLD_IMPL_NAME
    : AUX_FOLD_CHUNK_NAME<N>
        ::template result_<First,Last,State,ForwardOp>
{
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

} // namespace aux
} // namespace mpl
} // namespace boost

#   undef AUX_FOLD_IMPL_NAME
#   undef AUX_FOLD_CHUNK_NAME
#   undef AUX_ITER_FOLD_STEP

#undef BOOST_MPL_AUX_FOLD_IMPL_OP
#undef BOOST_MPL_AUX_FOLD_IMPL_NAME_PREFIX

///// iteration

#else

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    >
struct AUX_FOLD_IMPL_NAME<BOOST_PP_FRAME_ITERATION(1),First,Last,State,ForwardOp>
{
    typedef First iter0;
    typedef State state0;

    BOOST_MPL_PP_REPEAT(
          BOOST_PP_FRAME_ITERATION(1)
        , AUX_ITER_FOLD_STEP
        , unused
        )

    typedef BOOST_PP_CAT(state,BOOST_PP_FRAME_ITERATION(1)) state;
    typedef BOOST_PP_CAT(iter,BOOST_PP_FRAME_ITERATION(1)) iterator;
};

#else

template<>
struct AUX_FOLD_CHUNK_NAME<BOOST_PP_FRAME_ITERATION(1)>
{
    template<
          typename First
        , typename Last
        , typename State
        , typename ForwardOp
        >
    struct result_
    {
        typedef First iter0;
        typedef State state0;

        BOOST_MPL_PP_REPEAT(
              BOOST_PP_FRAME_ITERATION(1)
            , AUX_ITER_FOLD_STEP
            , unused
            )

        typedef BOOST_PP_CAT(state,BOOST_PP_FRAME_ITERATION(1)) state;
        typedef BOOST_PP_CAT(iter,BOOST_PP_FRAME_ITERATION(1)) iterator;
    };

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
    //: ETI workaround
    template<> struct result_<int,int,int,int>
    {
        typedef int state;
        typedef int iterator;
    };
#endif
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

#endif // BOOST_PP_IS_ITERATING
