// preprocessed version of 'boost/mpl/aux_/fold_impl.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {
namespace aux {

// forward declaration
template<
      nttp_long N
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct fold_impl;

template< nttp_long N >
struct fold_chunk;

template<>
struct fold_chunk<0>
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
        typedef state0 state;
        typedef iter0 iterator;
    };

    // ETI workaround
    template<> struct result_<int, int, int, int>
    {
        typedef int state;
        typedef int iterator;
    };

};

template<>
struct fold_chunk<1>
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
        typedef typename apply2<ForwardOp, state0, typename iter0::type>::type state1;
        typedef typename iter0::next iter1;
        

        typedef state1 state;
        typedef iter1 iterator;
    };

    // ETI workaround
    template<> struct result_<int, int, int, int>
    {
        typedef int state;
        typedef int iterator;
    };

};

template<>
struct fold_chunk<2>
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
        typedef typename apply2<ForwardOp, state0, typename iter0::type>::type state1;
        typedef typename iter0::next iter1;
        typedef typename apply2<ForwardOp, state1, typename iter1::type>::type state2;
        typedef typename iter1::next iter2;
        

        typedef state2 state;
        typedef iter2 iterator;
    };

    // ETI workaround
    template<> struct result_<int, int, int, int>
    {
        typedef int state;
        typedef int iterator;
    };

};

template<>
struct fold_chunk<3>
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
        typedef typename apply2<ForwardOp, state0, typename iter0::type>::type state1;
        typedef typename iter0::next iter1;
        typedef typename apply2<ForwardOp, state1, typename iter1::type>::type state2;
        typedef typename iter1::next iter2;
        typedef typename apply2<ForwardOp, state2, typename iter2::type>::type state3;
        typedef typename iter2::next iter3;
        

        typedef state3 state;
        typedef iter3 iterator;
    };

    // ETI workaround
    template<> struct result_<int, int, int, int>
    {
        typedef int state;
        typedef int iterator;
    };

};

template<>
struct fold_chunk<4>
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
        typedef typename apply2<ForwardOp, state0, typename iter0::type>::type state1;
        typedef typename iter0::next iter1;
        typedef typename apply2<ForwardOp, state1, typename iter1::type>::type state2;
        typedef typename iter1::next iter2;
        typedef typename apply2<ForwardOp, state2, typename iter2::type>::type state3;
        typedef typename iter2::next iter3;
        typedef typename apply2<ForwardOp, state3, typename iter3::type>::type state4;
        typedef typename iter3::next iter4;
        

        typedef state4 state;
        typedef iter4 iterator;
    };

    // ETI workaround
    template<> struct result_<int, int, int, int>
    {
        typedef int state;
        typedef int iterator;
    };

};

template< nttp_long N > 
struct fold_chunk
{
    template<
          typename First
        , typename Last
        , typename State
        , typename ForwardOp
        > 
    struct result_
    {
        typedef fold_impl<
              4
            , First
            , Last
            , State
            , ForwardOp
            > chunk_;

        typedef fold_impl<
              ( (N - 4) < 0 ? 0 : N - 4 )
            , typename chunk_::iterator
            , Last
            , typename chunk_::state
            , ForwardOp
            > res_;

        typedef typename res_::state state;
        typedef typename res_::iterator iterator;
    };
};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct fold_step;

template<
      typename Last
    , typename State
    >
struct fold_null_step
{
    typedef Last iterator;
    typedef State state;
};

template<> 
struct fold_chunk< -1 >
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
              typename is_same< First,Last >::type
            , fold_null_step< Last,State >
            , fold_step< First,Last,State,ForwardOp >
            >::type res_;

        typedef typename res_::state state;
        typedef typename res_::iterator iterator;
    };

    // ETI workaround
    template<> struct result_<int, int, int, int>
    {
        typedef int state;
        typedef int iterator;
    };

};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct fold_step
{
    typedef fold_chunk< -1 >::template result_<
          typename First::next
        , Last
        , typename apply2<ForwardOp, State, typename First::type>::type
        , ForwardOp
        > chunk_;

    typedef typename chunk_::state state;
    typedef typename chunk_::iterator iterator;
};

template<
      nttp_long N
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct fold_impl
    : fold_chunk<N>
        ::template result_< First,Last,State,ForwardOp >
{
};

} // namespace aux
} // namespace mpl
} // namespace boost

