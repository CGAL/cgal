// preprocessed version of 'boost/mpl/aux_/iter_fold_impl.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {
namespace aux {

// forward declaration
template<
      long N
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct iter_fold_impl;

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    >
struct iter_fold_impl< 0,First,Last,State,ForwardOp >
{
    typedef First iter0;
    typedef State state0;
    typedef state0 state;
    typedef iter0 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    >
struct iter_fold_impl< 1,First,Last,State,ForwardOp >
{
    typedef First iter0;
    typedef State state0;
    typedef typename ForwardOp::template apply< state0,iter0 >::type state1;
    typedef typename iter0::next iter1;
    

    typedef state1 state;
    typedef iter1 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    >
struct iter_fold_impl< 2,First,Last,State,ForwardOp >
{
    typedef First iter0;
    typedef State state0;
    typedef typename ForwardOp::template apply< state0,iter0 >::type state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< state1,iter1 >::type state2;
    typedef typename iter1::next iter2;
    

    typedef state2 state;
    typedef iter2 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    >
struct iter_fold_impl< 3,First,Last,State,ForwardOp >
{
    typedef First iter0;
    typedef State state0;
    typedef typename ForwardOp::template apply< state0,iter0 >::type state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< state1,iter1 >::type state2;
    typedef typename iter1::next iter2;
    typedef typename ForwardOp::template apply< state2,iter2 >::type state3;
    typedef typename iter2::next iter3;
    

    typedef state3 state;
    typedef iter3 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    >
struct iter_fold_impl< 4,First,Last,State,ForwardOp >
{
    typedef First iter0;
    typedef State state0;
    typedef typename ForwardOp::template apply< state0,iter0 >::type state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< state1,iter1 >::type state2;
    typedef typename iter1::next iter2;
    typedef typename ForwardOp::template apply< state2,iter2 >::type state3;
    typedef typename iter2::next iter3;
    typedef typename ForwardOp::template apply< state3,iter3 >::type state4;
    typedef typename iter3::next iter4;
    

    typedef state4 state;
    typedef iter4 iterator;
};

template<
      long N
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct iter_fold_impl
{
    typedef iter_fold_impl<
          4
        , First
        , Last
        , State
        , ForwardOp
        > chunk_;

    typedef iter_fold_impl<
          ( (N - 4) < 0 ? 0 : N - 4 )
        , typename chunk_::iterator
        , Last
        , typename chunk_::state
        , ForwardOp
        > res_;
        
    typedef typename res_::state state;
    typedef typename res_::iterator iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct iter_fold_impl< -1,First,Last,State,ForwardOp >
    : iter_fold_impl<
          -1
        , typename First::next
        , Last
        , typename ForwardOp::template apply< State,First >::type
        , ForwardOp
        >
{
};

template<
      typename Last
    , typename State
    , typename ForwardOp
    > 
struct iter_fold_impl< -1,Last,Last,State,ForwardOp >
{
    typedef State state;
    typedef Last iterator;
};

} // namespace aux
} // namespace mpl
} // namespace boost

