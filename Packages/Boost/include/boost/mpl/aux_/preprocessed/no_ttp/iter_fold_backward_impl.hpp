// preprocessed version of 'boost/mpl/aux_/iter_fold_backward_impl.hpp' header
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
    , typename BackwardOp
    , typename ForwardOp
    > 
struct iter_fold_backward_impl;

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    >
struct iter_fold_backward_impl< 0,First,Last,State,BackwardOp,ForwardOp >
{
    typedef First iter0;
    typedef State fwd_state0;
    typedef fwd_state0 bkwd_state0;
    typedef bkwd_state0 state;
    typedef iter0 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    >
struct iter_fold_backward_impl< 1,First,Last,State,BackwardOp,ForwardOp >
{
    typedef First iter0;
    typedef State fwd_state0;
    typedef typename ForwardOp::template apply< fwd_state0,iter0 >::type fwd_state1;
    typedef typename iter0::next iter1;
    

    typedef fwd_state1 bkwd_state1;
    typedef typename BackwardOp::template apply< bkwd_state1,iter0 >::type bkwd_state0;
    typedef bkwd_state0 state;
    typedef iter1 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    >
struct iter_fold_backward_impl< 2,First,Last,State,BackwardOp,ForwardOp >
{
    typedef First iter0;
    typedef State fwd_state0;
    typedef typename ForwardOp::template apply< fwd_state0,iter0 >::type fwd_state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< fwd_state1,iter1 >::type fwd_state2;
    typedef typename iter1::next iter2;
    

    typedef fwd_state2 bkwd_state2;
    typedef typename BackwardOp::template apply< bkwd_state2,iter1 >::type bkwd_state1;
    typedef typename BackwardOp::template apply< bkwd_state1,iter0 >::type bkwd_state0;
    

    typedef bkwd_state0 state;
    typedef iter2 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    >
struct iter_fold_backward_impl< 3,First,Last,State,BackwardOp,ForwardOp >
{
    typedef First iter0;
    typedef State fwd_state0;
    typedef typename ForwardOp::template apply< fwd_state0,iter0 >::type fwd_state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< fwd_state1,iter1 >::type fwd_state2;
    typedef typename iter1::next iter2;
    typedef typename ForwardOp::template apply< fwd_state2,iter2 >::type fwd_state3;
    typedef typename iter2::next iter3;
    

    typedef fwd_state3 bkwd_state3;
    typedef typename BackwardOp::template apply< bkwd_state3,iter2 >::type bkwd_state2;
    typedef typename BackwardOp::template apply< bkwd_state2,iter1 >::type bkwd_state1;
    typedef typename BackwardOp::template apply< bkwd_state1,iter0 >::type bkwd_state0;
    

    typedef bkwd_state0 state;
    typedef iter3 iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    >
struct iter_fold_backward_impl< 4,First,Last,State,BackwardOp,ForwardOp >
{
    typedef First iter0;
    typedef State fwd_state0;
    typedef typename ForwardOp::template apply< fwd_state0,iter0 >::type fwd_state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< fwd_state1,iter1 >::type fwd_state2;
    typedef typename iter1::next iter2;
    typedef typename ForwardOp::template apply< fwd_state2,iter2 >::type fwd_state3;
    typedef typename iter2::next iter3;
    typedef typename ForwardOp::template apply< fwd_state3,iter3 >::type fwd_state4;
    typedef typename iter3::next iter4;
    

    typedef fwd_state4 bkwd_state4;
    typedef typename BackwardOp::template apply< bkwd_state4,iter3 >::type bkwd_state3;
    typedef typename BackwardOp::template apply< bkwd_state3,iter2 >::type bkwd_state2;
    typedef typename BackwardOp::template apply< bkwd_state2,iter1 >::type bkwd_state1;
    typedef typename BackwardOp::template apply< bkwd_state1,iter0 >::type bkwd_state0;
    

    typedef bkwd_state0 state;
    typedef iter4 iterator;
};

template<
      long N
    , typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct iter_fold_backward_impl
{
    typedef First iter0;
    typedef State fwd_state0;
    typedef typename ForwardOp::template apply< fwd_state0,iter0 >::type fwd_state1;
    typedef typename iter0::next iter1;
    typedef typename ForwardOp::template apply< fwd_state1,iter1 >::type fwd_state2;
    typedef typename iter1::next iter2;
    typedef typename ForwardOp::template apply< fwd_state2,iter2 >::type fwd_state3;
    typedef typename iter2::next iter3;
    typedef typename ForwardOp::template apply< fwd_state3,iter3 >::type fwd_state4;
    typedef typename iter3::next iter4;
    

    typedef iter_fold_backward_impl<
          ( (N - 4) < 0 ? 0 : N - 4 )
        , iter4
        , Last
        , fwd_state4
        , BackwardOp
        , ForwardOp
        > nested_chunk;
        
    typedef typename nested_chunk::state bkwd_state4;
    typedef typename BackwardOp::template apply< bkwd_state4,iter3 >::type bkwd_state3;
    typedef typename BackwardOp::template apply< bkwd_state3,iter2 >::type bkwd_state2;
    typedef typename BackwardOp::template apply< bkwd_state2,iter1 >::type bkwd_state1;
    typedef typename BackwardOp::template apply< bkwd_state1,iter0 >::type bkwd_state0;
    

    typedef bkwd_state0 state;
    typedef typename nested_chunk::iterator iterator;
};

template<
      typename First
    , typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct iter_fold_backward_impl< -1,First,Last,State,BackwardOp,ForwardOp >
{
    typedef iter_fold_backward_impl<
          -1
        , typename First::next
        , Last
        , typename ForwardOp::template apply< State,First >::type
        , BackwardOp
        , ForwardOp
        > nested_step;

    typedef typename BackwardOp::template apply<typename nested_step::state, First>::type state;
    typedef typename nested_step::iterator iterator;
};

template<
      typename Last
    , typename State
    , typename BackwardOp
    , typename ForwardOp
    > 
struct iter_fold_backward_impl< -1,Last,Last,State,BackwardOp,ForwardOp >
{
    typedef State state;
    typedef Last iterator;
};

} // namespace aux
} // namespace mpl
} // namespace boost

