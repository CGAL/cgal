// preprocessed version of 'boost/mpl/aux_/fold_impl.hpp' header
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
struct fold_impl;

template<
      long N
    , typename First
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct fold_impl
{
    typedef fold_impl<
          -1
        , typename First::next
        , Last
        , typename apply2<ForwardOp, State, typename First::type>::type
        , ForwardOp
        > res_;

    typedef typename res_::state state;
    typedef typename res_::iterator iterator;
    typedef state type;
};

template<
      long N
    , typename Last
    , typename State
    , typename ForwardOp
    > 
struct fold_impl< N,Last,Last,State,ForwardOp >
{
    typedef State state;
    typedef Last iterator;
    typedef state type;
};

} // namespace aux
} // namespace mpl
} // namespace boost

