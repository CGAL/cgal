// preprocessed version of 'boost/mpl/apply.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template<
      typename F, typename T1 = void_, typename T2 = void_
    , typename T3 = void_, typename T4 = void_, typename T5 = void_
    >
struct apply;

template< typename F >
struct apply0 : F
{
};

template<
      typename F
    >
struct apply< F,void_,void_,void_,void_,void_ >
    : apply0<F>
{
};

template<
      typename F, typename T1
    >
struct apply1
    : F::template apply<
          T1
        >
{
};

template<
      typename F, typename T1
    >
struct apply< F,T1,void_,void_,void_,void_ >
    : apply1< F,T1 >
{
};

template<
      typename F, typename T1, typename T2
    >
struct apply2
    : F::template apply<
          T1, T2
        >
{
};

template<
      typename F, typename T1, typename T2
    >
struct apply< F,T1,T2,void_,void_,void_ >
    : apply2< F,T1,T2 >
{
};

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply3
    : F::template apply<
          T1, T2, T3
        >
{
};

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply< F,T1,T2,T3,void_,void_ >
    : apply3< F,T1,T2,T3 >
{
};

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply4
    : F::template apply<
          T1, T2, T3, T4
        >
{
};

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply< F,T1,T2,T3,T4,void_ >
    : apply4< F,T1,T2,T3,T4 >
{
};

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct apply5
    : F::template apply<
          T1, T2, T3, T4, T5
        >
{
};

// primary template (not a specialization!)
template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct apply
    : apply5< F,T1,T2,T3,T4,T5 >
{
};

} // namespace mpl
} // namespace boost

