// preprocessed version of 'boost/mpl/bind.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(bind1st)
BOOST_MPL_AUX_COMMON_NAME_WKND(bind2nd)

namespace aux {

template<
      typename T, typename U1, typename U2, typename U3, typename U4
    , typename U5
    >
struct resolve_bind_arg
{
    typedef T type;
};

template<
      typename T
    , typename Arg
    >
struct replace_unnamed_arg
{
    typedef Arg next_arg;
    typedef T type;
};

template<
      typename Arg
    >
struct replace_unnamed_arg< arg<-1>,Arg >
{
    typedef typename Arg::next next_arg;
    typedef Arg type;
};

} // namespace aux
template< typename F, typename T > struct bind1st;
template< typename F, typename T > struct bind2nd;

namespace aux {

template<
      int N, typename U1, typename U2, typename U3, typename U4, typename U5
    >
struct resolve_bind_arg< arg<N>,U1,U2,U3,U4,U5 >
{
    typedef typename apply5<mpl::arg< N>,U1,U2,U3,U4,U5 >::type type;
};

template<
      typename F, typename T, typename U1, typename U2, typename U3
    , typename U4, typename U5
    >
struct resolve_bind_arg< bind1st<F,T>,U1,U2,U3,U4,U5 >
{
    typedef bind1st< F,T > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

template<
      typename F, typename T, typename U1, typename U2, typename U3
    , typename U4, typename U5
    >
struct resolve_bind_arg< bind2nd<F,T>,U1,U2,U3,U4,U5 >
{
    typedef bind2nd< F,T > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux
BOOST_MPL_AUX_ARITY_SPEC(2, bind1st)
BOOST_MPL_AUX_ARITY_SPEC(2, bind2nd)

template<
      typename F
    >
struct bind0
{
    template<
          typename U1, typename U2, typename U3, typename U4, typename U5
        >
    struct apply
    {
     private:
        typedef aux::replace_unnamed_arg< F, mpl::arg< 1> > r0;
        typedef typename r0::type a0;
        typedef typename r0::next_arg n1;
        typedef typename aux::resolve_bind_arg< a0,U1,U2,U3,U4,U5 >::type f_;

     public:
        typedef typename apply0<f_>::type type;
    };
};

namespace aux {

template<
      typename F, typename U1, typename U2, typename U3, typename U4
    , typename U5
    >
struct resolve_bind_arg<
      bind0<F>,U1, U2, U3, U4, U5
    >
{
    typedef bind0<F> f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(1, bind0)

template<
      typename F, typename T1
    >
struct bind1
{
    template<
          typename U1, typename U2, typename U3, typename U4, typename U5
        >
    struct apply
    {
     private:
        typedef aux::replace_unnamed_arg< F, mpl::arg< 1> > r0;
        typedef typename r0::type a0;
        typedef typename r0::next_arg n1;
        typedef typename aux::resolve_bind_arg< a0,U1,U2,U3,U4,U5 >::type f_;

        typedef aux::replace_unnamed_arg< T1,n1 > r1;
        typedef typename r1::type a1;
        typedef typename r1::next_arg n2;
        typedef typename aux::resolve_bind_arg< a1,U1,U2,U3,U4,U5 >::type t1;

     public:
        typedef typename apply1< f_,t1 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename U1, typename U2, typename U3
    , typename U4, typename U5
    >
struct resolve_bind_arg<
      bind1< F,T1 >,U1, U2, U3, U4, U5
    >
{
    typedef bind1< F,T1 > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(2, bind1)

template<
      typename F, typename T1, typename T2
    >
struct bind2
{
    template<
          typename U1, typename U2, typename U3, typename U4, typename U5
        >
    struct apply
    {
     private:
        typedef aux::replace_unnamed_arg< F, mpl::arg< 1> > r0;
        typedef typename r0::type a0;
        typedef typename r0::next_arg n1;
        typedef typename aux::resolve_bind_arg< a0,U1,U2,U3,U4,U5 >::type f_;

        typedef aux::replace_unnamed_arg< T1,n1 > r1;
        typedef typename r1::type a1;
        typedef typename r1::next_arg n2;
        typedef typename aux::resolve_bind_arg< a1,U1,U2,U3,U4,U5 >::type t1;

        typedef aux::replace_unnamed_arg< T2,n2 > r2;
        typedef typename r2::type a2;
        typedef typename r2::next_arg n3;
        typedef typename aux::resolve_bind_arg< a2,U1,U2,U3,U4,U5 >::type t2;

     public:
        typedef typename apply2< f_,t1,t2 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename U1, typename U2
    , typename U3, typename U4, typename U5
    >
struct resolve_bind_arg<
      bind2< F,T1,T2 >,U1, U2, U3, U4, U5
    >
{
    typedef bind2< F,T1,T2 > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(3, bind2)

template<
      typename F, typename T1, typename T2, typename T3
    >
struct bind3
{
    template<
          typename U1, typename U2, typename U3, typename U4, typename U5
        >
    struct apply
    {
     private:
        typedef aux::replace_unnamed_arg< F, mpl::arg< 1> > r0;
        typedef typename r0::type a0;
        typedef typename r0::next_arg n1;
        typedef typename aux::resolve_bind_arg< a0,U1,U2,U3,U4,U5 >::type f_;

        typedef aux::replace_unnamed_arg< T1,n1 > r1;
        typedef typename r1::type a1;
        typedef typename r1::next_arg n2;
        typedef typename aux::resolve_bind_arg< a1,U1,U2,U3,U4,U5 >::type t1;

        typedef aux::replace_unnamed_arg< T2,n2 > r2;
        typedef typename r2::type a2;
        typedef typename r2::next_arg n3;
        typedef typename aux::resolve_bind_arg< a2,U1,U2,U3,U4,U5 >::type t2;

        typedef aux::replace_unnamed_arg< T3,n3 > r3;
        typedef typename r3::type a3;
        typedef typename r3::next_arg n4;
        typedef typename aux::resolve_bind_arg< a3,U1,U2,U3,U4,U5 >::type t3;

     public:
        typedef typename apply3< f_,t1,t2,t3 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename U1
    , typename U2, typename U3, typename U4, typename U5
    >
struct resolve_bind_arg<
      bind3< F,T1,T2,T3 >,U1, U2, U3, U4, U5
    >
{
    typedef bind3< F,T1,T2,T3 > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(4, bind3)

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct bind4
{
    template<
          typename U1, typename U2, typename U3, typename U4, typename U5
        >
    struct apply
    {
     private:
        typedef aux::replace_unnamed_arg< F, mpl::arg< 1> > r0;
        typedef typename r0::type a0;
        typedef typename r0::next_arg n1;
        typedef typename aux::resolve_bind_arg< a0,U1,U2,U3,U4,U5 >::type f_;

        typedef aux::replace_unnamed_arg< T1,n1 > r1;
        typedef typename r1::type a1;
        typedef typename r1::next_arg n2;
        typedef typename aux::resolve_bind_arg< a1,U1,U2,U3,U4,U5 >::type t1;

        typedef aux::replace_unnamed_arg< T2,n2 > r2;
        typedef typename r2::type a2;
        typedef typename r2::next_arg n3;
        typedef typename aux::resolve_bind_arg< a2,U1,U2,U3,U4,U5 >::type t2;

        typedef aux::replace_unnamed_arg< T3,n3 > r3;
        typedef typename r3::type a3;
        typedef typename r3::next_arg n4;
        typedef typename aux::resolve_bind_arg< a3,U1,U2,U3,U4,U5 >::type t3;

        typedef aux::replace_unnamed_arg< T4,n4 > r4;
        typedef typename r4::type a4;
        typedef typename r4::next_arg n5;
        typedef typename aux::resolve_bind_arg< a4,U1,U2,U3,U4,U5 >::type t4;

     public:
        typedef typename apply4< f_,t1,t2,t3,t4 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename U1, typename U2, typename U3, typename U4, typename U5
    >
struct resolve_bind_arg<
      bind4< F,T1,T2,T3,T4 >,U1, U2, U3, U4, U5
    >
{
    typedef bind4< F,T1,T2,T3,T4 > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(5, bind4)

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct bind5
{
    template<
          typename U1, typename U2, typename U3, typename U4, typename U5
        >
    struct apply
    {
     private:
        typedef aux::replace_unnamed_arg< F, mpl::arg< 1> > r0;
        typedef typename r0::type a0;
        typedef typename r0::next_arg n1;
        typedef typename aux::resolve_bind_arg< a0,U1,U2,U3,U4,U5 >::type f_;

        typedef aux::replace_unnamed_arg< T1,n1 > r1;
        typedef typename r1::type a1;
        typedef typename r1::next_arg n2;
        typedef typename aux::resolve_bind_arg< a1,U1,U2,U3,U4,U5 >::type t1;

        typedef aux::replace_unnamed_arg< T2,n2 > r2;
        typedef typename r2::type a2;
        typedef typename r2::next_arg n3;
        typedef typename aux::resolve_bind_arg< a2,U1,U2,U3,U4,U5 >::type t2;

        typedef aux::replace_unnamed_arg< T3,n3 > r3;
        typedef typename r3::type a3;
        typedef typename r3::next_arg n4;
        typedef typename aux::resolve_bind_arg< a3,U1,U2,U3,U4,U5 >::type t3;

        typedef aux::replace_unnamed_arg< T4,n4 > r4;
        typedef typename r4::type a4;
        typedef typename r4::next_arg n5;
        typedef typename aux::resolve_bind_arg< a4,U1,U2,U3,U4,U5 >::type t4;

        typedef aux::replace_unnamed_arg< T5,n5 > r5;
        typedef typename r5::type a5;
        typedef typename r5::next_arg n6;
        typedef typename aux::resolve_bind_arg< a5,U1,U2,U3,U4,U5 >::type t5;

     public:
        typedef typename apply5< f_,t1,t2,t3,t4,t5 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5, typename U1, typename U2, typename U3, typename U4
    , typename U5
    >
struct resolve_bind_arg<
      bind5< F,T1,T2,T3,T4,T5 >,U1, U2, U3, U4, U5
    >
{
    typedef bind5< F,T1,T2,T3,T4,T5 > f_;
    typedef typename apply5< f_,U1,U2,U3,U4,U5 >::type type;
};

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(6, bind5)

template< typename F, typename T >
struct bind1st
{
    template<
          typename U
        , typename U2, typename U3, typename U4, typename U5
        >
    struct apply
        : apply2< F,T,U >
    {
    };
};

template< typename F, typename T >
struct bind2nd
{
    template<
          typename U
        , typename U2, typename U3, typename U4, typename U5
        >
    struct apply
        : apply2< F,U,T >
    {
    };
};

} // namespace mpl
} // namespace boost

