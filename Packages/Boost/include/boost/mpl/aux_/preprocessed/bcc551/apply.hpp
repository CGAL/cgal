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
    static int const arity = 1; typedef F arg1;
 friend class apply0_rebind;
 typedef apply0_rebind rebind;
 };
 class apply0_rebind { public: template< typename U1 > struct apply { typedef typename apply0<U1>::type type;
 };
 
};

template<
      typename F
    >
struct apply< F,void_,void_,void_,void_,void_ >
    : apply0<F>
{
};

namespace aux {
template<
      int N, typename F, typename T1
    >
struct apply_impl1;
}

namespace aux {

template<
      typename F, typename T1
    >
struct apply_impl1<
          1
        , F
        , T1
        >
{
    typedef typename F::template apply<
          T1
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1
    >
struct apply_impl1<
          2
        , F
        , T1
        >
{
    typedef typename F::template apply<
          T1
        , void_
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1
    >
struct apply_impl1<
          3
        , F
        , T1
        >
{
    typedef typename F::template apply<
          T1
        , void_, void_
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1
    >
struct apply_impl1<
          4
        , F
        , T1
        >
{
    typedef typename F::template apply<
          T1
        , void_, void_, void_
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1
    >
struct apply_impl1<
          5
        , F
        , T1
        >
{
    typedef typename F::template apply<
          T1
        , void_, void_, void_, void_
        > type;
};

} // namespace aux

template<
      typename F, typename T1
    >
struct apply1
    : aux::apply_impl1<
          ::boost::mpl::aux::arity< F,1 >::value
        , F
        , T1
        >::type
{
    static int const arity = 2; typedef F arg1;
 typedef T1 arg2;
 friend class apply1_rebind;
 typedef apply1_rebind rebind;
 };
 class apply1_rebind { public: template< typename U1, typename U2 > struct apply { typedef typename apply1< U1,U2 >::type type;
 };
 
};

template<
      typename F, typename T1
    >
struct apply< F,T1,void_,void_,void_,void_ >
    : apply1< F,T1 >
{
};

namespace aux {
template<
      int N, typename F, typename T1, typename T2
    >
struct apply_impl2;
}

namespace aux {

template<
      typename F, typename T1, typename T2
    >
struct apply_impl2<
          2
        , F
        , T1, T2
        >
{
    typedef typename F::template apply<
          T1, T2
         
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1, typename T2
    >
struct apply_impl2<
          3
        , F
        , T1, T2
        >
{
    typedef typename F::template apply<
          T1, T2
        , void_
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1, typename T2
    >
struct apply_impl2<
          4
        , F
        , T1, T2
        >
{
    typedef typename F::template apply<
          T1, T2
        , void_, void_
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1, typename T2
    >
struct apply_impl2<
          5
        , F
        , T1, T2
        >
{
    typedef typename F::template apply<
          T1, T2
        , void_, void_, void_
        > type;
};

} // namespace aux

template<
      typename F, typename T1, typename T2
    >
struct apply2
    : aux::apply_impl2<
          ::boost::mpl::aux::arity< F,2 >::value
        , F
        , T1, T2
        >::type
{
    static int const arity = 3; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 friend class apply2_rebind;
 typedef apply2_rebind rebind;
 };
 class apply2_rebind { public: template< typename U1, typename U2, typename U3 > struct apply { typedef typename apply2< U1,U2,U3 >::type type;
 };
 
};

template<
      typename F, typename T1, typename T2
    >
struct apply< F,T1,T2,void_,void_,void_ >
    : apply2< F,T1,T2 >
{
};

namespace aux {
template<
      int N, typename F, typename T1, typename T2, typename T3
    >
struct apply_impl3;
}

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply_impl3<
          3
        , F
        , T1, T2, T3
        >
{
    typedef typename F::template apply<
          T1, T2, T3
         
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply_impl3<
          4
        , F
        , T1, T2, T3
        >
{
    typedef typename F::template apply<
          T1, T2, T3
        , void_
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply_impl3<
          5
        , F
        , T1, T2, T3
        >
{
    typedef typename F::template apply<
          T1, T2, T3
        , void_, void_
        > type;
};

} // namespace aux

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply3
    : aux::apply_impl3<
          ::boost::mpl::aux::arity< F,3 >::value
        , F
        , T1, T2, T3
        >::type
{
    static int const arity = 4; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 typedef T3 arg4;
 friend class apply3_rebind;
 typedef apply3_rebind rebind;
 };
 class apply3_rebind { public: template< typename U1, typename U2, typename U3, typename U4 > struct apply { typedef typename apply3< U1,U2,U3,U4 >::type type;
 };
 
};

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply< F,T1,T2,T3,void_,void_ >
    : apply3< F,T1,T2,T3 >
{
};

namespace aux {
template<
      int N, typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply_impl4;
}

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply_impl4<
          4
        , F
        , T1, T2, T3, T4
        >
{
    typedef typename F::template apply<
          T1, T2, T3, T4
         
        > type;
};

} // namespace aux

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply_impl4<
          5
        , F
        , T1, T2, T3, T4
        >
{
    typedef typename F::template apply<
          T1, T2, T3, T4
        , void_
        > type;
};

} // namespace aux

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply4
    : aux::apply_impl4<
          ::boost::mpl::aux::arity< F,4 >::value
        , F
        , T1, T2, T3, T4
        >::type
{
    static int const arity = 5; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 typedef T3 arg4;
 typedef T4 arg5;
 friend class apply4_rebind;
 typedef apply4_rebind rebind;
 };
 class apply4_rebind { public: template< typename U1, typename U2, typename U3, typename U4, typename U5 > struct apply { typedef typename apply4< U1,U2,U3,U4,U5 >::type type;
 };
 
};

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply< F,T1,T2,T3,T4,void_ >
    : apply4< F,T1,T2,T3,T4 >
{
};

namespace aux {
template<
      int N, typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct apply_impl5;
}

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct apply_impl5<
          5
        , F
        , T1, T2, T3, T4, T5
        >
{
    typedef typename F::template apply<
          T1, T2, T3, T4, T5
         
        > type;
};

} // namespace aux

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct apply5
    : aux::apply_impl5<
          ::boost::mpl::aux::arity< F,5 >::value
        , F
        , T1, T2, T3, T4, T5
        >::type
{
    static int const arity = 6; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 typedef T3 arg4;
 typedef T4 arg5;
 typedef T5 arg6;
 friend class apply5_rebind;
 typedef apply5_rebind rebind;
 };
 class apply5_rebind { public: template< typename U1, typename U2, typename U3, typename U4, typename U5, typename U6 > struct apply { typedef typename apply5< U1,U2,U3,U4,U5,U6 >::type type;
 };
 
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

