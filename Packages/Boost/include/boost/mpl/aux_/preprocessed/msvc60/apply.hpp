// preprocessed version of 'boost/mpl/apply.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

template< typename F >
struct apply0 : F
{
    enum { arity = 1 }; typedef F arg1;
 friend class apply0_rebind;
 typedef apply0_rebind rebind;
 };
 class apply0_rebind { public: template< typename U1 > struct apply : apply0<U1> { };
 
};

// workaround for the ETI bug
template<>
struct apply0<int>
{
    typedef int type;
};

namespace aux {

template< typename F>
struct msvc_apply1
{
    template< bool > struct f_ : F {};
    template<> struct f_<true>
    {
        template< typename P1 > struct apply
        {
        };
    };

    template< typename T1 > struct result_
        : f_< aux::msvc_never_true<F>::value >
            ::template apply<T1>
    {
    };
};

} // namespace aux

template<
      typename F, typename T1
    >
struct apply1
{
    typedef typename aux::msvc_apply1<F>::template result_<
          T1
    >::type type;
    
    enum { arity = 2 }; typedef F arg1;
 typedef T1 arg2;
 friend class apply1_rebind;
 typedef apply1_rebind rebind;
 };
 class apply1_rebind { public: template< typename U1, typename U2 > struct apply : apply1< U1,U2 > { };
 
};

// workaround for ETI bug
template<>
struct apply1< int,int >
{
    typedef int type;
};

namespace aux {

template< typename F>
struct msvc_apply2
{
    template< bool > struct f_ : F {};
    template<> struct f_<true>
    {
        template< typename P1, typename P2 > struct apply
        {
        };
    };

    template< typename T1, typename T2 > struct result_
        : f_< aux::msvc_never_true<F>::value >
            ::template apply< T1,T2 >
    {
    };
};

} // namespace aux

template<
      typename F, typename T1, typename T2
    >
struct apply2
{
    typedef typename aux::msvc_apply2<F>::template result_<
          T1, T2
    >::type type;
    
    enum { arity = 3 }; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 friend class apply2_rebind;
 typedef apply2_rebind rebind;
 };
 class apply2_rebind { public: template< typename U1, typename U2, typename U3 > struct apply : apply2< U1,U2,U3 > { };
 
};

// workaround for ETI bug
template<>
struct apply2< int,int,int >
{
    typedef int type;
};

namespace aux {

template< typename F>
struct msvc_apply3
{
    template< bool > struct f_ : F {};
    template<> struct f_<true>
    {
        template< typename P1, typename P2, typename P3 > struct apply
        {
        };
    };

    template< typename T1, typename T2, typename T3 > struct result_
        : f_< aux::msvc_never_true<F>::value >
            ::template apply< T1,T2,T3 >
    {
    };
};

} // namespace aux

template<
      typename F, typename T1, typename T2, typename T3
    >
struct apply3
{
    typedef typename aux::msvc_apply3<F>::template result_<
          T1, T2, T3
    >::type type;
    enum { arity = 4 }; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 typedef T3 arg4;
 friend class apply3_rebind;
 typedef apply3_rebind rebind;
 };
 class apply3_rebind { public: template< typename U1, typename U2, typename U3, typename U4 > struct apply : apply3< U1,U2,U3,U4 > { };
 
};

// workaround for ETI bug
template<>
struct apply3< int,int,int,int >
{
    typedef int type;
};

namespace aux {

template< typename F>
struct msvc_apply4
{
    template< bool > struct f_ : F {};
    template<> struct f_<true>
    {
        template<
              typename P1, typename P2, typename P3, typename P4
            >
        struct apply
        {
        };
    };

    template<
          typename T1, typename T2, typename T3, typename T4
        >
    struct result_
        : f_< aux::msvc_never_true<F>::value >
            ::template apply< T1,T2,T3,T4 >
    {
    };
};

} // namespace aux

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct apply4
{
    typedef typename aux::msvc_apply4<F>::template result_<
        T1, T2, T3, T4
    >::type type;
    enum { arity = 5 }; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 typedef T3 arg4;
 typedef T4 arg5;
 friend class apply4_rebind;
 typedef apply4_rebind rebind;
 };
 class apply4_rebind { public: template< typename U1, typename U2, typename U3, typename U4, typename U5 > struct apply : apply4< U1,U2,U3,U4,U5 > { };
 
};

// workaround for ETI bug
template<>
struct apply4< int,int,int,int,int >
{
    typedef int type;
};

namespace aux {

template< typename F>
struct msvc_apply5
{
    template< bool > struct f_ : F {};
    template<> struct f_<true>
    {
        template<
              typename P1, typename P2, typename P3, typename P4
            , typename P5
            >
        struct apply
        {
        };
    };

    template<
          typename T1, typename T2, typename T3, typename T4
        , typename T5
        >
    struct result_
        : f_< aux::msvc_never_true<F>::value >
            ::template apply< T1,T2,T3,T4,T5 >
    {
    };
};

} // namespace aux

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct apply5
{
    typedef typename aux::msvc_apply5<F>::template result_<
        T1, T2, T3, T4, T5
    >::type type;
    enum { arity = 6 }; typedef F arg1;
 typedef T1 arg2;
 typedef T2 arg3;
 typedef T3 arg4;
 typedef T4 arg5;
 typedef T5 arg6;
 friend class apply5_rebind;
 typedef apply5_rebind rebind;
 };
 class apply5_rebind { public: template< typename U1, typename U2, typename U3, typename U4, typename U5, typename U6 > struct apply : apply5< U1,U2,U3,U4,U5,U6 > { };
 
};

// workaround for ETI bug
template<>
struct apply5< int,int,int,int,int,int >
{
    typedef int type;
};

} // namespace mpl
} // namespace boost

