// preprocessed version of 'boost/mpl/bind.hpp' header
// see the original for copyright information

namespace boost {
namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(bind1st)
BOOST_MPL_AUX_COMMON_NAME_WKND(bind2nd)

namespace aux {
template< bool >
struct resolve_arg_impl
{
    template<
          typename T, typename U1, typename U2, typename U3
        , typename U4, typename U5
        >
    struct result_
    {
        typedef T type;
    };
};

template<>
struct resolve_arg_impl<true>
{
    template<
          typename T, typename U1, typename U2, typename U3
        , typename U4, typename U5
        >
    struct result_
    {
        typedef typename apply5< T,U1,U2,U3,U4,U5 >::type type;
    };
};

template< typename T > struct is_bind_template;

template<
      typename T, typename U1, typename U2, typename U3, typename U4
    , typename U5
    >
struct resolve_bind_arg
    : resolve_arg_impl< is_bind_template<T >::value >
            ::template result_< T,U1,U2,U3,U4,U5 >
{
};

} // namespace aux
template< typename F, typename T > struct bind1st;
template< typename F, typename T > struct bind2nd;

namespace aux {
template< nttp_int arity_ > struct bind_impl_chooser;

aux::no_tag is_bind_helper(...);
template< typename T > aux::no_tag is_bind_helper(protect<T>*);

template< nttp_int N >
aux::yes_tag is_bind_helper(arg<N>*);

template< typename F, typename T > aux::yes_tag is_bind_helper(bind1st< F,T >*);
template< typename F, typename T > aux::yes_tag is_bind_helper(bind2nd< F,T >*);

template< bool is_ref_ = true >
struct is_bind_template_impl
{
    template< typename T > struct result_
    {
        enum { value = false };
    };
};

template<>
struct is_bind_template_impl<false>
{
    template< typename T > struct result_
    {
        enum { value =
             sizeof(aux::is_bind_helper(static_cast<T*>(0))) ==
             sizeof(aux::yes_tag)
            };

    };
};

template< typename T > struct is_bind_template
    : is_bind_template_impl< ::boost::detail::is_reference_impl<T>::value >
        ::template result_<T>
{
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
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef typename aux::resolve_bind_arg< F,U1,U2,U3,U4,U5 >::type f_;

     public:
        typedef typename apply0<f_>::type type;
    };
};

namespace aux {

template<
      typename F
    >
aux::yes_tag
is_bind_helper(bind0<F>*);

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(1, bind0)

template<
      typename F, typename T1
    >
struct bind1
{
    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef typename aux::resolve_bind_arg< F,U1,U2,U3,U4,U5 >::type f_;
        typedef typename aux::resolve_bind_arg< T1,U1,U2,U3,U4,U5 >::type t1;

     public:
        typedef typename apply1< f_,t1 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1
    >
aux::yes_tag
is_bind_helper(bind1< F,T1 >*);

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(2, bind1)

template<
      typename F, typename T1, typename T2
    >
struct bind2
{
    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef typename aux::resolve_bind_arg< F,U1,U2,U3,U4,U5 >::type f_;
        typedef typename aux::resolve_bind_arg< T1,U1,U2,U3,U4,U5 >::type t1;
        typedef typename aux::resolve_bind_arg< T2,U1,U2,U3,U4,U5 >::type t2;

     public:
        typedef typename apply2< f_,t1,t2 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2
    >
aux::yes_tag
is_bind_helper(bind2< F,T1,T2 >*);

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(3, bind2)

template<
      typename F, typename T1, typename T2, typename T3
    >
struct bind3
{
    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef typename aux::resolve_bind_arg< F,U1,U2,U3,U4,U5 >::type f_;
        typedef typename aux::resolve_bind_arg< T1,U1,U2,U3,U4,U5 >::type t1;
        typedef typename aux::resolve_bind_arg< T2,U1,U2,U3,U4,U5 >::type t2;
        typedef typename aux::resolve_bind_arg< T3,U1,U2,U3,U4,U5 >::type t3;

     public:
        typedef typename apply3< f_,t1,t2,t3 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3
    >
aux::yes_tag
is_bind_helper(bind3< F,T1,T2,T3 >*);

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(4, bind3)

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
struct bind4
{
    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef typename aux::resolve_bind_arg< F,U1,U2,U3,U4,U5 >::type f_;
        typedef typename aux::resolve_bind_arg< T1,U1,U2,U3,U4,U5 >::type t1;
        typedef typename aux::resolve_bind_arg< T2,U1,U2,U3,U4,U5 >::type t2;
        typedef typename aux::resolve_bind_arg< T3,U1,U2,U3,U4,U5 >::type t3;
        typedef typename aux::resolve_bind_arg< T4,U1,U2,U3,U4,U5 >::type t4;

     public:
        typedef typename apply4< f_,t1,t2,t3,t4 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    >
aux::yes_tag
is_bind_helper(bind4< F,T1,T2,T3,T4 >*);

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(5, bind4)

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
struct bind5
{
    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef typename aux::resolve_bind_arg< F,U1,U2,U3,U4,U5 >::type f_;
        typedef typename aux::resolve_bind_arg< T1,U1,U2,U3,U4,U5 >::type t1;
        typedef typename aux::resolve_bind_arg< T2,U1,U2,U3,U4,U5 >::type t2;
        typedef typename aux::resolve_bind_arg< T3,U1,U2,U3,U4,U5 >::type t3;
        typedef typename aux::resolve_bind_arg< T4,U1,U2,U3,U4,U5 >::type t4;
        typedef typename aux::resolve_bind_arg< T5,U1,U2,U3,U4,U5 >::type t5;

     public:
        typedef typename apply5< f_,t1,t2,t3,t4,t5 >::type type;
    };
};

namespace aux {

template<
      typename F, typename T1, typename T2, typename T3, typename T4
    , typename T5
    >
aux::yes_tag
is_bind_helper(bind5< F,T1,T2,T3,T4,T5 >*);

} // namespace aux

BOOST_MPL_AUX_ARITY_SPEC(6, bind5)

template< typename F, typename T >
struct bind1st
{
    template<
          typename U
        , typename U2 = void_, typename U3 = void_, typename U4 = void_
        , typename U5 = void_
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
        , typename U2 = void_, typename U3 = void_, typename U4 = void_
        , typename U5 = void_
        >
    struct apply
        : apply2< F,U,T >
    {
    };
};

} // namespace mpl
} // namespace boost

