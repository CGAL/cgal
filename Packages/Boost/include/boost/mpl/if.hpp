
#ifndef BOOST_MPL_IF_HPP_INCLUDED
#define BOOST_MPL_IF_HPP_INCLUDED

// + file: boost/mpl/if.hpp
// + last modified: 17/sep/03

// Copyright (c) 2000-03 Boost.org
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.
//
// See http://www.boost.org/libs/mpl for documentation.

#include "boost/mpl/void.hpp"
#include "boost/mpl/aux_/value_wknd.hpp"
#include "boost/mpl/aux_/static_cast.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"
#include "boost/config.hpp"

#if !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)
#   include "boost/mpl/arg_fwd.hpp"
#endif

namespace boost {
namespace mpl {

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      bool C
    , typename T1
    , typename T2
    >
struct if_c
{
    typedef T1 type;
};

template<
      typename T1
    , typename T2
    >
struct if_c<false,T1,T2>
{
    typedef T2 type;
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(C)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    >
struct if_
{
 private:
    // agurt, 02/jan/03: two-step 'type' definition for the sake of aCC 
    typedef if_c<
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x561))
          BOOST_MPL_AUX_VALUE_WKND(C)::value
#else
          BOOST_MPL_AUX_STATIC_CAST(bool, BOOST_MPL_AUX_VALUE_WKND(C)::value)
#endif
        , T1
        , T2
        > almost_type_;
 
 public:
    typedef typename almost_type_::type type;
    
    BOOST_MPL_AUX_LAMBDA_SUPPORT(3,if_,(C,T1,T2))
};

#else

// no partial class template specialization

namespace aux {

template< bool C >
struct if_impl
{
    template< typename T1, typename T2 > struct result_
    {
        typedef T1 type;
    };
};

template<>
struct if_impl<false>
{
    template< typename T1, typename T2 > struct result_
    { 
        typedef T2 type;
    };
};

} // namespace aux

template<
      bool C_
    , typename T1
    , typename T2
    >
struct if_c
{
    typedef typename aux::if_impl< C_ >
        ::template result_<T1,T2>::type type;
};

// (almost) copy & paste in order to save one more 
// recursively nested template instantiation to user
template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(C_)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T1)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T2)
    >
struct if_
{
    enum { msvc70_wknd_ = C_::value };

    typedef typename aux::if_impl< BOOST_MPL_AUX_STATIC_CAST(bool, msvc70_wknd_) >
        ::template result_<T1,T2>::type type;

    BOOST_MPL_AUX_LAMBDA_SUPPORT(3,if_,(C_,T1,T2))
};

#endif // BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(3, if_)


#if !defined(BOOST_MPL_NO_FULL_LAMBDA_SUPPORT)

// Aleksey, check it out: lazy if_ evaluation in lambdas!
// I think this doesn't handle the case of
//
//    _1<foo<_2>, bar<_2>, baz<_2> >
//
// (or however it is that you express that... when the ordinary bind3
// computes the function based on the actual arguments). That leads me
// to think that some kind of true currying might be a better
// approach, e.g.:
//
//
//  boost::mpl::bind3<
//         boost::mpl::quote3<boost::mpl::if_>
//       , boost::mpl::bind1<boost::mpl::quote1<boost::is_reference>, boost::mpl::arg<1> >
//       , boost::mpl::arg<1>
//       , boost::mpl::bind1<boost::mpl::quote1<add_ptr>, boost::mpl::arg<1> > 
//  >::apply<...>
//
// becomes:
//
//  boost::mpl::bind<
//         boost::mpl::quote3<boost::mpl::if_>
//  >::bind<
//       , boost::mpl::bind1<boost::mpl::quote1<boost::is_reference>,
//         boost::mpl::arg<1> >
//  >::bind<
//         boost::mpl::arg<1>
//  >::bind<
//       boost::mpl::bind1<boost::mpl::quote1<add_ptr>, boost::mpl::arg<1> > 
//  >::apply<...>
//
// so that after the 2nd bind we have a different function depending
// on the result of is_reference.

template <class T1, class T2, class T3, class T4> struct bind3;
template <template <class T1, class T2, class T3> class F, class tag> struct quote3;

namespace aux
{
  template <
      typename T
    , BOOST_MPL_PP_PARAMS(BOOST_MPL_METAFUNCTION_MAX_ARITY, typename U)
  > struct resolve_bind_arg;

  template<
        typename T
      , typename Arg
      >
  struct replace_unnamed_arg;
}

template<
      typename T1, typename T2, typename T3
    >
struct bind3<quote3<if_, void_>, T1, T2, T3>
{
    template<
          typename U1 = void_, typename U2 = void_, typename U3 = void_
        , typename U4 = void_, typename U5 = void_
        >
    struct apply
    {
     private:
        typedef quote3<if_, void_> a0;
        typedef mpl::arg< 1> n1;
        
        typedef aux::replace_unnamed_arg< T1,n1 > r1;
        typedef typename r1::type a1;
        typedef typename r1::next_arg n2;
        typedef typename aux::resolve_bind_arg< a1,U1,U2,U3,U4,U5 >::type t1;

        typedef aux::replace_unnamed_arg< T2,n2 > r2;
        typedef typename r2::type a2;
        typedef typename r2::next_arg n3;
        typedef typename aux::resolve_bind_arg< a2,U1,U2,U3,U4,U5 > f2;

        typedef aux::replace_unnamed_arg< T3,n3 > r3;
        typedef typename r3::type a3;
        typedef typename r3::next_arg n4;
        typedef typename aux::resolve_bind_arg< a3,U1,U2,U3,U4,U5 > f3;

        typedef typename if_<t1,f2,f3>::type f_;
     public:
        typedef typename f_::type type;
    };
};
#endif

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_IF_HPP_INCLUDED
