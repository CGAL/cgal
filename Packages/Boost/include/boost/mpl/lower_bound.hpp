//-----------------------------------------------------------------------------
// boost mpl/lower_bound.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2001-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_LOWER_BOUND_HPP_INCLUDED
#define BOOST_MPL_LOWER_BOUND_HPP_INCLUDED

#include "boost/mpl/less.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

#if defined(__BORLANDC__) && (__BORLANDC__ <= 0x561 || !defined(BOOST_STRICT_CONFIG))
#   define BOOST_MPL_CFG_STRIPPED_DOWN_LOWER_BOUND_IMPL
#endif

#if !defined(BOOST_MPL_CFG_STRIPPED_DOWN_LOWER_BOUND_IMPL)
#   include "boost/mpl/minus.hpp"
#   include "boost/mpl/divides.hpp"
#   include "boost/mpl/size.hpp"
#   include "boost/mpl/advance.hpp"
#   include "boost/mpl/begin_end.hpp"
#   include "boost/mpl/integral_c.hpp"
#   include "boost/mpl/int.hpp"
#   include "boost/mpl/apply_if.hpp"
#   include "boost/mpl/apply.hpp"
#   include "boost/mpl/aux_/apply.hpp"
#   include "boost/mpl/aux_/deref_wknd.hpp"
#   include "boost/mpl/aux_/value_wknd.hpp"
#else
#   include "boost/mpl/not.hpp"
#   include "boost/mpl/find.hpp"
#   include "boost/mpl/bind.hpp"
#endif

#include "boost/config.hpp"

namespace boost {
namespace mpl {

#if defined(BOOST_MPL_CFG_STRIPPED_DOWN_LOWER_BOUND_IMPL)

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

// agurt 23/oct/02: has a wrong complexity etc., but at least it works
// feel free to contribute a better implementation!
template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    , typename Predicate = less<>
    , typename pred_ = typename lambda<Predicate>::type
    >
struct lower_bound
    : find_if< Sequence, bind1< not_<>, bind2<pred_,_,T> > >
{
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

#else

namespace aux {

template<
      typename Distance
    , typename Predicate
    , typename T
    , typename DeferredIterator
    >
struct lower_bound_step_impl;

template< 
      typename Distance
    , typename Predicate
    , typename T
    , typename DeferredIterator
    >
struct lower_bound_step
{
    typedef typename apply_if<
          Distance
        , lower_bound_step_impl<Distance,Predicate,T,DeferredIterator>
        , apply0<DeferredIterator>
        >::type type;
};
    
template<
      typename Distance
    , typename Predicate
    , typename T
    , typename DeferredIterator
    >
struct lower_bound_step_impl
{
    typedef typename divides< Distance, integral_c<long,2> >::type offset_;
    typedef typename DeferredIterator::type iter_;
    typedef typename advance< iter_,offset_ >::type middle_;
    typedef typename BOOST_MPL_AUX_APPLY2(
              Predicate
            , typename BOOST_MPL_AUX_DEREF_WNKD(middle_)
            , T
            )::type cond_;

    typedef typename minus< Distance, offset_, integral_c<long,1> >::type step_;
    typedef lower_bound_step< offset_,Predicate,T,DeferredIterator > step_forward_;
    typedef lower_bound_step< step_,Predicate,T,next<middle_> > step_backward_;
    typedef typename apply_if<
          cond_
        , step_backward_
        , step_forward_
        >::type type;
};


} // namespace aux

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(T)
    , typename Predicate = less<>
    >
struct lower_bound
{
 private:
    typedef typename lambda<Predicate>::type pred_;
    typedef typename size<Sequence>::type size_;

 public:
    typedef typename aux::lower_bound_step<
        size_,pred_,T,begin<Sequence>
        >::type type;
};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

#endif // BOOST_MPL_CFG_STRIPPED_DOWN_LOWER_BOUND_IMPL

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(2, lower_bound)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_LOWER_BOUND_HPP_INCLUDED
