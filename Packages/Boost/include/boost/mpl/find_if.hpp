//-----------------------------------------------------------------------------
// boost mpl/find_if.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_FIND_IF_HPP_INCLUDED
#define BOOST_MPL_FIND_IF_HPP_INCLUDED

#include "boost/mpl/aux_/find_if_pred.hpp"
#include "boost/mpl/arg.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/iter_fold_if.hpp"
#include "boost/mpl/protect.hpp"
#include "boost/mpl/aux_/common_name_wknd.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/lambda_support.hpp"

namespace boost {
namespace mpl {

BOOST_MPL_AUX_COMMON_NAME_WKND(find_if)

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_BEGIN

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Predicate)
    >
struct find_if
{
private:

    typedef typename lambda<Predicate>::type pred_;

    typedef typename iter_fold_if<
          Sequence
        , void
        , protect< arg<1> > // ignore
        , protect< aux::find_if_pred<pred_> >
        >::type result_;

public:

    typedef typename result_::second type;

    BOOST_MPL_AUX_LAMBDA_SUPPORT(2,find_if,(Sequence,Predicate))

};

BOOST_MPL_AUX_AGLORITHM_NAMESPACE_END

BOOST_MPL_AUX_ALGORITHM_VOID_SPEC(2,find_if)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_FIND_IF_HPP_INCLUDED
