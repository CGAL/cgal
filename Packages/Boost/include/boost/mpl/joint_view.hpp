
#ifndef BOOST_MPL_JOINT_VIEW_HPP_INCLUDED
#define BOOST_MPL_JOINT_VIEW_HPP_INCLUDED

// + file: boost/mpl/joint_view.hpp
// + last modified: 25/may/03

// Copyright (c) 2000-03
// Aleksey Gurtovoy
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

#include "boost/mpl/aux_/joint_iter.hpp"
#include "boost/mpl/plus.hpp"
#include "boost/mpl/size_fwd.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost {
namespace mpl {

namespace aux {
struct joint_view_tag;
}

template<>
struct size_traits< aux::joint_view_tag >
{
    template < typename JointView > struct algorithm
      : plus<
            size<typename JointView::sequence1_>
          , size<typename JointView::sequence2_>
          >
    {};
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence1_)
    , typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence2_)
    >
struct joint_view
{
    typedef typename begin<Sequence1_>::type first1_;
    typedef typename end<Sequence1_>::type last1_;
    typedef typename begin<Sequence2_>::type first2_;
    typedef typename end<Sequence2_>::type last2_;

 public:
    // agurt, 25/may/03: for the 'size_traits' implementation above
    typedef Sequence1_ sequence1_;
    typedef Sequence2_ sequence2_;

    typedef aux::joint_view_tag tag;
    typedef typename aux::joint_iter<first1_,last1_,first2_> begin;
    typedef typename aux::joint_iter<last1_,last1_,last2_> end;
};

BOOST_MPL_AUX_VOID_SPEC(2, joint_view)

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_JOINT_VIEW_HPP_INCLUDED
