//-----------------------------------------------------------------------------
// boost mpl/zip_view.hpp header file
// See http://www.boost.org for updates, documentation, and revision history.
//-----------------------------------------------------------------------------
//
// Copyright (c) 2000-02
// David Abrahams, Aleksey Gurtovoy
//
// Permission to use, copy, modify, distribute and sell this software
// and its documentation for any purpose is hereby granted without fee, 
// provided that the above copyright notice appears in all copies and 
// that both the copyright notice and this permission notice appear in 
// supporting documentation. No representations are made about the 
// suitability of this software for any purpose. It is provided "as is" 
// without express or implied warranty.

#ifndef BOOST_MPL_ZIP_VIEW_HPP_INCLUDED
#define BOOST_MPL_ZIP_VIEW_HPP_INCLUDED

#include "boost/mpl/transform.hpp"
#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/iterator_tag.hpp"
#include "boost/mpl/next.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/aux_/void_spec.hpp"

namespace boost { namespace mpl {

template< typename IteratorSeq >
struct zip_iterator
{
    typedef input_iter_tag_ category;
    typedef typename transform<
          IteratorSeq
        , apply0<_1>
        >::type type;

    typedef zip_iterator<
          typename transform<
                IteratorSeq
              , next<_1>
            >::type
        > next;
};

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequences)
    >
struct zip_view
{
 private:    
    typedef typename transform< Sequences, begin<_1> >::type first_ones_;
    typedef typename transform< Sequences, end<_1> >::type last_ones_;
    
 public:
    typedef nested_begin_end_tag tag;
    typedef zip_iterator<first_ones_> begin;
    typedef zip_iterator<last_ones_> end;
};

BOOST_MPL_AUX_VOID_SPEC(1, zip_view)

}} // namespace boost::mpl

#endif // BOOST_MPL_ZIP_VIEW_HPP_INCLUDED
