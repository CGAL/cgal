//-----------------------------------------------------------------------------
// boost mpl/for_each.hpp header file
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

#ifndef BOOST_MPL_FOR_EACH_HPP_INCLUDED
#define BOOST_MPL_FOR_EACH_HPP_INCLUDED

#include "boost/mpl/begin_end.hpp"
#include "boost/mpl/apply.hpp"
#include "boost/mpl/bool.hpp"
#include "boost/mpl/lambda.hpp"
#include "boost/mpl/identity.hpp"
#include "boost/mpl/aux_/unwrap.hpp"

#include "boost/type_traits/is_same.hpp"
#include "boost/utility/value_init.hpp"

namespace boost {
namespace mpl {

namespace aux {

template <bool done = true>
struct for_each_impl
{
    template<
          typename Iterator
        , typename LastIterator
        , typename TransformFunc
        , typename F
        >
    static void execute(
          Iterator*
        , LastIterator*
        , TransformFunc*
        , F
        )
    {
    }
};

template <>
struct for_each_impl<false>
{
    template<
          typename Iterator
        , typename LastIterator
        , typename TransformFunc
        , typename F
        >
    static void execute(
          Iterator*
        , LastIterator*
        , TransformFunc* 
        , F f
        )
    {
        typedef typename Iterator::type item;
        typedef typename apply1<TransformFunc,item>::type arg;
    
        // dwa 2002/9/10 -- make sure not to invoke undefined behavior
        // when we pass arg.
        value_initialized<arg> x;
        aux::unwrap(f, 0)(boost::get(x));
        
        typedef typename Iterator::next iter;
        for_each_impl<boost::is_same<iter,LastIterator>::value>::execute(
            (iter*)0, (LastIterator*)0, (TransformFunc*)0, f);
    }
};

} // namespace aux

// agurt, 17/mar/02: pointer default parameters are necessary to workaround 
// MSVC 6.5 function template signature's mangling bug
template<
      typename Sequence
    , typename TransformOp
    , typename F
    >
inline
void for_each(F f, Sequence* = 0, TransformOp* = 0)
{
    typedef typename begin<Sequence>::type first;
    typedef typename end<Sequence>::type last;
    typedef typename lambda<TransformOp>::type transform_op;

    aux::for_each_impl< boost::is_same<first,last>::value >::execute(
        (first*)0, (last*)0, (transform_op*)0, f);
}

template<
      typename Sequence
    , typename F
    >
inline
void for_each(F f, Sequence* = 0)
{
    for_each<Sequence, identity<> >(f);
}

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_FOR_EACH_HPP_INCLUDED
