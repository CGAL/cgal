//-----------------------------------------------------------------------------
// boost mpl/sequence_tag.hpp header file
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

#ifndef BOOST_MPL_SEQUENCE_TAG_HPP_INCLUDED
#define BOOST_MPL_SEQUENCE_TAG_HPP_INCLUDED

#include "boost/mpl/sequence_tag_fwd.hpp"
#include "boost/mpl/aux_/has_tag.hpp"
#include "boost/mpl/aux_/has_begin.hpp"
#include "boost/mpl/aux_/void_spec.hpp"
#include "boost/mpl/aux_/is_msvc_eti_arg.hpp"
#include "boost/mpl/aux_/config/eti.hpp"
#include "boost/mpl/aux_/yes_no.hpp"
#include "boost/mpl/aux_/config/workaround.hpp"

namespace boost { namespace mpl {

// agurt, 27/nov/02: have to use a simplistic 'sequence_tag' implementation
// on MSVC to avoid dreadful "internal structure overflow" error
#if BOOST_WORKAROUND(BOOST_MSVC, < 1300)

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    >
struct sequence_tag
{
    typedef typename Sequence::tag type;
};

#elif BOOST_WORKAROUND(BOOST_MSVC, == 1300)

// agurt, 07/feb/03: workaround for what seems to be MSVC 7.0-specific ETI issue

namespace aux {

template< bool >
struct sequence_tag_impl
{
    template< typename Sequence > struct result_
    {
        typedef typename Sequence::tag type;
    };
};

template<>
struct sequence_tag_impl<false>
{
    template< typename Sequence > struct result_
    {
        typedef int type;
    };
};

} // namespace aux

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    >
struct sequence_tag
    : aux::sequence_tag_impl< !aux::is_msvc_eti_arg<Sequence>::value >
        ::template result_<Sequence>
{
};

#else

namespace aux {

template< bool has_tag_, bool has_begin_ >
struct sequence_tag_impl
{
    // agurt 24/nov/02: MSVC 6.5 gets confused in 'sequence_tag_impl<true>' 
    // specialization below, if we name it 'result_' here
    template< typename Sequence > struct result2_;
};

#   define AUX_CLASS_SEQUENCE_TAG_SPEC(has_tag, has_begin, result_type) \
template<> struct sequence_tag_impl<has_tag,has_begin> \
{ \
    template< typename Sequence > struct result2_ \
    { \
        typedef result_type type; \
    }; \
}; \
/**/

AUX_CLASS_SEQUENCE_TAG_SPEC(true, true, typename Sequence::tag)
AUX_CLASS_SEQUENCE_TAG_SPEC(true, false, typename Sequence::tag)
AUX_CLASS_SEQUENCE_TAG_SPEC(false, true, nested_begin_end_tag)
AUX_CLASS_SEQUENCE_TAG_SPEC(false, false, non_sequence_tag)

#   undef AUX_CLASS_SEQUENCE_TAG_SPEC

} // namespace aux

template<
      typename BOOST_MPL_AUX_VOID_SPEC_PARAM(Sequence)
    >
struct sequence_tag
    : aux::sequence_tag_impl<
          ::boost::mpl::aux::has_tag<Sequence>::value
        , ::boost::mpl::aux::has_begin<Sequence>::value
        >::template result2_<Sequence>
{
};

#endif // BOOST_MSVC

#if defined(BOOST_MPL_MSVC_60_ETI_BUG)
template<> struct sequence_tag<int>
{
    typedef int type;
};
#endif

BOOST_MPL_AUX_VOID_SPEC(1, sequence_tag)

}} // namespace boost::mpl

#endif // BOOST_MPL_SEQUENCE_TAG_HPP_INCLUDED
