/*=============================================================================
    Copyright (c) 1999-2003 Jaakko Järvi
    Copyright (c) 2001-2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(FUSION_SEQUENCE_GET_HPP)
#define FUSION_SEQUENCE_GET_HPP

#include <boost/spirit/fusion/sequence/at.hpp>
#include <boost/spirit/fusion/sequence/detail/sequence_base.hpp>

namespace boost { namespace fusion
{
    ///////////////////////////////////////////////////////////////////////////
    //
    //  get function
    //
    //      Given a constant integer N and a sequence, returns a reference to
    //      the Nth element of the sequence. (N is a zero based index). Usage:
    //
    //          get<N>(seq)
    //
    //  This function is provided here for compatibility with the tuples TR1
    //  specification. This function forwards to at<N>(seq).
    //
    ///////////////////////////////////////////////////////////////////////////
    template <int N, typename Sequence>
    inline typename meta::at_c<Sequence const, N>::type
    get(sequence_base<Sequence> const& seq FUSION_GET_MSVC_WORKAROUND)
    {
        return at<N>(seq.cast());
    }

    template <int N, typename Sequence>
    inline typename meta::at_c<Sequence, N>::type
    get(sequence_base<Sequence>& seq FUSION_GET_MSVC_WORKAROUND)
    {
        return at<N>(seq.cast());
    }
}}

#endif

