/*=============================================================================
    Copyright (c) 1999-2003 Jaakko Järvi
    Copyright (c) 1999-2003 Jeremiah Willcock
    Copyright (c) 2001-2003 Joel de Guzman

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#if !defined(BOOST_IO_HPP)
#define BOOST_IO_HPP

#include <iostream>
#include <boost/spirit/fusion/sequence/detail/io.hpp>
#include <boost/spirit/fusion/sequence/as_fusion_sequence.hpp>

namespace boost { namespace fusion
{
    ///////////////////////////////////////////////////////////////////////////
    //
    //  Sequence I/O (<< and >> operators)
    //
    ///////////////////////////////////////////////////////////////////////////
    template <typename OStream, typename Sequence>
    inline OStream&
    operator<<(OStream& os, Sequence const& seq)
    {
        detail::print_sequence(os,
            as_fusion_sequence<Sequence const>::convert(seq));
        return os;
    }

    template <typename IStream, typename Sequence>
    inline IStream&
    operator>>(IStream& is, Sequence& seq)
    {
        detail::read_sequence(is,
            as_fusion_sequence<Sequence>::convert(seq));
        return is;
    }
}}

#endif
