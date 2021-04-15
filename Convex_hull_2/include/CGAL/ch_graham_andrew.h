// Copyright (c) 1999  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra


#ifndef CGAL_CH_GRAHAM_ANDREW_H
#define CGAL_CH_GRAHAM_ANDREW_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <iterator>


namespace CGAL {

// computes the sorted sequence of extreme points which are not left
// of $pq$ and reports this sequence in a range starting at |result|,
// where $p$ is the value of |first| and $q$ is the value of |last| $-1$.
// The sequence reported starts with $p$, point $q$ is omitted.
// {\it Precondition:} The points in [|first|,|last|) are sorted with respect
// to $pq$ and the range [|first|,|last|) contains at least two different
// points.
// {\sc traits}: uses |Traits::Left_turn_2| and |Traits::Equal_2| operating on the
// point type |Traits::Point_2|.
template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
                       BidirectionalIterator last,
                       OutputIterator        result,
                       const Traits& ch_traits );

template <class BidirectionalIterator, class OutputIterator>
inline
OutputIterator
ch_graham_andrew_scan( BidirectionalIterator first,
                       BidirectionalIterator last,
                       OutputIterator        result )
{
    typedef std::iterator_traits<BidirectionalIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_graham_andrew_scan( first, last, result, Kernel());
}

template <class BidirectionalIterator, class OutputIterator, class Traits>
OutputIterator
ch__ref_graham_andrew_scan( BidirectionalIterator first,
                            BidirectionalIterator last,
                            OutputIterator&       result,
                            const Traits&         ch_traits);


// same as |convex_hull_2(first,last,result)|.
// {\sc traits}: uses |Traits::Point_2|, |Traits::Left_turn_2|
// and |Traits::Less_xy_2|.
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_graham_andrew( InputIterator  first,
                  InputIterator  last,
                  OutputIterator result,
                  const Traits&  ch_traits );

template <class InputIterator, class OutputIterator>
inline
OutputIterator
ch_graham_andrew( InputIterator  first,
                  InputIterator  last,
                  OutputIterator result )
{
    typedef std::iterator_traits<InputIterator>   ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return ch_graham_andrew( first, last, result, Kernel());
}



template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_lower_hull_scan( InputIterator  first,
                    InputIterator  last,
                    OutputIterator result,
                    const Traits&  ch_traits);

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
ch_upper_hull_scan( InputIterator  first,
                    InputIterator  last,
                    OutputIterator result,
                    const Traits&  ch_traits);

} //namespace CGAL

#include <CGAL/Convex_hull_2/ch_graham_andrew_impl.h>

#endif // CGAL_CH_GRAHAM_ANDREW_H
