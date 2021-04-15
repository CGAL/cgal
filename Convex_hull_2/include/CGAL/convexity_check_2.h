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


#ifndef CGAL_CONVEXITY_CHECK_2_H
#define CGAL_CONVEXITY_CHECK_2_H

#include <CGAL/license/Convex_hull_2.h>


#include <CGAL/basic.h>
#include <iterator>

namespace CGAL {

// returns true, if the point elements in [|first|,|last|) form a
// counterclockwise oriented strongly convex polygon. Strongly means,
// there are no three collinear points.
// {\sc traits}: uses |Traits::Left_turn_2|, |Traits::Equal_2| and |Traits::Less_xy_2|.
template <class ForwardIterator, class Traits>
bool
is_ccw_strongly_convex_2( ForwardIterator first, ForwardIterator last,
                          const Traits& ch_traits);

template <class ForwardIterator>
inline
bool
is_ccw_strongly_convex_2( ForwardIterator first, ForwardIterator last )
{
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return is_ccw_strongly_convex_2( first, last, Kernel());
}




// returns true, if the point elements in [|first|,|last|) form a
// clockwise oriented strongly convex polygon. Strongly means, there are
// no three collinear points.
// {\sc traits}: uses |Traits::Left_turn_2|, |Traits::Equal_2| and |Traits::Less_xy_2|.
template <class ForwardIterator, class Traits>
bool
is_cw_strongly_convex_2( ForwardIterator first, ForwardIterator last,
                         const Traits& ch_traits);

template <class ForwardIterator>
inline
bool
is_cw_strongly_convex_2( ForwardIterator first, ForwardIterator last )
{
    typedef std::iterator_traits<ForwardIterator> ITraits;
    typedef typename ITraits::value_type          value_type;
    typedef CGAL::Kernel_traits<value_type>       KTraits;
    typedef typename KTraits::Kernel              Kernel;
    return is_cw_strongly_convex_2( first, last, Kernel());
}


// returns true, if all points in [|first1|,|last1|) are
// not right of the lines defined by consecutive points in the range
// [|first2|,|last2|), where the range is considered as a cycle.
// {\sc traits}: uses |Traits::Left_turn_2|.
template <class ForwardIterator1, class ForwardIterator2, class Traits>
bool
ch_brute_force_check_2(ForwardIterator1 first1, ForwardIterator1 last1,
                       ForwardIterator2 first2, ForwardIterator2 last2,
                       const Traits& ch_traits);

template <class ForwardIterator1, class ForwardIterator2>
inline
bool
ch_brute_force_check_2(ForwardIterator1 first1, ForwardIterator1 last1,
                       ForwardIterator2 first2, ForwardIterator2 last2)
{
    typedef std::iterator_traits<ForwardIterator1> ITraits;
    typedef typename ITraits::value_type           value_type;
    typedef CGAL::Kernel_traits<value_type>        KTraits;
    typedef typename KTraits::Kernel               Kernel;
    return ch_brute_force_check_2( first1, last1, first2, last2, Kernel() );
}


// returns true, if all points in [|first1|,|last1|) are
// not right of the lines defined by consecutive points in the range
// [|first2|,|last2|).
// {\sc traits}: uses |Traits::Left_turn_2|.
template <class ForwardIterator1, class ForwardIterator2, class Traits>
bool
ch_brute_force_chain_check_2(ForwardIterator1 first1,
                             ForwardIterator1 last1,
                             ForwardIterator2 first2,
                             ForwardIterator2 last2,
                             const Traits& ch_traits);

template <class ForwardIterator1, class ForwardIterator2>
inline
bool
ch_brute_force_chain_check_2(ForwardIterator1 first1, ForwardIterator1 last1,
                             ForwardIterator2 first2, ForwardIterator2 last2)
{
    typedef std::iterator_traits<ForwardIterator1> ITraits;
    typedef typename ITraits::value_type           value_type;
    typedef CGAL::Kernel_traits<value_type>        KTraits;
    typedef typename KTraits::Kernel               Kernel;
    return ch_brute_force_chain_check_2( first1, last1, first2, last2,
                                         Kernel());
}

} //namespace CGAL

#include <CGAL/Convex_hull_2/convexity_check_2_impl.h>

#endif // CGAL_CONVEXITY_CHECK_2_H
