// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/convexity_check_2.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CONVEXITY_CHECK_2_H
#define CGAL_CONVEXITY_CHECK_2_H

#include <CGAL/basic.h>
#include <iterator>

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/convexity_check_2.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CONVEXITY_CHECK_2_H

