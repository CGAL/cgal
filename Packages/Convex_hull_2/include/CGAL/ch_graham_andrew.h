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
// file          : include/CGAL/ch_graham_andrew.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_CH_GRAHAM_ANDREW_H
#define CGAL_CH_GRAHAM_ANDREW_H

#include <CGAL/basic.h>
#include <iterator>


CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/ch_graham_andrew.C>
#endif // CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION

#endif // CGAL_CH_GRAHAM_ANDREW_H

