// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : random_selection.h
// chapter       : $CGAL_Chapter: Geometric Object Generators $
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : generators.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Copy n randomly chosen items
// ============================================================================

#ifndef CGAL_RANDOM_SELECTION_H
#define CGAL_RANDOM_SELECTION_H 1
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_PROTECT_ITERATOR
#include <iterator>
#define CGAL_PROTECT_ITERATOR
#endif
#ifndef CGAL_RANDOM_H
#include <CGAL/Random.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class RandomAccessIterator, class Size, class OutputIterator,
          class Random>
OutputIterator random_selection( RandomAccessIterator first,
                                 RandomAccessIterator last,
                                 Size n,
                                 OutputIterator result,
                                 Random& rnd)
    // choose a random item from the range [`first',`last') and write it
    // to `first2', each item from the range with equal probability.
    // Repeat this n times, thus writing n items to `first2'. A single
    // random number is needed from `rnd' for each item. Returns the
    // value of `first2' after inserting the n items.
{
    int m = int(last - first);
    for ( Size i = 0; i < n; i++) {
        *result++ = first[ rnd(m)];
    }
    return result;
}

template <class RandomAccessIterator, class Size, class OutputIterator>
OutputIterator random_selection( RandomAccessIterator first,
                                 RandomAccessIterator last,
                                 Size n,
                                 OutputIterator result)
{
    return random_selection( first, last, n, result, default_random);
}

CGAL_END_NAMESPACE    
#endif // CGAL_RANDOM_SELECTION_H //
// EOF //
