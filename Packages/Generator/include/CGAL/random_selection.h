// Copyright (c) 1997  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

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
