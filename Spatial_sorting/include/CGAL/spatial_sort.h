// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Christophe Delage

#ifndef CGAL_SPATIAL_SORT_H
#define CGAL_SPATIAL_SORT_H

#include <CGAL/basic.h>

#include <CGAL/Hilbert_sort_2.h>
#include <CGAL/Hilbert_sort_3.h>

#include <CGAL/Multiscale_sort.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

    template <class RandomAccessIterator, class Kernel>
    void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                       const Kernel &k, typename Kernel::Point_2 *)
    {
        typedef Hilbert_sort_2<Kernel> Sort;
        (Multiscale_sort<Sort> (Sort (k, 4), 16, 0.25)) (begin, end);
    }

    template <class RandomAccessIterator, class Kernel>
    void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                       const Kernel &k, typename Kernel::Point_3 *)
    {
        typedef Hilbert_sort_3<Kernel> Sort;
        (Multiscale_sort<Sort> (Sort (k, 8), 64, 0.125)) (begin, end);
    }
}

template <class RandomAccessIterator, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;

    CGALi::spatial_sort (begin, end, k, static_cast<value_type *> (0));
}

template <class RandomAccessIterator>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    spatial_sort (begin, end, Kernel());
}

CGAL_END_NAMESPACE

#endif//CGAL_SPATIAL_SORT_H
