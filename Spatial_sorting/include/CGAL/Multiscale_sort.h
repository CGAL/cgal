// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_MULTISCALE_SORT_H
#define CGAL_MULTISCALE_SORT_H

#include <CGAL/basic.h>
#include <iterator>
#include <cstddef>

namespace CGAL {

template <class Sort>
class Multiscale_sort
{
    Sort _sort;
    std::ptrdiff_t _threshold;
    double _ratio;

public:
    Multiscale_sort (const Sort &sort = Sort(), std::ptrdiff_t threshold = 1, double ratio = 0.5)
        : _sort (sort), _threshold (threshold), _ratio (ratio)
    {
        CGAL_precondition (0. <= ratio && ratio <= 1.);
    }

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
	typedef typename std::iterator_traits<RandomAccessIterator>::difference_type difference_type;
        RandomAccessIterator middle = begin;
        if (end - begin >= _threshold) {
            middle = begin + difference_type ((end - begin) * _ratio);
            this->operator() (begin, middle);
        }
        _sort (middle, end);
    }
};

} // namespace CGAL

#endif // CGAL_MULTISCALE_SORT_H
