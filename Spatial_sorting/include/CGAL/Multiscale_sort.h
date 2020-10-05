// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Christophe Delage

#ifndef CGAL_MULTISCALE_SORT_H
#define CGAL_MULTISCALE_SORT_H

#include <CGAL/config.h>
#include <CGAL/assertions.h>
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
      middle = begin + difference_type (double(end - begin) * _ratio);
      this->operator() (begin, middle);
    }
    _sort (middle, end);
  }
};

} // namespace CGAL

#endif // CGAL_MULTISCALE_SORT_H
