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

#ifndef CGAL_MULTISCALE_SORT_H
#define CGAL_MULTISCALE_SORT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template <class Sort>
class Multiscale_sort
{
    Sort _sort;
    int _threshold;
    double _ratio;

public:
    Multiscale_sort (const Sort &sort = Sort(), int threshold = 1, double ratio = 0.5)
        : _sort (sort), _threshold (threshold), _ratio (ratio)
    {
        CGAL_precondition (0. <= ratio && ratio <= 1.);
    }

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
        RandomAccessIterator middle = begin;
        if (end - begin >= _threshold) {
            middle = begin + int ((end - begin) * _ratio);
            this->operator() (begin, middle);
        }
        _sort (middle, end);
    }
};

CGAL_END_NAMESPACE

#endif//CGAL_MULTISCALE_SORT_H
