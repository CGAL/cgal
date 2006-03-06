#ifndef CGAL_HILBERT_SORT_BASE_H
#define CGAL_HILBERT_SORT_BASE_H

#include <CGAL/basic.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

    template <class RandomAccessIterator, class Cmp>
    RandomAccessIterator
    hilbert_split (RandomAccessIterator begin, RandomAccessIterator end,
                   Cmp cmp = Cmp ())
    {
        if (begin >= end) return begin;

        RandomAccessIterator middle = begin + (end - begin) / 2;
        std::nth_element (begin, middle, end, cmp);
        return middle;
    }
}

CGAL_END_NAMESPACE

#endif//CGAL_HILBERT_SORT_BASE_H
