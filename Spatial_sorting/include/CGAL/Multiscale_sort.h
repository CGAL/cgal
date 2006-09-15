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
