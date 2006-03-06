#ifndef CGAL_MULTISCALE_SORT_H
#define CGAL_MULTISCALE_SORT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template <class Sort>
class Multiscale_sort
{
    Sort sort;
    int threshold;
    int ratio;

public:
    Multiscale_sort (int _ratio = 8, int _threshold = 64, const Sort &_sort = Sort())
        : sort (_sort), threshold (_threshold), ratio (_ratio)
    {}

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
        RandomAccessIterator middle = begin;
        if (end - begin >= threshold) {
            middle = begin + (end - begin) / ratio;
            this->operator() (begin, middle);
        }
        sort (middle, end);
    }
};

CGAL_END_NAMESPACE

#endif//CGAL_MULTISCALE_SORT_H
