#ifndef CGAL_BRIO_SORT_H
#define CGAL_BRIO_SORT_H

#include <CGAL/basic.h>

#include <CGAL/Hilbert_sort_2.h>
#include <CGAL/Hilbert_sort_3.h>

CGAL_BEGIN_NAMESPACE

template <class Sort, class RandomAccessIterator>
class BRIO_sort
{
    Sort sort;
    int threshold;
    int ratio;
public:
    BRIO_sort (int _ratio = 8, int _threshold = 64, const Sort &_sort = Sort())
        : sort (_sort), threshold (_threshold), ratio (_ratio) 
    {}
    
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

template <class RandomAccessIterator, class Kernel>
void brio_sort_2 (RandomAccessIterator begin, RandomAccessIterator end, Kernel k)
{
    typedef Hilbert_sort_2<Kernel, RandomAccessIterator> Sort;
    (BRIO_sort<Sort,RandomAccessIterator> (4, 16, Sort (k))) (begin, end);
}

template <class RandomAccessIterator>
void brio_sort_2 (RandomAccessIterator begin, RandomAccessIterator end)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    brio_sort_2 (begin, end, Kernel ());
}

template <class RandomAccessIterator, class Kernel>
void brio_sort_3 (RandomAccessIterator begin, RandomAccessIterator end, Kernel k)
{
    typedef Hilbert_sort_3<Kernel, RandomAccessIterator> Sort;
    (BRIO_sort<Sort,RandomAccessIterator> (8, 64, Sort (k))) (begin, end);
}

template <class RandomAccessIterator>
void brio_sort_3 (RandomAccessIterator begin, RandomAccessIterator end)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    brio_sort_3 (begin, end, Kernel ());
}

CGAL_END_NAMESPACE


#endif//CGAL_BRIO_SORT_H
