#ifndef CGAL_HILBERT_SORT_2_H
#define CGAL_HILBERT_SORT_2_H

#include <CGAL/basic.h>

#include <CGAL/Hilbert_sort_base.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
    template <class K, int x, bool up> struct Hilbert_cmp_2;


    template <class K, int x>
    struct Hilbert_cmp_2<K,x,true>
        : public std::binary_function<typename K::Point_2,
                                      typename K::Point_2, bool>
    {
        typedef typename K::Point_2 Point;
        K k;
        Hilbert_cmp_2 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        { 
            return Hilbert_cmp_2<K,x,false> (k) (q, p);
        }
    };
    
    template <class K>
    struct Hilbert_cmp_2<K,0,false>
        : public std::binary_function<typename K::Point_2,
                                      typename K::Point_2, bool>
    {
        typedef typename K::Point_2 Point;
        K k;
        Hilbert_cmp_2 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        { 
            return k.less_x_2_object() (p, q);
        }
    };
    
    template <class K>
    struct Hilbert_cmp_2<K,1,false>
        : public std::binary_function<typename K::Point_2,
                                      typename K::Point_2, bool>
    {
        typedef typename K::Point_2 Point;
        K k;
        Hilbert_cmp_2 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        { 
            return k.less_y_2_object() (p, q);
        }
    };
};

template <class K, class RandomAccessIterator>
class Hilbert_sort_2
{
public:
    typedef RandomAccessIterator Iterator;
    typedef K Kernel;
    typedef typename Kernel::Point_2 Point;
    
private:
    Kernel k;

    template <int x, bool up> struct Cmp : public CGALi::Hilbert_cmp_2<Kernel,x,up>
    { Cmp (const Kernel &k) : CGALi::Hilbert_cmp_2<Kernel,x,up> (k) {} };

public:
    Hilbert_sort_2 (const Kernel &_k = Kernel()) : k(_k) {}

    template <int x, bool upx, bool upy>
    void sort (Iterator begin, Iterator end) const
    {
        const int y = (x + 1) % 2;
        if (end - begin <= 8) return;

        Iterator m0 = begin, m4 = end;

        Iterator m2 = CGALi::hilbert_split (m0, m4, Cmp< x,  upx> (k));
        Iterator m1 = CGALi::hilbert_split (m0, m2, Cmp< y,  upy> (k));
        Iterator m3 = CGALi::hilbert_split (m2, m4, Cmp< y, !upy> (k));

        if (end - begin <= 8*4) return;

        sort<y, upy, upx> (m0, m1);
        sort<x, upx, upy> (m1, m2);
        sort<x, upx, upy> (m2, m3);
        sort<y,!upy,!upx> (m3, m4);
    }

    void operator() (Iterator begin, Iterator end) const
    {
        sort <0, false, false> (begin, end);
    }
};

template <class RandomAccessIterator, class Kernel>
void hilbert_sort_2 (RandomAccessIterator begin, RandomAccessIterator end,
                     Kernel k)
{
    (Hilbert_sort_2<Kernel,RandomAccessIterator> (k)) (begin, end);
}

template <class RandomAccessIterator>
void hilbert_sort_2 (RandomAccessIterator begin, RandomAccessIterator end)
{
    typedef std::iterator_traits<RandomAccessIterator> ITraits;
    typedef typename ITraits::value_type               value_type;
    typedef CGAL::Kernel_traits<value_type>            KTraits;
    typedef typename KTraits::Kernel                   Kernel;

    hilbert_sort_2 (begin, end, Kernel());
}

CGAL_END_NAMESPACE

#endif//CGAL_HILBERT_SORT_2_H
