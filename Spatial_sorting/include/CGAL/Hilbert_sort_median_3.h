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
// $URL: svn+ssh://odevil@scm.gforge.inria.fr/svn/cgal/trunk/Spatial_sorting/include/CGAL/Hilbert_sort_3.h $
// $Id: Hilbert_sort_3.h 51456 2009-08-24 17:10:04Z spion $
//
// Author(s)     : Christophe Delage

#ifndef CGAL_HILBERT_SORT_MEDIAN_3_H
#define CGAL_HILBERT_SORT_MEDIAN_3_H

#include <CGAL/config.h>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_base.h>

namespace CGAL {

namespace internal {
    template <class K, int x, bool up> struct Hilbert_cmp_3;

    template <class K, int x>
    struct Hilbert_cmp_3<K,x,true>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        {
            return Hilbert_cmp_3<K,x,false> (k) (q, p);
        }
    };

    template <class K>
    struct Hilbert_cmp_3<K,0,false>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        {
            return k.less_x_3_object() (p, q);
        }
    };

    template <class K>
    struct Hilbert_cmp_3<K,1,false>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        {
            return k.less_y_3_object() (p, q);
        }
    };

    template <class K>
    struct Hilbert_cmp_3<K,2,false>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        Hilbert_cmp_3 (const K &_k = K()) : k(_k) {}
        bool operator() (const Point &p, const Point &q) const
        {
            return k.less_z_3_object() (p, q);
        }
    };
}

template <class K>
class Hilbert_sort_median_3
{
public:
    typedef K Kernel;
    typedef typename Kernel::Point_3 Point;

private:
    Kernel _k;
    std::ptrdiff_t _limit;

    template <int x, bool up> struct Cmp : public internal::Hilbert_cmp_3<Kernel,x,up>
    { Cmp (const Kernel &k) : internal::Hilbert_cmp_3<Kernel,x,up> (k) {} };

public:
    Hilbert_sort_median_3 (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
        : _k(k), _limit (limit)
    {}

    template <int x, bool upx, bool upy, bool upz, class RandomAccessIterator>
    void sort (RandomAccessIterator begin, RandomAccessIterator end) const
    {
        const int y = (x + 1) % 3, z = (x + 2) % 3;
        if (end - begin <= _limit) return;

        RandomAccessIterator m0 = begin, m8 = end;

        RandomAccessIterator m4 = internal::hilbert_split (m0, m8, Cmp< x,  upx> (_k));
        RandomAccessIterator m2 = internal::hilbert_split (m0, m4, Cmp< y,  upy> (_k));
        RandomAccessIterator m1 = internal::hilbert_split (m0, m2, Cmp< z,  upz> (_k));
        RandomAccessIterator m3 = internal::hilbert_split (m2, m4, Cmp< z, !upz> (_k));
        RandomAccessIterator m6 = internal::hilbert_split (m4, m8, Cmp< y, !upy> (_k));
        RandomAccessIterator m5 = internal::hilbert_split (m4, m6, Cmp< z,  upz> (_k));
        RandomAccessIterator m7 = internal::hilbert_split (m6, m8, Cmp< z, !upz> (_k));

        sort<z, upz, upx, upy> (m0, m1);
        sort<y, upy, upz, upx> (m1, m2);
        sort<y, upy, upz, upx> (m2, m3);
        sort<x, upx,!upy,!upz> (m3, m4);
        sort<x, upx,!upy,!upz> (m4, m5);
        sort<y,!upy, upz,!upx> (m5, m6);
        sort<y,!upy, upz,!upx> (m6, m7);
        sort<z,!upz,!upx, upy> (m7, m8);
    }

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
        sort <0, false, false, false> (begin, end);
    }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MEDIAN_3_H
