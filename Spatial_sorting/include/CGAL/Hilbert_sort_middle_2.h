// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
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
// Author(s)     :  Olivier Devillers

#ifndef CGAL_HILBERT_SORT_MIDDLE_2_H
#define CGAL_HILBERT_SORT_MIDDLE_2_H

#include <CGAL/basic.h>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_middle_base.h>
#include <CGAL/Polygon_2_algorithms.h> 

namespace CGAL {

namespace internal {
    template <class K, int x, bool up> struct Fixed_hilbert_cmp_2;

    template <class K, int x>
    struct Fixed_hilbert_cmp_2<K,x,true>
        : public std::binary_function<typename K::Point_2,
                                      typename K::Point_2, bool>
    {
        typedef typename K::Point_2 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_2 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        { 
	  return ! Fixed_hilbert_cmp_2<K,x,false> (value, k) (p);
        }
    };
    
    template <class K>
    struct Fixed_hilbert_cmp_2<K,0,false>
        : public std::binary_function<typename K::Point_2,
                                      typename K::Point_2, bool>
    {
        typedef typename K::Point_2 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_2 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        { 
	  return to_double(k.compute_x_2_object()(p)) < value;
        }
    };
    
    template <class K>
    struct Fixed_hilbert_cmp_2<K,1,false>
        : public std::binary_function<typename K::Point_2,
                                      typename K::Point_2, bool>
    {
        typedef typename K::Point_2 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_2 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        { 
	  return to_double(k.compute_y_2_object()(p)) < value;
        }
    };
}


template <class K>
class Hilbert_sort_middle_2
{
public:
    typedef K Kernel;
    typedef typename Kernel::Point_2 Point;
    
private:
    Kernel _k;
    std::ptrdiff_t _limit;

    template <int x, bool up> struct Cmp : public internal::Fixed_hilbert_cmp_2<Kernel,x,up>
      { Cmp (double v, const Kernel &k) : internal::Fixed_hilbert_cmp_2<Kernel,x,up> (v, k) {} };

public:
    Hilbert_sort_middle_2 (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
        : _k(k), _limit (limit)
    {}

    template <int x, bool upx, bool upy, class RandomAccessIterator>
    void sort (RandomAccessIterator begin, RandomAccessIterator end,
	       double xmin, double ymin, double xmax, double ymax) const
    {
        const int y = (x + 1) % 2;
        if (end - begin <= _limit) return;

	double xmed= (xmin+xmax)/2;
	double ymed= (ymin+ymax)/2;

        RandomAccessIterator m0 = begin, m4 = end;

        RandomAccessIterator m2 = 
	  internal::fixed_hilbert_split (m0, m4, Cmp< x,  upx> (xmed,_k));
        RandomAccessIterator m1 = 
	  internal::fixed_hilbert_split (m0, m2, Cmp< y,  upy> (ymed,_k));
        RandomAccessIterator m3 = 
	  internal::fixed_hilbert_split (m2, m4, Cmp< y, !upy> (ymed,_k));

        sort<y, upy, upx> (m0, m1, ymin, xmin, ymed, xmed);
        sort<x, upx, upy> (m1, m2, xmin, ymed, xmed, ymax);
        sort<x, upx, upy> (m2, m3, xmed, ymed, xmax, ymax);
        sort<y,!upy,!upx> (m3, m4, ymed, xmax, ymin, xmed);
    }

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
      Bbox_2 box=bbox_2(begin, end);
      sort <0, false, false> (begin, end, 
			      box.xmin(), box.ymin(), box.xmax(), box.ymax());
    }
};


} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_2_H
