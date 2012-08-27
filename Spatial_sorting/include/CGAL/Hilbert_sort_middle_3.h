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

#ifndef CGAL_HILBERT_SORT_MIDDLE_3_H
#define CGAL_HILBERT_SORT_MIDDLE_3_H

#include <CGAL/basic.h>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_middle_base.h>

namespace CGAL {

namespace internal {
    template <class K, int x, bool up> struct Fixed_hilbert_cmp_3;

    template <class K, int x>
    struct Fixed_hilbert_cmp_3<K,x,true>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
	  return ! Fixed_hilbert_cmp_3<K,x,false> (value,k) (p);
        }
    };

    template <class K>
    struct Fixed_hilbert_cmp_3<K,0,false>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
	  return to_double(k.compute_x_3_object()(p)) < value;
        }
    };

    template <class K>
    struct Fixed_hilbert_cmp_3<K,1,false>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
	  return to_double(k.compute_y_3_object()(p)) < value;
        }
    };

    template <class K>
    struct Fixed_hilbert_cmp_3<K,2,false>
        : public std::binary_function<typename K::Point_3,
                                      typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
	double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k = K()) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
	  return to_double(k.compute_z_3_object()(p)) < value ;
        }
    };
}

template <class K>
class Hilbert_sort_middle_3
{
public:
    typedef K Kernel;
    typedef typename Kernel::Point_3 Point;

private:
    Kernel _k;
    std::ptrdiff_t _limit;

    template <int x, bool up> struct Cmp : public internal::Fixed_hilbert_cmp_3<Kernel,x,up>
      { Cmp (double v,const Kernel &k) : internal::Fixed_hilbert_cmp_3<Kernel,x,up> (v,k) {} };

public:
    Hilbert_sort_middle_3 (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
        : _k(k), _limit (limit)
    {}

    template <int x, bool upx, bool upy, bool upz, class RandomAccessIterator>
      void sort (RandomAccessIterator begin, RandomAccessIterator end,
		 double xmin, double ymin, double zmin, 
		 double xmax, double ymax, double zmax) const
    {
        const int y = (x + 1) % 3, z = (x + 2) % 3;
        if (end - begin <= _limit) return;

	double xmed= (xmin+xmax)/2;
	double ymed= (ymin+ymax)/2;
	double zmed= (zmin+zmax)/2;


        RandomAccessIterator m0 = begin, m8 = end;

        RandomAccessIterator m4 = 
	  internal::fixed_hilbert_split (m0, m8, Cmp< x,  upx> (xmed,_k));
        RandomAccessIterator m2 = 
	  internal::fixed_hilbert_split (m0, m4, Cmp< y,  upy> (ymed,_k));
        RandomAccessIterator m6 = 
	  internal::fixed_hilbert_split (m4, m8, Cmp< y, !upy> (ymed,_k));
        RandomAccessIterator m1 = 
	  internal::fixed_hilbert_split (m0, m2, Cmp< z,  upz> (zmed,_k));
        RandomAccessIterator m3 = 
	  internal::fixed_hilbert_split (m2, m4, Cmp< z, !upz> (zmed,_k));
        RandomAccessIterator m5 = 
	  internal::fixed_hilbert_split (m4, m6, Cmp< z,  upz> (zmed,_k));
        RandomAccessIterator m7 = 
	  internal::fixed_hilbert_split (m6, m8, Cmp< z, !upz> (zmed,_k));


        sort<z, upz, upx, upy> (m0, m1, zmin, xmin, ymin, zmed, xmed, ymed);
        sort<y, upy, upz, upx> (m1, m2, ymin, zmed, xmin, ymed, zmax, xmed);
        sort<y, upy, upz, upx> (m2, m3, ymed, zmed, xmin, ymax, zmax, xmed);
        sort<x, upx,!upy,!upz> (m3, m4, xmin, ymax, zmed, xmed, ymed, zmin);
        sort<x, upx,!upy,!upz> (m4, m5, xmed, ymax, zmed, xmax, ymed, zmin);
        sort<y,!upy, upz,!upx> (m5, m6, ymax, zmed, xmax, ymed, zmax, xmed);
        sort<y,!upy, upz,!upx> (m6, m7, ymed, zmed, xmax, ymin, zmax, xmed);
        sort<z,!upz,!upx, upy> (m7, m8, zmed, xmax, ymin, zmin, xmed, ymed);
    }

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
      K k;
      double xmin=to_double(k.compute_x_3_object()(*begin)),
	     ymin=to_double(k.compute_y_3_object()(*begin)),
	     zmin=to_double(k.compute_z_3_object()(*begin)),
  	     xmax=xmin,
	     ymax=ymin,
	     zmax=zmin;
      for(RandomAccessIterator it=begin+1; it<end; ++it){
	if ( to_double(k.compute_x_3_object()(*it)) < xmin) 
	  xmin = to_double(k.compute_x_3_object()(*it));
	if ( to_double(k.compute_y_3_object()(*it)) < ymin) 
	  ymin = to_double(k.compute_y_3_object()(*it));
	if ( to_double(k.compute_z_3_object()(*it)) < zmin) 
	  zmin = to_double(k.compute_z_3_object()(*it));
	if ( to_double(k.compute_x_3_object()(*it)) > xmax) 
	  xmax = to_double(k.compute_x_3_object()(*it));
	if ( to_double(k.compute_y_3_object()(*it)) > ymax) 
	  ymax = to_double(k.compute_y_3_object()(*it));
	if ( to_double(k.compute_z_3_object()(*it)) > zmax) 
	  zmax = to_double(k.compute_z_3_object()(*it));
      }
      
      sort <0, false, false, false> (begin, end, xmin,ymin,zmin,xmax,ymax,zmax);
    }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_3_H
