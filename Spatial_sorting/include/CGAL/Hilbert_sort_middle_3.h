// Copyright (c) 2011  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :  Olivier Devillers

#ifndef CGAL_HILBERT_SORT_MIDDLE_3_H
#define CGAL_HILBERT_SORT_MIDDLE_3_H

#include <CGAL/config.h>
#include <cmath>
#include <functional>
#include <cstddef>
#include <CGAL/Hilbert_sort_middle_base.h>
#include <CGAL/Spatial_sorting/internal/Spatial_sorting_helper.h>

namespace CGAL {

namespace internal {
    template <class K, int x, bool up> struct Fixed_hilbert_cmp_3;

    template <class K, int x>
    struct Fixed_hilbert_cmp_3<K,x,true>
        : public CGAL::cpp98::binary_function<typename K::Point_3,
                                              typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
          return ! Fixed_hilbert_cmp_3<K,x,false> (value,k) (p);
        }
    };

    template <class K>
    struct Fixed_hilbert_cmp_3<K,0,false>
        : public CGAL::cpp98::binary_function<typename K::Point_3,
                                              typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
          return to_double(k.compute_x_3_object()(p)) < value;
        }
    };

    template <class K>
    struct Fixed_hilbert_cmp_3<K,1,false>
        : public CGAL::cpp98::binary_function<typename K::Point_3,
                                              typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k) : k(_k),value(v) {}
        bool operator() (const Point &p) const
        {
          return to_double(k.compute_y_3_object()(p)) < value;
        }
    };

    template <class K>
    struct Fixed_hilbert_cmp_3<K,2,false>
        : public CGAL::cpp98::binary_function<typename K::Point_3,
                                              typename K::Point_3, bool>
    {
        typedef typename K::Point_3 Point;
        K k;
        double value;
        Fixed_hilbert_cmp_3 (double v, const K &_k) : k(_k),value(v) {}
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
    Hilbert_sort_middle_3 (const Kernel &k, std::ptrdiff_t limit = 1)
      : _k(k), _limit (limit)
    {}

    // Hilbert sort use double representations so two points not strictly equal will treat identically by hilbert sort
    bool have_equal_representations(Point a, Point b) const{
            return to_double(_k.compute_x_3_object()(a)) == to_double(_k.compute_x_3_object()(b)) &&
                   to_double(_k.compute_y_3_object()(a)) == to_double(_k.compute_y_3_object()(b)) &&
                   to_double(_k.compute_z_3_object()(a)) == to_double(_k.compute_z_3_object()(b));
    }

    template <int x, bool upx, bool upy, bool upz, class RandomAccessIterator>
    void skip_empty_cells (RandomAccessIterator begin, RandomAccessIterator end,
                           double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, Bbox_3 bb) const
    {
      const int y = (x + 1) % 3, z = (x + 2) % 3;
      double xmed= (xmin+xmax)/2;
      double ymed= (ymin+ymax)/2;
      double zmed= (zmin+zmax)/2;

      const bool x_degen = (xmin == xmed || xmed == xmax);
      const bool y_degen = (ymin == ymed || ymed == ymax);
      const bool z_degen = (zmin == zmed || zmed == zmax);

      if(x_degen && y_degen && z_degen)
        return;

      const bool x_lower = (x==0) ? xmed<=bb.xmin() : ((x==1) ? xmed<=bb.ymin() : xmed<=bb.zmin());
      const bool x_upper = (x==0) ? xmed>=bb.xmax() : ((x==1) ? xmed>=bb.ymax() : xmed>=bb.zmax());
      const bool y_lower = (x==0) ? ymed<=bb.ymin() : ((x==1) ? ymed<=bb.zmin() : ymed<=bb.xmin());
      const bool y_upper = (x==0) ? ymed>=bb.ymax() : ((x==1) ? ymed>=bb.zmax() : ymed>=bb.xmax());
      const bool z_lower = (x==0) ? zmed<=bb.zmin() : ((x==1) ? zmed<=bb.xmin() : zmed<=bb.ymin());
      const bool z_upper = (x==0) ? zmed>=bb.zmax() : ((x==1) ? zmed>=bb.xmax() : zmed>=bb.ymax());

      if( (!x_degen && !x_lower && !x_upper) || (!y_degen && !y_lower && !y_upper) || (!z_degen && !z_lower && !z_upper) ){
        // Not contain in one quarter
        sort<x, upx, upy, upz> (begin, end, xmin, ymin, zmin, xmax, ymax, zmax);
        return;
      }

      CGAL_assertion(!x_upper || !x_lower || (x==0 ? (bb.xmax() == bb.xmin()) : x==1 ? (bb.ymax() == bb.ymin()) : (bb.zmax() == bb.zmin())));
      CGAL_assertion(!y_upper || !y_lower || (y==0 ? (bb.xmax() == bb.xmin()) : y==1 ? (bb.ymax() == bb.ymin()) : (bb.zmax() == bb.zmin())));
      CGAL_assertion(!z_upper || !z_lower || (z==0 ? (bb.xmax() == bb.xmin()) : z==1 ? (bb.ymax() == bb.ymin()) : (bb.zmax() == bb.zmin())));

      if(x_lower == upx)
      {
        if(y_lower == upy)
        {
          if(z_lower == upz)
            skip_empty_cells<z, upz, upx, upy> (begin, end, zmin, xmin, ymin, zmed, xmed, ymed, bb);
          else
            skip_empty_cells<y, upy, upz, upx> (begin, end, ymin, zmed, xmin, ymed, zmax, xmed, bb);
        }
        else
        {
          if(z_lower == upz)
            skip_empty_cells<x, upx, !upy, !upz> (begin, end, xmin, ymax, zmed, xmed, ymed, zmin, bb);
          else
            skip_empty_cells<y, upy, upz, upx> (begin, end, ymed, zmed, xmin, ymax, zmax, xmed, bb);
        }
      }
      else
      {
        if(y_upper != upy)
        {
          if(z_upper != upz)
            skip_empty_cells<z, !upz, !upx, upy> (begin, end, zmed, xmax, ymin, zmin, xmed, ymed, bb);
          else
            skip_empty_cells<y, !upy, upz, !upx> (begin, end, ymed, zmed, xmax, ymin, zmax, xmed, bb);
        }
        else
        {
          if(z_upper != upz)
            skip_empty_cells<x, upx, !upy, !upz> (begin, end, xmed, ymax, zmed, xmax, ymed, zmin, bb);
          else
            skip_empty_cells<y, !upy, upz, !upx> (begin, end, ymax, zmed, xmax, ymed, zmax, xmed, bb);
        }
      }
    }

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

        bool all_contained_in_one_corner = true;
        if(m1 != m0 && m1 != m8) all_contained_in_one_corner = false;
        if(m2 != m0 && m2 != m8) all_contained_in_one_corner = false;
        if(m3 != m0 && m3 != m8) all_contained_in_one_corner = false;
        if(m4 != m0 && m4 != m8) all_contained_in_one_corner = false;
        if(m5 != m0 && m5 != m8) all_contained_in_one_corner = false;
        if(m6 != m0 && m6 != m8) all_contained_in_one_corner = false;
        if(m7 != m0 && m7 != m8) all_contained_in_one_corner = false;

        if(all_contained_in_one_corner){
          // If all elements are equal, we can return.
          bool all_equal = true;
          for(auto it=begin; it!=end-1; ++it)
            if(!have_equal_representations(*it, *(it+1))){
              all_equal = false;
              break;
            }
          if(all_equal) return;
          Bbox_3 bb=bbox_3(begin, end, Spatial_sorting::internal::Kernel_wrapper_3<Kernel>(_k));
          skip_empty_cells<x, upx, upy, upz> (begin, end, xmin, ymin, zmin, xmax, ymax, zmax, bb);
        }
        else
        {
          sort<z, upz, upx, upy> (m0, m1, zmin, xmin, ymin, zmed, xmed, ymed);
          sort<y, upy, upz, upx> (m1, m2, ymin, zmed, xmin, ymed, zmax, xmed);
          sort<y, upy, upz, upx> (m2, m3, ymed, zmed, xmin, ymax, zmax, xmed);
          sort<x, upx,!upy,!upz> (m3, m4, xmin, ymax, zmed, xmed, ymed, zmin);
          sort<x, upx,!upy,!upz> (m4, m5, xmed, ymax, zmed, xmax, ymed, zmin);
          sort<y,!upy, upz,!upx> (m5, m6, ymax, zmed, xmax, ymed, zmax, xmed);
          sort<y,!upy, upz,!upx> (m6, m7, ymed, zmed, xmax, ymin, zmax, xmed);
          sort<z,!upz,!upx, upy> (m7, m8, zmed, xmax, ymin, zmin, xmed, ymed);
        }
    }

    template <class RandomAccessIterator>
    void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
    {
      if (begin == end) return;
      Bbox_3 bb = bbox_3(begin, end, Spatial_sorting::internal::Kernel_wrapper_3<Kernel>(_k));
      sort <0, false, false, false> (begin, end, bb.xmin(),bb.ymin(),bb.zmin(),bb.xmax(),bb.ymax(),bb.zmax());
    }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_3_H
