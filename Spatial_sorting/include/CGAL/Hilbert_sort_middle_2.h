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

#ifndef CGAL_HILBERT_SORT_MIDDLE_2_H
#define CGAL_HILBERT_SORT_MIDDLE_2_H

#include <CGAL/config.h>
#include <functional>
#include <cstddef>
#include <utility>
#include <CGAL/Hilbert_sort_middle_base.h>
#include <CGAL/FPU.h>
#include <CGAL/number_utils.h>
#include <CGAL/Spatial_sorting/internal/Spatial_sorting_helper.h>

namespace CGAL {

namespace internal {
template <class K, int x, bool up> struct Fixed_hilbert_cmp_2;

template <class K, int x>
struct Fixed_hilbert_cmp_2<K,x,true>
    : public CGAL::cpp98::binary_function<typename K::Point_2,
    typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  double value;
  Fixed_hilbert_cmp_2 (double v, const K &_k) : k(_k),value(v) {}
  bool operator() (const Point &p) const
  {
    return ! Fixed_hilbert_cmp_2<K,x,false> (value, k) (p);
  }
};

template <class K>
struct Fixed_hilbert_cmp_2<K,0,false>
    : public CGAL::cpp98::binary_function<typename K::Point_2,
    typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  double value;
  Fixed_hilbert_cmp_2 (double v, const K &_k) : k(_k),value(v) {}
  bool operator() (const Point &p) const
  {
    return to_double(k.compute_x_2_object()(p)) < value;
  }
};

template <class K>
struct Fixed_hilbert_cmp_2<K,1,false>
    : public CGAL::cpp98::binary_function<typename K::Point_2,
    typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  double value;
  Fixed_hilbert_cmp_2 (double v, const K &_k) : k(_k),value(v) {}
  bool operator() (const Point &p) const
  {
    return to_double(k.compute_y_2_object()(p)) < value;
  }
};

} // namespace internal

template <class K>
class Hilbert_sort_middle_2
{
public:
  typedef K Kernel;
  typedef typename Kernel::Point_2 Point;

private:
  Kernel _k;
  std::ptrdiff_t _limit;

  template <int x, bool up>
  struct Cmp
    : public internal::Fixed_hilbert_cmp_2<Kernel,x,up>
  {
    Cmp (double v, const Kernel &k) : internal::Fixed_hilbert_cmp_2<Kernel,x,up> (v, k) {}
  };

public:
  Hilbert_sort_middle_2 (const Kernel &k, std::ptrdiff_t limit = 1)
    : _k(k), _limit (limit)
  {}

  // Hilbert sort use double representations so two points not strictly equal will treat identically by hilbert sort
  bool have_equal_representations(Point a, Point b) const{
    return (to_double(_k.compute_x_2_object()(a)) == to_double(_k.compute_x_2_object()(b))) &&
           (to_double(_k.compute_y_2_object()(a)) == to_double(_k.compute_y_2_object()(b)));
  }

  template <class RandomAccessIterator>
  void dispatch_sort (int x, bool upx, bool upy,
                      RandomAccessIterator begin, RandomAccessIterator end,
                      double xmin, double ymin, double xmax, double ymax) const
  {
    switch ((x << 2) | (int(upx) << 1) | int(upy)) {
      case 0:  sort<0, false, false>(begin, end, xmin, ymin, xmax, ymax); break;
      case 1:  sort<0, false, true >(begin, end, xmin, ymin, xmax, ymax); break;
      case 2:  sort<0, true,  false>(begin, end, xmin, ymin, xmax, ymax); break;
      case 3:  sort<0, true,  true >(begin, end, xmin, ymin, xmax, ymax); break;
      case 4:  sort<1, false, false>(begin, end, xmin, ymin, xmax, ymax); break;
      case 5:  sort<1, false, true >(begin, end, xmin, ymin, xmax, ymax); break;
      case 6:  sort<1, true,  false>(begin, end, xmin, ymin, xmax, ymax); break;
      default: sort<1, true,  true >(begin, end, xmin, ymin, xmax, ymax); break;
    }
  }

  // Descends through the cells containing all the points without recursion:
  // reaching a cell of denormal size can take ~1075 levels.
  template <class RandomAccessIterator>
  void skip_empty_cells (int x, bool upx, bool upy,
                         RandomAccessIterator begin, RandomAccessIterator end,
                         double xmin, double ymin, double xmax, double ymax, const Bbox_2& bb) const
  {
    const double lo[2] = { bb.xmin(), bb.ymin() };
    const double hi[2] = { bb.xmax(), bb.ymax() };
    for(;;)
    {
      const int y = (x + 1) % 2;
      const double xmed= CGAL_IA_FORCE_TO_DOUBLE((xmin+xmax)/2);
      const double ymed= CGAL_IA_FORCE_TO_DOUBLE((ymin+ymax)/2);
      const bool x_degen = (xmin == xmed || xmed == xmax);
      const bool y_degen = (ymin == ymed || ymed == ymax);
      const bool x_lower = xmed <= lo[x], x_upper = xmed >= hi[x];
      const bool y_lower = ymed <= lo[y], y_upper = ymed >= hi[y];

      if(x_degen && y_degen)
        return;

      if( (!x_degen && !x_lower && !x_upper) || (!y_degen && !y_lower && !y_upper)){
        // Not contained in one quarter
        dispatch_sort(x, upx, upy, begin, end, xmin, ymin, xmax, ymax);
        return;
      }

      CGAL_assertion(!x_upper || !x_lower || lo[x] == hi[x]);
      CGAL_assertion(!y_upper || !y_lower || lo[y] == hi[y]);

      if( x_lower == upx )
      {
        if( y_lower == upy ) {          // first quarter
          std::swap(xmin, ymin);
          xmax = ymed;
          ymax = xmed;
          std::swap(upx, upy);
          x = y;
        } else {                        // second quarter
          ymin = ymed;
          xmax = xmed;
        }
      }
      else
      {
        if( y_lower == upy ) {          // fourth quarter
          xmin = ymed;
          std::swap(xmax, ymin);
          ymax = xmed;
          std::swap(upx, upy);
          upx = !upx;
          upy = !upy;
          x = y;
        } else {                        // third quarter
          xmin = xmed;
          ymin = ymed;
        }
      }
    }
  }

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end,
             double xmin, double ymin, double xmax, double ymax) const
  {
    const int y = (x + 1) % 2;
    if (end - begin <= _limit) return;

    double xmed= CGAL_IA_FORCE_TO_DOUBLE((xmin+xmax)/2);
    double ymed= CGAL_IA_FORCE_TO_DOUBLE((ymin+ymax)/2);

    RandomAccessIterator m0 = begin, m4 = end;

    RandomAccessIterator m2 = internal::fixed_hilbert_split (m0, m4, Cmp< x,  upx> (xmed,_k));
    RandomAccessIterator m1 = internal::fixed_hilbert_split (m0, m2, Cmp< y,  upy> (ymed,_k));
    RandomAccessIterator m3 = internal::fixed_hilbert_split (m2, m4, Cmp< y, !upy> (ymed,_k));

    if( ((m2 == m4) && ((m1 == m0) || (m1 == m2))) ||
        ((m2 == m0) && ((m3 == m2) || (m3 == m4))) )
    {
      // All elements are in one quarter
      bool all_equal = true;
      for(auto it=begin; it!=end-1; ++it)
        if(!have_equal_representations(*it, *(it+1))){
          all_equal = false;
          break;
        }
      if(all_equal) return;
      Bbox_2 bb=bbox_2(begin, end, Spatial_sorting::internal::Kernel_wrapper_2<Kernel>(_k));
      skip_empty_cells (x, upx, upy, begin, end, xmin, ymin, xmax, ymax, bb);
    }
    else
    {
      sort<y, upy, upx> (m0, m1, ymin, xmin, ymed, xmed);
      sort<x, upx, upy> (m1, m2, xmin, ymed, xmed, ymax);
      sort<x, upx, upy> (m2, m3, xmed, ymed, xmax, ymax);
      sort<y,!upy,!upx> (m3, m4, ymed, xmax, ymin, xmed);
    }
  }

  template <class RandomAccessIterator>
  void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    if (begin == end) return;
    Bbox_2 box=bbox_2(begin, end, Spatial_sorting::internal::Kernel_wrapper_2<Kernel>(_k));
    sort <0, false, false> (begin, end, box.xmin(), box.ymin(), box.xmax(), box.ymax());
  }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_2_H
