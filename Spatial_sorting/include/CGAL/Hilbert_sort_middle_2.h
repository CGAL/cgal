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
#include <CGAL/Hilbert_sort_middle_base.h>
#include <CGAL/number_utils.h>

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

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void skip_empty_cells (RandomAccessIterator begin, RandomAccessIterator end,
                         double xmin, double ymin, double xmax, double ymax, Bbox_2 bb) const
  {
    const int y = (x + 1) % 2;
    double xmed= (xmin+xmax)/2;
    double ymed= (ymin+ymax)/2;
    const bool x_lower = (x==0) ? xmed<bb.xmin() : xmed<bb.ymin();
    const bool x_upper = (x==0) ? xmed>bb.xmax() : xmed>bb.ymax();
    const bool y_lower = (x==0) ? ymed<bb.ymin() : ymed<bb.xmin();
    const bool y_upper = (x==0) ? ymed>bb.ymax() : ymed>bb.xmax();
    if( (!x_lower && !x_upper) || (!y_lower && !y_upper)){
      // Not contain in one quarter
      sort<x, upx, upy> (begin, end, xmin, ymin, xmax, ymax);
      return;
    }

    if( x_upper != upx )
    {
      if( y_upper != upy )
        skip_empty_cells<y, upy, upx> (begin, end, ymin, xmin, ymed, xmed, bb);
      else
        skip_empty_cells<x, upx, upy> (begin, end, xmin, ymed, xmed, ymax, bb);
    }
    else
    {
      if( y_upper != upy )
        skip_empty_cells<y, !upx, !upy> (begin, end, ymed, xmax, ymin, xmed, bb);
      else
        skip_empty_cells<x, upx, upy> (begin, end, xmed, ymed, xmax, ymax, bb);
    }
  }

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end,
             double xmin, double ymin, double xmax, double ymax) const
  {
    const int y = (x + 1) % 2;
    if (end - begin <= _limit) return;

    double xmed= (xmin+xmax)/2;
    double ymed= (ymin+ymax)/2;

    RandomAccessIterator m0 = begin, m4 = end;

    RandomAccessIterator m2 = internal::fixed_hilbert_split (m0, m4, Cmp< x,  upx> (xmed,_k));
    RandomAccessIterator m1 = internal::fixed_hilbert_split (m0, m2, Cmp< y,  upy> (ymed,_k));
    RandomAccessIterator m3 = internal::fixed_hilbert_split (m2, m4, Cmp< y, !upy> (ymed,_k));

    if( ((m2 == m4) && ((m1 == m0) || (m1 == m2))) ||
        ((m2 == m0) && ((m3 == m2) || (m3 == m4))) )
    {
      // All elements are in one quarter
      if(_k.equal_2_object()(*begin, *end))
        return; // Return if all elements are equal
      Bbox_2 bb=bbox_2(begin, end, _k);
      if constexpr(x == 0)
        skip_empty_cells<x, upx, upy> (begin, end, xmin, ymin, xmax, ymax, bb);
      else
        skip_empty_cells<x, upx, upy> (begin, end, xmin, ymin, xmax, ymax, bb);
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
    Bbox_2 box=bbox_2(begin, end, _k);
    sort <0, false, false> (begin, end, box.xmin(), box.ymin(), box.xmax(), box.ymax());
  }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MIDDLE_2_H
