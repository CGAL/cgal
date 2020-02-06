// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Christophe Delage

#ifndef CGAL_HILBERT_SORT_MEDIAN_2_H
#define CGAL_HILBERT_SORT_MEDIAN_2_H

#include <CGAL/assertions.h>
#include <CGAL/config.h>
#include <CGAL/tags.h>
#include <CGAL/Hilbert_sort_base.h>

#include <boost/type_traits/is_convertible.hpp>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_invoke.h>
#endif

#include <functional>
#include <cstddef>

namespace CGAL {

namespace internal {

template <class K, int x, bool up> struct Hilbert_cmp_2;

template <class K, int x>
struct Hilbert_cmp_2<K,x,true>
  : public CGAL::cpp98::binary_function<typename K::Point_2, typename K::Point_2, bool>
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
    : public CGAL::cpp98::binary_function<typename K::Point_2, typename K::Point_2, bool>
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
  : public CGAL::cpp98::binary_function<typename K::Point_2, typename K::Point_2, bool>
{
  typedef typename K::Point_2 Point;
  K k;
  Hilbert_cmp_2 (const K &_k = K()) : k(_k) {}
  bool operator() (const Point &p, const Point &q) const
  {
    return k.less_y_2_object() (p, q);
  }
};

} // namespace internal

template <class K, class ConcurrencyTag>
class Hilbert_sort_median_2
{
public:
  typedef Hilbert_sort_median_2<K, ConcurrencyTag> Self;
  typedef K Kernel;
  typedef typename Kernel::Point_2 Point;

private:
  Kernel _k;
  std::ptrdiff_t _limit;

  template <int x, bool up>
  struct Cmp
    : public internal::Hilbert_cmp_2<Kernel,x,up>
  {
    Cmp (const Kernel &k) : internal::Hilbert_cmp_2<Kernel,x,up> (k) {}
  };

public:
  Hilbert_sort_median_2 (const Kernel &k = Kernel(), std::ptrdiff_t limit = 1)
    : _k(k), _limit (limit)
  {}

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void recursive_sort (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    const int y = (x + 1) % 2;
    if (end - begin <= _limit)
      return;

    RandomAccessIterator m0 = begin, m4 = end;

    RandomAccessIterator m2 = internal::hilbert_split (m0, m4, Cmp< x,  upx> (_k));
    RandomAccessIterator m1 = internal::hilbert_split (m0, m2, Cmp< y,  upy> (_k));
    RandomAccessIterator m3 = internal::hilbert_split (m2, m4, Cmp< y, !upy> (_k));

    recursive_sort<y, upy, upx> (m0, m1);
    recursive_sort<x, upx, upy> (m1, m2);
    recursive_sort<x, upx, upy> (m2, m3);
    recursive_sort<y,!upy,!upx> (m3, m4);
  }

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  struct Recursive_sort
  {
    const Self& hs;
    RandomAccessIterator begin,end;

    Recursive_sort(const Self& hs, RandomAccessIterator begin, RandomAccessIterator end)
      : hs(hs), begin(begin), end(end)
    {}

    void operator()() const
    {
      hs.recursive_sort<x,upx,upy>(begin,end);
    }
  };

  template <class RandomAccessIterator, class Comp>
  struct Hilbert_split
  {
    RandomAccessIterator& res;
    RandomAccessIterator begin, end;
    const Comp& comp;

    Hilbert_split(RandomAccessIterator& res,
                  RandomAccessIterator begin, RandomAccessIterator end,
                  const Comp& comp)
      : res(res), begin(begin), end(end), comp(comp)
    {}

    void operator()() const
    {
      res = internal::hilbert_split(begin, end, comp);
    }
  };

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end, Parallel_tag) const
  {
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_USE(begin);
    CGAL_USE(end);
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    const int y = (x + 1) % 2;
    if (std::distance(begin,end) <= _limit)
      return;

    RandomAccessIterator m0 = begin, m4 = end;

    if(std::distance(begin,end) > 8192){ // 2^13, empirically a good cutoff
      RandomAccessIterator m1, m2, m3;
      m2 = internal::hilbert_split (m0, m4, Cmp< x,  upx> (_k));

      tbb::parallel_invoke(Hilbert_split<RandomAccessIterator, Cmp< y,  upy> > (m1, m0, m2, Cmp< y,  upy> (_k)),
                           Hilbert_split<RandomAccessIterator, Cmp< y, !upy> > (m3, m2, m4, Cmp< y, !upy> (_k)));

      tbb::parallel_invoke(Recursive_sort<y, upy, upx, RandomAccessIterator> (*this, m0, m1),
                           Recursive_sort<x, upx, upy, RandomAccessIterator> (*this, m1, m2),
                           Recursive_sort<x, upx, upy, RandomAccessIterator> (*this, m2, m3),
                           Recursive_sort<y,!upy,!upx, RandomAccessIterator> (*this, m3, m4));
    } else {
      recursive_sort<0, false, false>(begin, end);
    }
#endif
  }

  template <int x, bool upx, bool upy, class RandomAccessIterator>
  void sort (RandomAccessIterator begin, RandomAccessIterator end, Sequential_tag) const
  {
    recursive_sort<0, false, false>(begin, end);
  }

  template <class RandomAccessIterator>
  void operator() (RandomAccessIterator begin, RandomAccessIterator end) const
  {
    sort <0, false, false> (begin, end, ConcurrencyTag());
  }
};

} // namespace CGAL

#endif//CGAL_HILBERT_SORT_MEDIAN_2_H
