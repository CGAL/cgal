// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Olivier Devillers
//               : Pedro Machado Manhaes de Castro


#ifndef CGAL_SPATIAL_SORT_ON_SPHERE_H
#define CGAL_SPATIAL_SORT_ON_SPHERE_H

#include <CGAL/spatial_sort.h>
#include <CGAL/hilbert_sort_on_sphere.h>
#include <CGAL/Hilbert_policy_tags.h>

namespace CGAL {

namespace internal {

template <class RandomAccessIterator, class PolicyTag, class Kernel,
          class FT = typename Kernel::FT,
          class Point = typename Kernel::Point_3>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
                             const Kernel &k,
                             const Hilbert_policy<PolicyTag> /*policy*/,
                             const FT sq_r,
                             const Point &p,
                             std::ptrdiff_t threshold_hilbert,
                             std::ptrdiff_t threshold_multiscale,
                             double ratio)
{
  typedef Hilbert_sort_on_sphere_3<Kernel, Hilbert_policy<PolicyTag>, Point> Sort;
  typedef std::iterator_traits<RandomAccessIterator> ITraits;
  typedef typename ITraits::difference_type Diff_t;

  boost::rand48 random;
  boost::random_number_generator<boost::rand48, Diff_t> rng(random);
  CGAL::cpp98::random_shuffle(begin, end, rng);

  if (threshold_hilbert==0) threshold_hilbert=4;
  if (threshold_multiscale==0) threshold_multiscale=16;
  if (ratio==0.0) ratio=0.25;

  (Multiscale_sort<Sort> (Sort (k, sq_r, p, threshold_hilbert),
                          threshold_multiscale, ratio)) (begin, end);
}

} // end of namespace internal

template <class RandomAccessIterator, class PolicyTag,
          class Kernel = typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel,
          class FT = typename Kernel::FT,
          class Point = typename Kernel::Point_3>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
                             const Hilbert_policy<PolicyTag> policy,
                             const FT sq_r = FT(1),
                             const Point &p = Point(0,0,0),
                             std::ptrdiff_t threshold_hilbert = 0,
                             std::ptrdiff_t threshold_multiscale = 0,
                             const double ratio = 0.)
{
  internal::spatial_sort_on_sphere (begin, end, Kernel(), policy, sq_r, p,
                                    threshold_hilbert, threshold_multiscale, ratio);
}

template <class RandomAccessIterator, class Kernel,
          class FT = typename Kernel::FT,
          class Point = typename Kernel::Point_3>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
                             const Kernel &k,
                             const FT sq_r = FT(1),
                             const Point &p = Point(0,0,0),
                             std::ptrdiff_t threshold_hilbert = 0,
                             std::ptrdiff_t threshold_multiscale = 0,
                             const double ratio = 0.)
{
  internal::spatial_sort_on_sphere (begin, end, k,
                                    Hilbert_sort_median_policy(), sq_r, p,
                                    threshold_hilbert, threshold_multiscale, ratio);
}

template <class RandomAccessIterator,
          class Kernel = typename CGAL::Kernel_traits<typename std::iterator_traits<RandomAccessIterator>::value_type>::Kernel,
          class FT = typename Kernel::FT,
          class Point = typename Kernel::Point_3>
void spatial_sort_on_sphere (RandomAccessIterator begin, RandomAccessIterator end,
                             const FT sq_r = FT(1),
                             const Point &p = Point(0,0,0),
                             std::ptrdiff_t threshold_hilbert = 0,
                             std::ptrdiff_t threshold_multiscale = 0,
                             const double ratio = 0.)
{
  internal::spatial_sort_on_sphere (begin, end, Kernel(),
                                    Hilbert_sort_median_policy(), sq_r, p,
                                    threshold_hilbert, threshold_multiscale, ratio);
}

} // end of namespace CGAL

#endif // CGAL_SPATIAL_SORT_ON_SPHERE_H
