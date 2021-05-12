// Copyright (c) 2007-2011  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Christophe Delage

#ifndef CGAL_SPATIAL_SORT_H
#define CGAL_SPATIAL_SORT_H

#include <CGAL/config.h>

#include <CGAL/hilbert_sort.h>
#include <CGAL/Multiscale_sort.h>

#include <boost/random/random_number_generator.hpp>
#include <CGAL/algorithm.h>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>

#include <algorithm>

namespace CGAL {


namespace internal {

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
                   Policy /*policy*/,
                   typename Kernel::Point_2 *,
                   std::ptrdiff_t threshold_hilbert,
                   std::ptrdiff_t threshold_multiscale,
                   double ratio)
{
  typedef std::iterator_traits<RandomAccessIterator> Iterator_traits;
  typedef typename Iterator_traits::difference_type Diff_t;
  typedef Hilbert_sort_2<Kernel, Policy, ConcurrencyTag> Sort;
  boost::rand48 random;
  boost::random_number_generator<boost::rand48, Diff_t> rng(random);
  CGAL::cpp98::random_shuffle(begin,end,rng);

  if (threshold_hilbert==0) threshold_hilbert=4;
  if (threshold_multiscale==0) threshold_multiscale=16;
  if (ratio==0.0) ratio=0.25;

  (Multiscale_sort<Sort> (Sort (k, threshold_hilbert), threshold_multiscale, ratio)) (begin, end);
}

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
                   Policy /*policy*/,
                   typename Kernel::Point_3 *,
                   std::ptrdiff_t threshold_hilbert,
                   std::ptrdiff_t threshold_multiscale,
                   double ratio)
{
  typedef std::iterator_traits<RandomAccessIterator> Iterator_traits;
  typedef typename Iterator_traits::difference_type Diff_t;
  typedef Hilbert_sort_3<Kernel, Policy, ConcurrencyTag> Sort;
  boost::rand48 random;
  boost::random_number_generator<boost::rand48, Diff_t> rng(random);
  CGAL::cpp98::random_shuffle(begin,end, rng);

  if (threshold_hilbert==0) threshold_hilbert=8;
  if (threshold_multiscale==0) threshold_multiscale=64;
  if (ratio==0.0) ratio=0.125;

  (Multiscale_sort<Sort> (Sort (k, threshold_hilbert), threshold_multiscale, ratio)) (begin, end);
}

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
                   Policy /*policy*/,
                   typename Kernel::Point_d *,
                   std::ptrdiff_t threshold_hilbert,
                   std::ptrdiff_t threshold_multiscale,
                   double ratio)
{
  typedef std::iterator_traits<RandomAccessIterator> Iterator_traits;
  typedef typename Iterator_traits::difference_type Diff_t;
  typedef Hilbert_sort_d<Kernel, Policy> Sort;
  boost::rand48 random;
  boost::random_number_generator<boost::rand48, Diff_t> rng(random);
  CGAL::cpp98::random_shuffle(begin,end, rng);

  if (threshold_hilbert==0) threshold_hilbert=10;
  if (threshold_multiscale==0) threshold_multiscale=500;
  if (ratio==0.0) ratio=0.05;

  (Multiscale_sort<Sort> (Sort (k, threshold_hilbert), threshold_multiscale, ratio)) (begin, end);
}

} //namespace internal

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Policy, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
                   Policy policy,
                   std::ptrdiff_t threshold_hilbert=0,
                   std::ptrdiff_t threshold_multiscale=0,
                   double ratio=0.0)
{
  typedef std::iterator_traits<RandomAccessIterator> ITraits;
  typedef typename ITraits::value_type               value_type;

  internal::spatial_sort<ConcurrencyTag>(begin, end, k, policy, static_cast<value_type *> (nullptr),
                                         threshold_hilbert,threshold_multiscale,ratio);
}

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   Hilbert_sort_median_policy policy,
                   std::ptrdiff_t threshold_hilbert=0,
                   std::ptrdiff_t threshold_multiscale=0,
                   double ratio=0.0)
{
  typedef std::iterator_traits<RandomAccessIterator> ITraits;
  typedef typename ITraits::value_type               value_type;
  typedef CGAL::Kernel_traits<value_type>            KTraits;
  typedef typename KTraits::Kernel                   Kernel;

  spatial_sort<ConcurrencyTag> (begin, end, Kernel(), policy,
                                threshold_hilbert,threshold_multiscale,ratio);
}

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   Hilbert_sort_middle_policy policy,
                   std::ptrdiff_t threshold_hilbert=0,
                   std::ptrdiff_t threshold_multiscale=0,
                   double ratio=0.0)
{
  typedef std::iterator_traits<RandomAccessIterator> ITraits;
  typedef typename ITraits::value_type               value_type;
  typedef CGAL::Kernel_traits<value_type>            KTraits;
  typedef typename KTraits::Kernel                   Kernel;

  spatial_sort<ConcurrencyTag> (begin, end, Kernel(), policy,
                                threshold_hilbert,threshold_multiscale,ratio);
}


template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator, class Kernel>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   const Kernel &k,
                   std::ptrdiff_t threshold_hilbert=0,
                   std::ptrdiff_t threshold_multiscale=0,
                   double ratio=0.0)
{
  spatial_sort<ConcurrencyTag> (begin, end, k,
                                Hilbert_sort_median_policy(),
                                threshold_hilbert,threshold_multiscale,ratio);
}

template <class ConcurrencyTag = Sequential_tag, class RandomAccessIterator>
void spatial_sort (RandomAccessIterator begin, RandomAccessIterator end,
                   std::ptrdiff_t threshold_hilbert=0,
                   std::ptrdiff_t threshold_multiscale=0,
                   double ratio=0.0)
{
  spatial_sort<ConcurrencyTag> (begin, end,
                                Hilbert_sort_median_policy(),
                                threshold_hilbert,threshold_multiscale,ratio);
}

} // namespace CGAL

#endif//CGAL_SPATIAL_SORT_H
