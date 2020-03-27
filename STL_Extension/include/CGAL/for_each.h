// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_FOR_EACH_H
#define CGAL_FOR_EACH_H

#include <CGAL/iterator.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#endif // CGAL_LINKED_WITH_TBB

/*
  CGAL::for_each<ConcurrencyTag>(Range, Function) does the following:

  - if Sequential_tag is used, apply Function to all elements of Range
  - if Parallel_tag is used:
    * static assert that TBB is available
    * if Range has random access iterators, use a TBB parallel_for
    loop to apply Function to all elements of Range
    * if Range doesn't have random access iterators, first copy the
    iterators in a vector and then use a TBB parallel_for loop to
    apply Function to all elements of Range

  The loop is interrupted if `functor` returns false (it carries on
  until the end otherwise).
*/

namespace CGAL {

namespace internal {

template <typename RangeRef, typename IteratorCategory>
void for_each (RangeRef range,
               const std::function<bool(typename std::iterator_traits
                                        <typename Range_iterator_type<RangeRef>::type>::reference)>& functor,
               const Sequential_tag&,
               IteratorCategory)
{
  for (typename std::iterator_traits
         <typename Range_iterator_type<RangeRef>::type>::reference r : range)
    if (!functor(r))
      break;
}

#ifdef CGAL_LINKED_WITH_TBB
template <typename RangeRef, typename Fct, typename IteratorCategory>
void for_each (RangeRef range,
               const Fct& functor,
               const Parallel_tag&,
               IteratorCategory)
{
  std::size_t range_size = std::distance (range.begin(), range.end());

  std::vector<typename Range_iterator_type<RangeRef>::type> iterators;
  iterators.reserve (range_size);
  for (typename Range_iterator_type<RangeRef>::type it = range.begin(); it != range.end(); ++ it)
    iterators.push_back (it);

  tbb::parallel_for (tbb::blocked_range<std::size_t>(0, range_size),
                     [&](const tbb::blocked_range<std::size_t>& r)
                     {
                       for (std::size_t i = r.begin(); i != r.end(); ++ i)
                         if (!functor (*(iterators[i])))
                           break;
                     });
}

template <typename RangeRef, typename Fct>
void for_each (RangeRef range,
               const Fct& functor,
               const Parallel_tag&,
               std::random_access_iterator_tag)
{
  std::size_t range_size = std::distance (range.begin(), range.end());

  tbb::parallel_for (tbb::blocked_range<std::size_t>(0, range_size),
                     [&](const tbb::blocked_range<std::size_t>& r)
                     {
                       for (std::size_t i = r.begin(); i != r.end(); ++ i)
                         if (!functor (*(range.begin() + i)))
                           break;
                     });
}
#endif

} // namespace internal

template <typename ConcurrencyTag, typename Range>
void for_each (const Range& range,
               const std::function<bool(typename std::iterator_traits
                                        <typename Range::const_iterator>::reference)>& functor)
{
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#endif

  internal::for_each<const Range&>
    (range, functor,
     ConcurrencyTag(),
     typename std::iterator_traits<typename Range::const_iterator>::iterator_category());
}

template <typename ConcurrencyTag, typename Range>
void for_each (Range& range,
               const std::function<bool(typename std::iterator_traits
                                        <typename Range::iterator>::reference)>& functor)
{
#ifndef CGAL_LINKED_WITH_TBB
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#endif

  internal::for_each<Range&>
    (range, functor,
     ConcurrencyTag(),
     typename std::iterator_traits<typename Range::iterator>::iterator_category());
}

} // namespace CGAL


#endif // CGAL_FOR_EACH_H
