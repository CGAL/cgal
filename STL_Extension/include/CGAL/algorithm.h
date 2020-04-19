// Copyright (c) 2003
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_ALGORITHM_H
#define CGAL_ALGORITHM_H

#include <CGAL/config.h>
#include <CGAL/utils.h>
#include <CGAL/enum.h>
#include <algorithm>
#include <iosfwd>
#include <iostream>

#include <boost/random/random_number_generator.hpp>
#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>

namespace CGAL {

// Not documented
template <class T> inline
bool
are_sorted(const T & a, const T & b, const T & c)
{
  return a <= b && b <= c;
}

// Not documented
template <class T, class Compare> inline
bool
are_sorted(const T & a, const T & b, const T & c, Compare cmp)
{
  return !cmp(b, a) && !cmp(c, b);
}

// Not documented
template <class T> inline
bool
are_strictly_sorted(const T & a, const T & b, const T & c)
{
  return a < b && b < c;
}

// Not documented
template <class T, class Compare> inline
bool
are_strictly_sorted(const T & a, const T & b, const T & c, Compare cmp)
{
  return cmp(a, b) && cmp(b, c);
}

// Not documented
// Checks that b is in the interval [min(a, c) , max(a, c)].
template <class T> inline
bool
are_ordered(const T & a, const T & b, const T & c)
{
  const T& min = (CGAL::min)(a, c);
  const T& max = (CGAL::max)(a, c);
  return min <= b && b <= max;
}

// Not documented
// Checks that b is in the interval [min(a, c) , max(a, c)].
template <class T, class Compare> inline
bool
are_ordered(const T & a, const T & b, const T & c, Compare cmp)
{
  const T& min = (std::min)(a, c, cmp);
  const T& max = (std::max)(a, c, cmp);
  return !cmp(b, min) && !cmp(max, b);
}

// Not documented
// Checks that b is in the interval ]min(a, c) , max(a, c)[.
template <class T> inline
bool
are_strictly_ordered(const T & a, const T & b, const T & c)
{
  const T& min = (CGAL::min)(a, c);
  const T& max = (CGAL::max)(a, c);
  return min < b && b < max;
}

// Not documented
// Checks that b is in the interval ]min(a, c) , max(a, c)[.
template <class T, class Compare> inline
bool
are_strictly_ordered(const T & a, const T & b, const T & c, Compare cmp)
{
  const T& min = (std::min)(a, c, cmp);
  const T& max = (std::max)(a, c, cmp);
  return cmp(min, b) && cmp(b, max);
}


#ifndef CGAL_NO_DEPRECATED_CODE
template <class ForwardIterator>
inline
CGAL_DEPRECATED
ForwardIterator
successor( ForwardIterator it )
{
  return ++it;
}

template <class BidirectionalIterator>
inline
CGAL_DEPRECATED
BidirectionalIterator
predecessor( BidirectionalIterator it )
{
  return --it;
}
#endif // CGAL_NO_DEPRECATED_CODE

template < class ForwardIterator >
std::pair< ForwardIterator, ForwardIterator >
min_max_element(ForwardIterator first, ForwardIterator last)
{
  typedef std::pair< ForwardIterator, ForwardIterator > FP;
  FP result(first, first);
  if (first != last)
    while (++first != last) {
      if (*first < *result.first)
        result.first = first;
      if (*result.second < *first)
        result.second = first;
    }
  return result;
}

template < class ForwardIterator, class CompareMin, class CompareMax >
std::pair< ForwardIterator, ForwardIterator >
min_max_element(ForwardIterator first,
                ForwardIterator last,
                CompareMin comp_min,
                CompareMax comp_max)
{
  typedef std::pair< ForwardIterator, ForwardIterator > FP;
  FP result(first, first);
  if (first != last)
    while (++first != last) {
      if (comp_min(*first, *result.first))
        result.first = first;
      if (comp_max(*result.second, *first))
        result.second = first;
    }
  return result;
}

template < class ForwardIterator, class Predicate >
ForwardIterator
min_element_if(ForwardIterator first,
               ForwardIterator last,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (*first < *result && pred(*first))
        result = first;
  return result;
}

template < class ForwardIterator, class Compare, class Predicate >
ForwardIterator
min_element_if(ForwardIterator first,
               ForwardIterator last,
               Compare comp,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (comp(*first, *result) && pred(*first))
        result = first;
  return result;
}

template < class ForwardIterator, class Predicate >
ForwardIterator
max_element_if(ForwardIterator first,
               ForwardIterator last,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (*result < *first && pred(*first))
        result = first;
  return result;
}

template < class ForwardIterator, class Compare, class Predicate >
ForwardIterator
max_element_if(ForwardIterator first,
               ForwardIterator last,
               Compare comp,
               Predicate pred)
{
  ForwardIterator result = first = std::find_if(first, last, pred);
  if (first != last)
    while (++first != last)
      if (comp(*result, *first) && pred(*first))
        result = first;
  return result;
}



/*! \brief lexicographic comparison of the two ranges using the \a cmp
    function object.

    Compares the two ranges \c [first1,last1) and \c [first2,last2)
    lexicographically and returns one of the \c CGAL::Comparison_result enum
    values respectively:
      - \c CGAL::SMALLER
      - \c CGAL::EQUAL
      - \c CGAL::LARGER

    \pre The \c value_type of \a InputIterator1 must be convertible
    to the \c first_argument_type of \c BinaryFunction.
    The \c value_type of \a InputIterator2 must be convertible
    to the \c second_argument_type of \c BinaryFunction.
    The \c result_type of \c BinaryFunction must be convertible to
    \c CGAL::Comparison_result.
 */

template <class InputIterator1, class InputIterator2, class BinaryFunction>
CGAL::Comparison_result
lexicographical_compare_three_valued( InputIterator1 first1, InputIterator1 last1,
             InputIterator2 first2, InputIterator2 last2,
             BinaryFunction cmp) {
    while ( first1 != last1 && first2 != last2) {
        CGAL::Comparison_result result = cmp( *first1, *first2);
        if ( result != CGAL::EQUAL)
            return result;
        ++first1;
        ++first2;
    }
    if ( first1 != last1)
        return CGAL::LARGER;
    if ( first2 != last2)
        return CGAL::SMALLER;
    return CGAL::EQUAL;
}

/*! \brief output iterator range to a stream, with separators

    The iterator range \c [first,beyond) is written
    to \c os (obeying CGAL I/O modes). Each element is bracketed by
    \c pre and \c post (default: empty string). Adjacent values are
    spearated by \c sep (default: ", ", i.e. comma space).
    The stream \c os is returned in its new state after output.

    Example:
<PRE>
    int a[] = {1, 2, 3};
    output_range(std::cout, a, a+3, ":", "(", ")");
</PRE>
    produces \c (1):(2):(3)
 */
template <class InputIterator>
std::ostream&
output_range(std::ostream& os,
             InputIterator first, InputIterator beyond,
             const char* sep = ", ", const char* pre = "", const char* post = "")
{
    InputIterator it = first;
    if (it != beyond) {
        os << pre << oformat(*it) << post;
        while (++it != beyond) os << sep << pre << oformat(*it) << post;
    }
    return os;
}


namespace cpp98 {

// Reimplementation of std::random_shuffle, for the use of Spatial_sorting.
// We want an implementation of random_shuffle that produces the same
// result on all platforms, for a given seeded random generator.
template <class RandomAccessIterator,
          class RandomGenerator>
void
random_shuffle(RandomAccessIterator begin, RandomAccessIterator end,
               RandomGenerator& random)
{
  if(begin == end) return;
  for(RandomAccessIterator it = begin + 1; it != end; ++it)
  {
    std::iter_swap( it, begin + random( (it - begin) + 1 ) );
    // The +1 inside random is because random(N) gives numbers in the open
    // interval [0, N[
  }
}


template <class RandomAccessIterator>
void
random_shuffle(RandomAccessIterator begin, RandomAccessIterator end)
{
  typedef std::iterator_traits<RandomAccessIterator> Iterator_traits;
  typedef typename Iterator_traits::difference_type Diff_t;
  boost::rand48 random;
  boost::random_number_generator<boost::rand48, Diff_t> rng(random);
  CGAL::cpp98::random_shuffle(begin,end, rng);
}

} // namespace cpp98

namespace internal {
namespace algorithm {

// Implementation of the algorithm described here:
//   https://en.wikipedia.org/w/index.php?title=Selection_algorithm&oldid=480099620#Partition-based_general_selection_algorithm
template <class RandomAccessIterator, class Compare>
RandomAccessIterator
partition(RandomAccessIterator left,
          RandomAccessIterator right,  // points to the last element of the sequence
          RandomAccessIterator pivot_it,
          Compare& compare)
{
  std::iter_swap(pivot_it, right); // move pivot to the right
  RandomAccessIterator result = left;
  for(RandomAccessIterator it = left; it != right; ++it) {
    if(compare(*it, *right)) {
      std::iter_swap(result, it);
      ++result;
    }
  }
  std::iter_swap(right, result);
  return result;
}

} // end namespace algorithm
} // end namespace internal

// Reimplementation of std::nth_element, for the use of Spatial_sorting.
// We want an implementation of nth_element that produces the same result
// on all platforms.
template <class RandomAccessIterator, class Compare>
void nth_element(RandomAccessIterator left,
                 RandomAccessIterator nth,
                 RandomAccessIterator right,
                 Compare& comp)
{
  if(left == right) return;
  --right; // 'right' points to the last element of the sequence
  if(left == right) return; // exit if there is only one element
  while(true) {
    RandomAccessIterator pivot_it = left + ((right - left) / 2);
    RandomAccessIterator new_pivot_it =
      internal::algorithm::partition(left, right, pivot_it, comp);
    if(new_pivot_it == nth) return;
    if(nth < new_pivot_it)
      right = new_pivot_it - 1;
    else
      left = new_pivot_it + 1;
  } // end while(true)
}

// Not a standard function, but close enough
template <class InputIterator, class Size, class OutputIterator, class UnaryFun>
OutputIterator transform_n( InputIterator first, Size n, OutputIterator result, UnaryFun fun)
{
  // copies the first `n' transformed items from `first' to `result'.
  // Returns the value of `result' after inserting the `n' items.
  while( n--) {
    *result = fun(*first);
    first++;
    result++;
  }
  return result;
}
} //namespace CGAL

#endif // CGAL_ALGORITHM_H
