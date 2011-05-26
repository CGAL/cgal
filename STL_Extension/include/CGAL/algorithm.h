// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion

#ifndef CGAL_ALGORITHM_H
#define CGAL_ALGORITHM_H

#include <CGAL/basic.h>
#include <algorithm>
#include <iosfwd>

namespace CGAL {

namespace cpp0x {
#ifndef CGAL_CFG_NO_CPP0X_COPY_N
  using std::copy_n;
#else
  using CGAL::copy_n;
#endif
} // cpp0x

// copy_n is usually in the STL as well, but not in the official
// standard. We provide our own copy_n.  It is planned for C++0x.

template <class InputIterator, class Size, class OutputIterator>
OutputIterator copy_n( InputIterator first,
                       Size n,
                       OutputIterator result)
{
  // copies the first `n' items from `first' to `result'. Returns
  // the value of `result' after inserting the `n' items.
  while( n--) {
    *result = *first;
    first++;
    result++;
  }
  return result;
}


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


template <class ForwardIterator>
inline
ForwardIterator
successor( ForwardIterator it )
{
  return ++it;
}

template <class BidirectionalIterator>
inline
BidirectionalIterator
predecessor( BidirectionalIterator it )
{
  return --it;
}

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

} //namespace CGAL

#endif // CGAL_ALGORITHM_H
