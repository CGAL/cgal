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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@mpi-sb.mpg.de>

#ifndef CGAL_ALGORITHM_H
#define CGAL_ALGORITHM_H 1

#include <CGAL/basic.h>
#include <CGAL/copy_n.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

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


CGAL_END_NAMESPACE

#endif // CGAL_ALGORITHM_H //
// EOF //
