// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>

#ifndef CGAL_ENVELOPE_SET_DIVIDORS_H
#define CGAL_ENVELOPE_SET_DIVIDORS_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/Random.h>

namespace CGAL
{
namespace Envelope_3
{

//! All even-indexed elements are stored in the first sequence and all the
//  odd-indexed are stored in the second sequence.
class Arbitrary_dividor
{
public:

  template <class InputIterator, class OutputIterator>
  void operator()(InputIterator begin, InputIterator end,
                  OutputIterator set1, OutputIterator set2)
  {
    bool set_first = true;
    for(; begin != end; ++begin, set_first = !set_first)
      if (set_first)
        *set1++ = *begin;
      else
        *set2++ = *begin;

  }

};

//! The last element is stored in the second sequence and all the other (n-1)
//  elments are stored in the first sequence.
class Incremental_dividor
{
public:

  template <class InputIterator, class OutputIterator>
  void operator()(InputIterator begin, InputIterator end,
                  OutputIterator set1, OutputIterator set2)
  {
    std::size_t i, n = std::distance(begin, end);
    // all items but the last go to the first set
    for(i=0; i<(n-1); ++i, ++begin)
      *set1++ = *begin;
    // the last item goes to the second set
    *set2++ = *begin;
  }

};

//! The elements are divided exactly in the middle.
class Middle_dividor
{
public:

  template <class InputIterator, class OutputIterator>
  void operator()(InputIterator begin, InputIterator end,
                  OutputIterator set1, OutputIterator set2)
  {
    std::size_t i, n = std::distance(begin, end);

    for(i=0; i < n; ++i, ++begin)
    {
      if (i < n/2)
        *set1++ = *begin;
      else
        *set2++ = *begin;
    }
  }

};

//! The elements are divided randomly
class Random_dividor
{
public:

  template <class InputIterator, class OutputIterator>
  void operator()(InputIterator begin, InputIterator end,
                  OutputIterator set1, OutputIterator set2)
  {
    while (begin != end)
    {

      if (CGAL::get_default_random().get_bool())
        *set1++ = *begin++;
      else
        *set2++ = *begin++;
    }
  }

};


} // namespace Envelope_3
} // namespace CGAL

#endif
