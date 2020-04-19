// Copyright (c) 2001-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_FILTERED_PREDICATE_H
#define CGAL_FILTERED_PREDICATE_H

#include <string>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>

namespace CGAL {

// This template class is a wrapper that implements the filtering for any
// predicate (dynamic filters with IA).

// TODO :
// - each predicate in the default kernel should define a tag that says if it
//   wants to be filtered or not (=> all homogeneous predicate define this
//   tag).  We could even test-suite that automatically.  It makes a strong
//   new requirement on the kernel though...
//   Could be done with a traits mechanism ?
//   A default template could use the current IA, but other tags or whatever
//   could specify no filtering at all, or static filtering...
// - same thing for constructions => virtual operator() ?
// - similarly, constructions should have a tag saying if they can throw or
//   not, or we let all this up to the compiler optimizer to figure out ?
// - Some caching could be done at the Point_2 level.


template <class EP, class AP, class C2E, class C2A, bool Protection = true>
class Filtered_predicate
{
  C2E c2e;
  C2A c2a;
  EP  ep;
  AP  ap;

  typedef typename AP::result_type  Ares;

public:

  typedef AP    Approximate_predicate;
  typedef EP    Exact_predicate;
  typedef C2E   To_exact_converter;
  typedef C2A   To_approximate_converter;

  typedef typename EP::result_type  result_type;
  // AP::result_type must be convertible to EP::result_type.

  Filtered_predicate()
  {}

  // These constructors are used for constructive predicates.
  // You should try to avoid constructive predicates, as they will construct
  // the exact values systematically (in the ctor), rather than lazily.
  template <class O>
  Filtered_predicate(const O &o1)
    : c2e(), c2a(), ep(c2e(o1)), ap(c2a(o1))
  {}

  template <class O1, class O2>
  Filtered_predicate(const O1 &o1, const O2 &o2)
    : c2e(), c2a(), ep(c2e(o1), c2e(o2)), ap(c2a(o1), c2a(o2))
  {}

  explicit Filtered_predicate(const EP&  e, const AP&  a)
    : ep(e), ap(a)
  {}

  template <typename... Args>
  result_type
  operator()(const Args&... args) const;
};

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <typename... Args>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const Args&... args) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    {
      Protect_FPU_rounding<Protection> p;
      try
        {
          Ares res = ap(c2a(args)...);
          if (is_certain(res))
            return get_certain(res);
        }
      catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(args)...);
}

} //namespace CGAL

#endif // CGAL_FILTERED_PREDICATE_H
