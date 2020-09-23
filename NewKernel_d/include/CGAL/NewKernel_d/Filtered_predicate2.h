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

#ifndef CGAL_FILTERED_PREDICATE2_H
#define CGAL_FILTERED_PREDICATE2_H

#include <string>
#include <CGAL/config.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Profile_counter.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <boost/preprocessor.hpp>

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


template <class K, class EP, class AP, class C2E, class C2A, bool Protection = true>
class Filtered_predicate2
{
//TODO: pack (at least use a tuple)
//FIXME: is it better to store those, or just store enough to recreate them
//(i.e. possibly references to the kernels)?
  CGAL_NO_UNIQUE_ADDRESS EP  ep;
  CGAL_NO_UNIQUE_ADDRESS AP  ap;
  CGAL_NO_UNIQUE_ADDRESS C2E c2e;
  CGAL_NO_UNIQUE_ADDRESS C2A c2a;

public:

  typedef AP    Approximate_predicate;
  typedef EP    Exact_predicate;
  typedef C2E   To_exact_converter;
  typedef C2A   To_approximate_converter;

  // FIXME: should use result_of, see emails by Nico
  typedef typename EP::result_type  result_type;
  // AP::result_type must be convertible to EP::result_type.

  Filtered_predicate2()
  {}

  Filtered_predicate2(const K& k)
    : ep(k.exact_kernel()), ap(k.approximate_kernel()), c2e(k,k.exact_kernel()), c2a(k,k.approximate_kernel())
  {}

  template <typename... Args>
  result_type
  operator()(Args&&... args) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    {
      Protect_FPU_rounding<Protection> p;
      try
        {
          // No forward here, the arguments may still be needed
          auto res = ap(c2a(args)...);
          if (is_certain(res))
            return get_certain(res);
        }
      catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(std::forward<Args>(args))...);
  }
};

} //namespace CGAL

#endif // CGAL_FILTERED_PREDICATE2_H
