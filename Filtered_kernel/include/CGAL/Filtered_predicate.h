// Copyright (c) 2001-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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


#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
  template <typename... Args>
  result_type
  operator()(const Args&... args) const;
#else

  template <class A1>
  result_type
  operator()(const A1 &a1) const;

  template <class A1, class A2>
  result_type
  operator()(const A1 &a1, const A2 &a2) const;

  template <class A1, class A2, class A3>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const;

  template <class A1, class A2, class A3, class A4>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const;

  template <class A1, class A2, class A3, class A4, class A5>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9) const;

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9, class A10>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9, const A10 &a10) const;

  // Idem for more than 10 arguments.  Do it on demand.

#endif
};

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

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
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(args)...);
}

#else

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  
	  Ares res = ap(c2a(a1));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7), c2a(a8));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7),
              c2e(a8));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7), c2a(a8), c2a(a9));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7),
              c2e(a8), c2e(a9));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9, class A10>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9, const A10 &a10) const
{
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    {
      Protect_FPU_rounding<Protection> p;
      try
	{
	  Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6),
			c2a(a7), c2a(a8), c2a(a9), c2a(a10));
	  if (is_certain(res))
	    return get_certain(res);
	}
      catch (Uncertain_conversion_exception) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7),
              c2e(a8), c2e(a9), c2e(a10));
}

#endif

} //namespace CGAL

#endif // CGAL_FILTERED_PREDICATE_H
