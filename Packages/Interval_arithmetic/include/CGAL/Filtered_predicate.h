// Copyright (c) 2001-2005  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_FILTER_PREDICATE_H
#define CGAL_FILTER_PREDICATE_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

// This template class is a wrapper that implements the filtering for any
// predicate (dynamic filters with IA).

// TODO :
// - We should specialize it for the static filters versions when possible...
// - each predicate in the default kernel should define a tag that says if it
//   wants to be filtered or not (=> all homogeneous predicate define this
//   tag).  We could even test-suite that automatically.  It makes a strong
//   new requirement on the kernel though...
//   Could be done with a traits mecanism ?
//   A default template could use the current IA, but other tags or whatever
//   could specify no filtering at all, or static filtering...
// - same thing for constructions => virtual operator() ?
// - similarly, constructions should have a tag saying if they can throw or
//   not, or we let all this up to the compiler optimizer to figure out ?
// - Some caching should be done at the Point_2 level.


template <class EP, class AP, class C2E, class C2A, bool Protection = true>
class Filtered_predicate
{
  EP  ep;
  AP  ap;
  C2E c2e;
  C2A c2a;

  typedef typename AP::result_type  Ares;

public:

  typedef AP    Approximate_predicate;
  typedef EP    Exact_predicate;
  typedef C2E   To_exact_converter;
  typedef C2A   To_approximate_converter;

  typedef typename EP::result_type  result_type;
  typedef typename EP::Arity        Arity;
  // AP::result_type must be convertible to EP::result_type.

  Filtered_predicate()
  {}

  // These constructors are used for constructive predicates.
  // You should try to avoid constructive predicates, as they will construct
  // the exact values systematically (in the ctor), rather than lazily.
  template <class O>
  Filtered_predicate(const O &o1)
    : ep(c2e(o1)), ap(c2a(o1))
  {}

  template <class O1, class O2>
  Filtered_predicate(const O1 &o1, const O2 &o2)
    : ep(c2e(o1), c2e(o2)), ap(c2a(o1), c2a(o2))
  {}

  template <class A1>
  result_type
  operator()(const A1 &a1) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1));
  }
#endif

  template <class A1, class A2>
  result_type
  operator()(const A1 &a1, const A2 &a2) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2));
  }
#endif

  template <class A1, class A2, class A3>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3));
  }
#endif

  template <class A1, class A2, class A3, class A4>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6, const A7 &a7) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4),
                    c2a(a5), c2a(a6), c2a(a7));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7));
  }
#endif

  // Idem for more than 7 arguments.  Do it on demand.
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1) const
{
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2) const
{
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
{
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
{
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
{
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5));
}

template <class EP, class AP, class C2E, class C2A, bool Protection>
  template <class A1, class A2, class A3, class A4, class A5, class A6>
typename Filtered_predicate<EP,AP,C2E,C2A,Protection>::result_type
Filtered_predicate<EP,AP,C2E,C2A,Protection>::
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
{
    try
    {
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6));
}
#endif

CGAL_END_NAMESPACE

#endif // CGAL_FILTER_PREDICATE_H
