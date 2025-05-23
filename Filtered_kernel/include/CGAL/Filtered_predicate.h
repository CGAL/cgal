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

#include <type_traits>

namespace CGAL {

// This template class is a wrapper that implements the filtering for any
// predicate (dynamic filters with IA).

// TODO :
// - each predicate in the default kernel should define a tag that says if it
//   wants to be filtered or not (=> all homogeneous predicates define this
//   tag).  We could even test-suite that automatically.  It makes a strong
//   new requirement on the kernel though...
//   Could be done with a traits mechanism ?
//   A default template could use the current IA, but other tags or whatever
//   could specify no filtering at all, or static filtering...
// - same thing for constructions => virtual operator() ?
// - similarly, constructions should have a tag saying if they can throw or
//   not, or we let all this up to the compiler optimizer to figure out ?
// - Some caching could be done at the Point_2 level.

// Protection is undocumented and currently always true, meaning that it
// assumes a default rounding mode of round-to-nearest. false would correspond
// to a default of round-towards-infinity, so interval arithmetic does not
// require protection but regular code may.

template <class EP, class AP, class C2E, class C2A, bool Protection = true>
class Filtered_predicate
{
  C2E c2e;
  C2A c2a;
  EP  ep;
  AP  ap;

public:
  // AP's result type must be convertible to EP's result type.
  typedef AP    Approximate_predicate;
  typedef EP    Exact_predicate;
  typedef C2E   To_exact_converter;
  typedef C2A   To_approximate_converter;

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
  auto
  operator()(const Args&... args) const
  {
    typedef typename Remove_needs_FT<CGAL::cpp20::remove_cvref_t<decltype(ep(c2e(args)...))> >::Type result_type;

#ifndef CGAL_EPICK_NO_INTERVALS
    typedef typename Remove_needs_FT<CGAL::cpp20::remove_cvref_t<decltype(ap(c2a(args)...))> >::Type Ares;

    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    {
      Protect_FPU_rounding<Protection> p;
      try
        {
          Ares res = ap(c2a(args)...);
          if (is_certain(res))
            return result_type(get_certain(res));
        }
      catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
#endif // CGAL_EPICK_NO_INTERVALS
    return result_type(ep(c2e(args)...));
  }
};

template <class EP_RT, class EP_FT, class AP, class C2E_RT, class C2E_FT, class C2A, bool Protection = true>
class Filtered_predicate_RT_FT
{
  C2E_RT c2e_rt;
  C2E_FT c2e_ft;
  C2A c2a;
  EP_RT ep_rt;
  EP_FT ep_ft;
  AP ap;

private:
  // Detect if the predicate's result type has been wrapped with the `Needs_FT` class
  template <typename... Args>
  struct Call_operator_needs_FT
  {
    template <typename T>
    struct is_Needs_FT : std::false_type { };
    template <typename ...T>
    struct is_Needs_FT<Needs_FT<T...> > : std::true_type { };

    typedef CGAL::cpp20::remove_cvref_t<decltype(ap(c2a(std::declval<const Args&>())...))> Actual_approx_res;
    enum { value = is_Needs_FT<Actual_approx_res>::value };
  };

  // If there is no `Needs_FT` in the result, then we can use an RT-based exact predicate
  template <typename... Args,
            std::enable_if_t<Call_operator_needs_FT<Args...>::value>* = nullptr>
  decltype(auto) exact_call(const Args&... args) const { return ep_ft(c2e_ft(args)...); }

  template <typename... Args,
            std::enable_if_t<! Call_operator_needs_FT<Args...>::value>* = nullptr>
  decltype(auto) exact_call(const Args&... args) const { return ep_rt(c2e_rt(args)...); }

public:
  // ## Important note
  //
  // If you want to remove of rename that member function template `needs_FT`,
  // please also change the lines with
  // `CGAL_GENERATE_MEMBER_DETECTOR(needs_FT);`
  // or `has_needs_FT<typename R::Compare_distance_3>` in
  // the file `Kernel_23/test/Kernel_23/include/CGAL/_test_new_3.h`.
  template <typename... Args>
  bool needs_FT(const Args&...) const { return Call_operator_needs_FT<Args...>::value; }

  template <typename... Args>
  auto
  operator()(const Args&... args) const
  {
    typedef typename Remove_needs_FT<CGAL::cpp20::remove_cvref_t<decltype(exact_call(args...))> >::Type result_type;

#ifndef CGAL_EPICK_NO_INTERVALS
    typedef typename Remove_needs_FT<CGAL::cpp20::remove_cvref_t<decltype(ap(c2a(args)...))> >::Type Ares;

    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    {
      Protect_FPU_rounding<Protection> p;
      try
        {
          Ares res = ap(c2a(args)...);
          if (is_certain(res))
            return result_type(get_certain(res));
        }
      catch (Uncertain_conversion_exception&) {}
    }
    CGAL_BRANCH_PROFILER_BRANCH(tmp);
    Protect_FPU_rounding<!Protection> p(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
#endif // CGAL_EPICK_NO_INTERVALS
    return result_type(exact_call(args...));
  }
};

} // namespace CGAL

#endif // CGAL_FILTERED_PREDICATE_H
