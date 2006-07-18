// Author(s)     : Menelaos Karavelas

#ifndef CGAL_MY_FILTERED_PREDICATE_H
#define CGAL_MY_FILTERED_PREDICATE_H

#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>

CGAL_BEGIN_NAMESPACE


template <class EP, class AP, class C2E, class C2A, bool Protection = true>
class My_filtered_predicate
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


  template <class A1>
  result_type
  operator()(const A1 &a1) const
  {
    Protect_FPU_rounding<Protection> P;
    Ares res = ap(c2a(a1));
    if ( !is_indeterminate(res) ) { return res; }

    Protect_FPU_rounding<!Protection> P1(CGAL_FE_TONEAREST);
    return ep(c2e(a1));
  }


  template <class A1, class A2>
  result_type
  operator()(const A1 &a1, const A2 &a2) const
  {
    Protect_FPU_rounding<Protection> P;
    Ares res = ap(c2a(a1), c2a(a2));
    if ( !is_indeterminate(res) ) { return res; }

    Protect_FPU_rounding<!Protection> P1(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2));
  }


  template <class A1, class A2, class A3>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
  {
    Protect_FPU_rounding<Protection> P;
    Ares res = ap(c2a(a1), c2a(a2), c2a(a3));
    if ( !is_indeterminate(res) ) { return res; }

    Protect_FPU_rounding<!Protection> P1(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3));
  }


  template <class A1, class A2, class A3, class A4>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4) const
  {
    Protect_FPU_rounding<Protection> P;
    Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4));
    if ( !is_indeterminate(res) ) { return res; }

    Protect_FPU_rounding<!Protection> P1(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4));
  }

  template <class A1, class A2, class A3, class A4, class A5>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5) const
  {
    Protect_FPU_rounding<Protection> P;
    Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5));
    if ( !is_indeterminate(res) ) { return res; }

    Protect_FPU_rounding<!Protection> P1(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5));
  }


  template <class A1, class A2, class A3, class A4, class A5, class A6>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
	     const A5 &a5, const A6 &a6) const
  {
    Protect_FPU_rounding<Protection> P;
    Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4), c2a(a5), c2a(a6));
    if ( !is_indeterminate(res) ) { return res; }

    Protect_FPU_rounding<!Protection> P1(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6));
  }

  //==============================================================
  //==============================================================
  //==============================================================
  //==============================================================
  // MK:: MODIFIED UP TP HERE
  //==============================================================
  //==============================================================
  //==============================================================
  //==============================================================
#if 0
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
      CGAL_PROFILER(std::string("calls to : ") + std::string(CGAL_PRETTY_FUNCTION));
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4),
                    c2a(a5), c2a(a6), c2a(a7));
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("failures of : ") + std::string(CGAL_PRETTY_FUNCTION));
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4), c2e(a5), c2e(a6), c2e(a7));
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      CGAL_PROFILER(std::string("calls to : ") + std::string(__PRETTY_FUNCTION__));
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4),
                    c2a(a5), c2a(a6), c2a(a7), c2a(a8) );
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("failures of : ") + std::string(__PRETTY_FUNCTION__));
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4),
              c2e(a5), c2e(a6), c2e(a7), c2e(a8) );
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      CGAL_PROFILER(std::string("calls to : ") + std::string(__PRETTY_FUNCTION__));
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4),
                    c2a(a5), c2a(a6), c2a(a7), c2a(a8), c2a(a9) );
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("failures of : ") + std::string(__PRETTY_FUNCTION__));
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4),
              c2e(a5), c2e(a6), c2e(a7), c2e(a8), c2e(a9) );
  }
#endif

  template <class A1, class A2, class A3, class A4, class A5, class A6,
            class A7, class A8, class A9, class A10>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3, const A4 &a4,
             const A5 &a5, const A6 &a6, const A7 &a7, const A8 &a8,
             const A9 &a9, const A10 &a10) const
#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  ;
#else
  {
    try
    {
      CGAL_PROFILER(std::string("calls to : ") + std::string(__PRETTY_FUNCTION__));
      Protect_FPU_rounding<Protection> P;
      Ares res = ap(c2a(a1), c2a(a2), c2a(a3), c2a(a4),
                    c2a(a5), c2a(a6), c2a(a7), c2a(a8), c2a(a9), c2a(a10) );
      if (! is_indeterminate(res))
        return res;
    }
    catch (Interval_nt_advanced::unsafe_comparison) {}
    CGAL_PROFILER(std::string("failures of : ") + std::string(__PRETTY_FUNCTION__));
    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    return ep(c2e(a1), c2e(a2), c2e(a3), c2e(a4),
              c2e(a5), c2e(a6), c2e(a7), c2e(a8), c2e(a9), c2e(a10) );
  }
#endif
  // Idem for more than 10 arguments.  Do it on demand.
#endif // #if 0
};


CGAL_END_NAMESPACE

#endif // CGAL_MY_FILTERED_PREDICATE_H
