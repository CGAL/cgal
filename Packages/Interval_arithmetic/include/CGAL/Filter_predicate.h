// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Filter_predicate.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval Arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_FILTER_PREDICATE_H
#define CGAL_FILTER_PREDICATE_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

// This template class is a generic wrapper that implements the filtering
// for all predicates.

// TODO :
// - need to write converters for kernel objects To_Approx() and To_Exact().
//   - CGAL_IA_PROTECTED... this should generalized somehow, if possible.
// - each predicate in the default kernel should define a tag that says if it
//   wants to be filtered or not (=> all homogeneous predicate define this tag).
//   We could even test-suite that automatically.
// - same thing for constructions => virtual operator() ?
// - similarly, constructions should have a tag saying if they can throw or
//   not, or we let all this up to the compiler optimizer to figure out ?
// - Note : using the same principle, we can have a kernel checker easily...
//   [able to generalize and replace PM's checker]
// - Potential caching is done at the Point_2 level.

// Maybe all the ..._object() functions could be implemented in a base class ?
// And easily shared by all kernels ?
// But then we have to duplicate the types ?

// Simple function to check if 2 arguments have the same type,
// ie ~ if 2 types are ~ equal, somehow.
template <class T>
same_type_checker(const T&, const T&)
{}

template <class EP, class AP, class EC, class AC>
class Filtered_predicate
{
  EP Exact_predicate;
  AP Approx_predicate;
  EC To_Exact;
  AC To_Approx;

public:

  typedef typename IP::result_type  result_type;

  Filtered_predicate()
  {
    same_type_checker(IP::result_type(), EP::result_type());
  }

  template <class A1>
  result_type
  operator()(const A1 &a1) const
  {
    try
    {
      Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
      return Approx_predicate(To_Approx(a1));
    }
    catch (Interval_nt_advanced::unsafe_comparison)
    {
      Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
      return Exact_predicate(To_Exact(a1));
    }
  }

  template <class A1, class A2>
  result_type
  operator()(const A1 &a1, const A2 &a2) const
  {
    try
    {
      Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
      return Approx_predicate(To_Approx(a1), To_Approx(a2));
    }
    catch (Interval_nt_advanced::unsafe_comparison)
    {
      Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
      return Exact_predicate(To_Exact(a1), To_Exact(a2));
    }
  }

  template <class A1, class A2, class A3>
  result_type
  operator()(const A1 &a1, const A2 &a2, const A3 &a3) const
  {
    try
    {
      Protect_FPU_rounding<CGAL_IA_PROTECTED> Protection;
      return Approx_predicate(To_Approx(a1), To_Approx(a2), To_Approx(a3));
    }
    catch (Interval_nt_advanced::unsafe_comparison)
    {
      Protect_FPU_rounding<!CGAL_IA_PROTECTED> Protection(CGAL_FE_TONEAREST);
      return Exact_predicate(To_Exact(a1), To_Exact(a2), To_Exact(a3));
    }
  }

  // idem for 4, 5, 6... arguments.  Do it as needed.
};

CGAL_END_NAMESPACE

#endif // CGAL_FILTER_PREDICATE_H
