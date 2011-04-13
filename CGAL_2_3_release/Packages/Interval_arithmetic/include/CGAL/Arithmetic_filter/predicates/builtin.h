// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : include/CGAL/Arithmetic_filter/predicates/builtin.h
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// ======================================================================

#ifndef CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H
#define CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H

CGAL_BEGIN_NAMESPACE

namespace NTS { // because the script transforms CGAL_NTS sign() into ...

struct Static_Filtered_sign_1
{
  static double _bound;
  static double _epsilon_0;
  static Sign update_epsilon( const Static_filter_error &a, double & epsilon_0)
  {
      epsilon_0 = a.error();
      return ZERO;
  }

  static void new_bound (const double b) // , const double error = 0)
  {
    _bound = b;
    (void) update_epsilon(b,_epsilon_0);
  }

  static Sign epsilon_variant( const Restricted_double &a,
                                const double & epsilon_0)
  {    // return compare_SAF(a,0,epsilon);
    if (to_double(a)> epsilon_0) return POSITIVE;
    if (to_double(a)<-epsilon_0) return NEGATIVE;
    if (to_double(a)==0 && epsilon_0==0) return ZERO;
    Interval_nt_advanced::number_of_failures++;
    throw Interval_nt_advanced::unsafe_comparison();
  }
};

struct Static_Filtered_compare_2
{
  static double _bound;
  static double _epsilon_0;
  static Comparison_result update_epsilon(
          const Static_filter_error &a,
          const Static_filter_error &b,
          double & epsilon_0)
  {
      Static_filter_error c = a-b;
      epsilon_0 = c.error();
      return EQUAL;
  }

  static void new_bound (const double b) // , const double error = 0)
  {
    _bound = b;
    (void) update_epsilon(b,b,_epsilon_0);
  }

  static Comparison_result epsilon_variant(
          const Restricted_double &a,
          const Restricted_double &b,
          const double & epsilon_0)
  {
    if (to_double(a) > to_double(b)+epsilon_0) return LARGER;
    if (to_double(a) < to_double(b)-epsilon_0) return SMALLER;
    if (to_double(a)==to_double(b) && epsilon_0==0) return EQUAL;
    Interval_nt_advanced::number_of_failures++;
    throw Interval_nt_advanced::unsafe_comparison();
  }
};

struct Static_Filtered_lexicographical_sign_2
{
  static double _bound;
  static double _epsilon_0;
  static Sign update_epsilon(
          const Static_filter_error &,
          const Static_filter_error &,
          double &)
  {
      // Not finished.
      return ZERO;
  }

  static void new_bound (const double b) // , const double error = 0)
  {
    _bound = b;
    (void) update_epsilon(b,b,_epsilon_0);
  }

  static Sign epsilon_variant(
          const Restricted_double &,
          const Restricted_double &,
          const double &)
  {
    // Not finished.
    Interval_nt_advanced::number_of_failures++;
    throw Interval_nt_advanced::unsafe_comparison();
  }
};

} // namespace NTS

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H
