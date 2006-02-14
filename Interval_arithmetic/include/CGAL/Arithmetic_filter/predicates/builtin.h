// Copyright (c) 1999  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H
#define CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H

CGAL_BEGIN_NAMESPACE

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

  static Uncertain<Sign> epsilon_variant( const Restricted_double &a,
                                          const double & epsilon_0)
  {    // return compare_SAF(a,0,epsilon);
    if (to_double(a)> epsilon_0) return POSITIVE;
    if (to_double(a)<-epsilon_0) return NEGATIVE;
    if (to_double(a)==0 && epsilon_0==0) return ZERO;
    return Uncertain<Sign>::indeterminate();
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

  static Uncertain<Comparison_result> epsilon_variant(
          const Restricted_double &a,
          const Restricted_double &b,
          const double & epsilon_0)
  {
    if (to_double(a) > to_double(b)+epsilon_0) return LARGER;
    if (to_double(a) < to_double(b)-epsilon_0) return SMALLER;
    if (to_double(a)==to_double(b) && epsilon_0==0) return EQUAL;
    return Uncertain<Comparison_result>::indeterminate();
  }
};

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H
