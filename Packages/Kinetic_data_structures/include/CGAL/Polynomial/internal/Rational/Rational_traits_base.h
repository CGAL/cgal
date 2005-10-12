// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_RATIONAL_TRAITS_BASE_H
#define CGAL_POLYNOMIAL_INTERNAL_RATIONAL_TRAITS_BASE_H

#include <CGAL/Polynomial/basic.h>


#include <CGAL/Polynomial/internal/Rational/Sign_at_rational.h>
#include <CGAL/Polynomial/internal/Rational/Sign_above_rational.h>
#include <CGAL/Polynomial/internal/Rational/Sign_below_rational.h>
#include <CGAL/Polynomial/internal/Rational/Compare_isolated_roots_in_interval.h>
#include <CGAL/Polynomial/internal/Rational/Derivative.h>
#include <CGAL/Polynomial/internal/Rational/Construct_function.h>
#include <CGAL/Polynomial/internal/Rational/Are_negations.h>
#include <CGAL/Polynomial/internal/Rational/Invert_variable.h>
#include <CGAL/Polynomial/internal/Rational/Negate_variable.h>
#include <CGAL/Polynomial/internal/Rational/Map_rational_interval_to_positive.h>
#include <CGAL/Polynomial/internal/Rational/Rational_translate_zero.h>
#include <CGAL/Polynomial/internal/Rational/Shift_power.h>
#include <CGAL/Polynomial/internal/Rational/Rational_multiplicity.h>
#include <CGAL/Polynomial/internal/Rational/Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Standard_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Sign_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Root_bound_evaluator.h>
#include <CGAL/Polynomial/internal/Rational/Pseudo_quotient.h>
#include <CGAL/Polynomial/internal/Rational/Pseudo_remainder.h>
#include <CGAL/Polynomial/internal/Rational/Quotient.h>
#include <CGAL/Polynomial/internal/Rational/Remainder.h>



CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Poly>
struct Rational_traits_base {
  typedef Rational_traits_base<Poly> This;
  typedef Poly Function;
  typedef typename Function::NT NT;
  
  //! The sign of a polynomial at a rational
  typedef Sign_at_rational<Function> Sign_at;
  Sign_at sign_at_object(const Function &f) const {
    return Sign_at(f);
  }

  //! The sign of a polynomial at a rational
  typedef Sign_above_rational<This> Sign_above;
  Sign_above sign_above_object(const Function &f) const {
    return Sign_above(f);
  }

 //! The sign of a polynomial at a rational
  typedef Sign_below_rational<This> Sign_below;
  Sign_below sign_below_object(const Function &f) const {
    return Sign_below(f);
  }
 
  //! Use strum sequences to compare two roots in an interval
  typedef internal::Compare_isolated_roots_in_interval<This>
  Compare_isolated_roots_in_interval;
  Compare_isolated_roots_in_interval compare_isolated_roots_in_interval_object(const Function &p0, 
									       const Function &p1) const {
    return Compare_isolated_roots_in_interval(p0, p1, *this);
  }


  // constructions


  //! Compute the quotient
  typedef internal::Quotient<Function> Quotient;
  Quotient quotient_object() const {
    return Quotient();
  }

  //! Compute the remainder
  typedef internal::Remainder<Function> Remainder;
  Remainder remainder_object() const {
    return Remainder();
  }

  //! Compute the pseudo quotient
  typedef internal::Pseudo_quotient<Function> Pseudo_quotient;
  Pseudo_quotient pseudo_quotient_object() const {
    return Pseudo_quotient();
  }

  //! Compute the pseudo remainder
  typedef internal::Pseudo_remainder<Function> Pseudo_remainder;
  Pseudo_remainder pseudo_remainder_object() const {
    return Pseudo_remainder();
  }

  //! Return true of two polynomials are negations of one another
  typedef Are_negations<Function> Are_negations;
  Are_negations are_negations_object() const {
    return Are_negations();
  }

  //! The the sturm sequence
  typedef Sturm_sequence<This> Sturm_sequence;
  Sturm_sequence Sturm_sequence_object(const Function &f, const Function &g) const {
    return Sturm_sequence(f, g, *this);
  }

  //! Compute the derivative
  typedef internal::Derivative<Function> Differentiate;
  Differentiate differentiate_object() const {
    return Differentiate();
  }


  //! The the sturm sequence
  typedef Sign_Sturm_sequence<Sturm_sequence> Sign_Sturm_sequence;
  Sign_Sturm_sequence sign_Sturm_sequence_object(const Function &f, const Function &g) const {
    return Sign_Sturm_sequence(f, g, *this);
  }

  //! The the standard sequence
  typedef Standard_sequence<Sturm_sequence> Standard_sequence;
  Standard_sequence standard_sequence_object(const Function &f) const {
    return Standard_sequence(f, *this);
  }

  //! A bound on the size of roots
  typedef Root_bound_evaluator<Function> Root_bound;
  Root_bound root_bound_object() const {
    return Root_bound();
  }
  
  //! The multiplicity of a rational number
  typedef Rational_multiplicity<This> Multiplicity;
  Multiplicity multiplicity_object(const Function &f) const {
    return Multiplicity(f, *this);
  }

  //! f(x) -> x^d f(1/x)
  typedef Invert_variable<Function> Invert_variable;
  Invert_variable invert_variable_object() const {
    return Invert_variable();
  }

  //! f(x) -> f(-x)
  typedef Negate_variable<Function> Negate_variable;
  Negate_variable negate_variable_object() const {
    return Negate_variable();
  }

  //! Map an interval to positive reals
  typedef Map_rational_interval_to_positive<This> Map_rational_interval_to_positive;
  Map_rational_interval_to_positive map_rational_interval_to_positive_object(const Function &f) const {
    return Map_rational_interval_to_positive(f, *this);
  }

  //! Map an interval to positive reals
  typedef Map_rational_interval_to_positive_2<This> Map_rational_interval_to_positive_2;
  Map_rational_interval_to_positive_2 map_rational_interval_to_positive_2_object(const NT &a, const NT &b) const {
    return Map_rational_interval_to_positive_2(a,b, *this);
  }

  //! Translates zero by a rational number
  typedef Rational_translate_zero<Function> Rational_translate_zero;
  Rational_translate_zero rational_translate_zero_object(const NT &p) const {
    return Rational_translate_zero(p);
  }

  //! multiply by x^m for some m. 
  typedef internal::Shift_power<Function> Shift_power;
  Shift_power shift_power_object(unsigned int p) const {
    return Shift_power(p);
  }

  //! Construct a function
  typedef internal::Construct_function<Function> Construct_function;
  Construct_function construct_function_object() const {
    return Construct_function();
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
