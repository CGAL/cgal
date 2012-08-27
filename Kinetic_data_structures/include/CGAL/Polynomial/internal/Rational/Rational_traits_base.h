// Copyright (c) 2005  Stanford University (USA).
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

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Poly>
struct Rational_traits_base
{
  typedef Rational_traits_base<Poly> This;
  typedef Poly Function;
  typedef typename Function::NT FT;

  //! The sign of a polynomial at a rational
  typedef Sign_at_rational<Function> Sign_at;
  Sign_at sign_at_object() const
  {
    return Sign_at();
  }
  
  //! The sign of a polynomial at a rational
  typedef Sign_above_rational<This> Sign_after;
  Sign_after sign_after_object() const
  {
    return Sign_after(*this);
  }

  //! The sign of a polynomial at a rational
  typedef Sign_below_rational<This> Sign_before;
  Sign_before sign_before_object() const
  {
    return Sign_before(*this);
  }

  //! Use strum sequences to compare two roots in an interval
  typedef internal::Compare_isolated_roots_in_interval<This>
  Compare_isolated_roots_in_interval;
  Compare_isolated_roots_in_interval compare_isolated_roots_in_interval_object(const Function &p0,
									       const Function &p1) const
  {
    return Compare_isolated_roots_in_interval(p0, p1, *this);
  }

  // constructions

  //! Compute the quotient
  typedef CGAL::POLYNOMIAL::internal::Quotient<Function> Quotient;
  Quotient quotient_object() const
  {
    return Quotient();
  }

  //! Compute the remainder
  typedef CGAL::POLYNOMIAL::internal::Remainder<Function> Remainder;
  Remainder remainder_object() const
  {
    return Remainder();
  }

  //! Compute the pseudo quotient
  typedef internal::Pseudo_quotient<Function> Pseudo_quotient;
  Pseudo_quotient pseudo_quotient_object() const
  {
    return Pseudo_quotient();
  }

  //! Compute the pseudo remainder
  typedef internal::Pseudo_remainder<Function> Pseudo_remainder;
  Pseudo_remainder pseudo_remainder_object() const
  {
    return Pseudo_remainder();
  }

  //! Return true of two polynomials are negations of one another
  typedef CGAL::POLYNOMIAL::internal::Are_negations<Function> Are_negations;
  Are_negations are_negations_object() const
  {
    return Are_negations();
  }

  //! The the sturm sequence
  typedef CGAL::POLYNOMIAL::internal::Sturm_sequence<This> Sturm_sequence;
  Sturm_sequence Sturm_sequence_object(const Function &f, const Function &g) const
  {
    return Sturm_sequence(f, g, *this);
  }

  //! Compute the derivative
  typedef internal::Derivative<Function> Differentiate;
  Differentiate differentiate_object() const
  {
    return Differentiate();
  }

  //! The the sturm sequence
  typedef CGAL::POLYNOMIAL::internal::Sign_Sturm_sequence<Sturm_sequence> Sign_Sturm_sequence;
  Sign_Sturm_sequence sign_Sturm_sequence_object(const Function &f, const Function &g) const
  {
    return Sign_Sturm_sequence(f, g, *this);
  }

  //! The the standard sequence
  typedef CGAL::POLYNOMIAL::internal::Standard_sequence<Sturm_sequence> Standard_sequence;
  Standard_sequence standard_sequence_object(const Function &f) const
  {
    return Standard_sequence(f, *this);
  }

  //! A bound on the size of roots
  typedef Root_bound_evaluator<Function> Root_bound;
  Root_bound root_bound_object() const
  {
    return Root_bound();
  }

  //! The multiplicity of a rational number
  typedef Rational_multiplicity<This> Multiplicity;
  Multiplicity multiplicity_object() const
  {
    return Multiplicity(*this);
  }

  //! f(x) -> x^d f(1/x)
  /*typedef CGAL::POLYNOMIAL::internal::Invert_variable<Function> Invert_variable;
  Invert_variable invert_variable_object() const
  {
    return Invert_variable();
    }*/

  //! f(x) -> f(-x)
  typedef CGAL::POLYNOMIAL::internal::Negate_variable<Function> Negate_variable;
  Negate_variable negate_variable_object() const
  {
    return Negate_variable();
  }
  /*
  //! Map an interval to positive reals
  typedef CGAL::POLYNOMIAL::internal::Map_rational_interval_to_positive<This> Map_rational_interval_to_positive;
  Map_rational_interval_to_positive map_rational_interval_to_positive_object(const Function &f) const
  {
    return Map_rational_interval_to_positive(f, *this);
  }

  //! Map an interval to positive reals
  typedef CGAL::POLYNOMIAL::internal::Map_rational_interval_to_positive_2<This> Map_rational_interval_to_positive_2;
  Map_rational_interval_to_positive_2 map_rational_interval_to_positive_2_object(const FT &a, const FT &b) const
  {
    return Map_rational_interval_to_positive_2(a,b, *this);
  }

  //! Translates zero by a rational number
  typedef CGAL::POLYNOMIAL::internal::Rational_translate_zero<Function> Rational_translate_zero;
  Rational_translate_zero rational_translate_zero_object(const FT &p) const
  {
    return Rational_translate_zero(p);
  }

  //! multiply by x^m for some m.
  typedef internal::Shift_power<Function> Shift_power;
  Shift_power shift_power_object(unsigned int p) const
  {
    return Shift_power(p);
    }*/

  //! Construct a function
  typedef internal::Construct_function<Function> Construct_function;
  Construct_function construct_function_object() const
  {
    return Construct_function();
  }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
