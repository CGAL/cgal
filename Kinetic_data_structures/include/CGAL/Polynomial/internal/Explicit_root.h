// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
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
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_EXPLICIT_ROOT_H
#define CGAL_POLYNOMIAL_EXPLICIT_ROOT_H
#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
//! a root that is represented explicitly
/*!
  Not, the number type must have std::numeric_limits defined and
  infinity() or max() must be something reasonable. So exact types
  won't work at the moment.

  Coef_nt is the type of the polynomial coefficients, which can be disjoint from
  the root storage type, if desired. Not sure if this is useful as I am not sure
  you can get too far with an integer type for the coefficients, and anything
  else can be used to store the roots.
*/
template <class NT >
class Explicit_root
{
  typedef Explicit_root<NT> This;
public:

  //! Set it to an invalid value
  Explicit_root(): value_(0), is_inf_(true), mult_(0){
#ifndef NDEBUG
    approximation_=compute_double();
#endif
  }

  template <class CNT>
  Explicit_root(const CNT &v, int mult=1): is_inf_(false), mult_(mult) {
    // Protect CORE::Expr from initialization with doubles (that is what this whole class is about anyway).
    if (std::numeric_limits<CNT>::has_infinity && CGAL::abs(v) == std::numeric_limits<CNT>::infinity()) {
      is_inf_=true;
      if (v > 0) value_=1;
      else value_=1;
    }
    else {
      value_= NT(v);
    }
#ifndef NDEBUG
    approximation_=compute_double();
#endif
  }

  //! Should be protected
  double compute_double() const
  {
    if (is_inf_) {
      if (CGAL::sign(value_)==CGAL::POSITIVE) {
	return infinity<double>();
      }
      else {
	return -infinity<double>();
      }
    }
    return CGAL_POLYNOMIAL_TO_DOUBLE(value_);
  }

  bool is_even_multiplicity() const
  {
    return mult_%2==0;
  }

  //! Do not use, should be protected
  std::pair<double, double> compute_interval() const
  {
    if (is_inf_) {
      return std::pair<double, double>(compute_double(), compute_double());
    }
    return CGAL_POLYNOMIAL_TO_INTERVAL(value_);
  }

  typedef NT Representation;
  const Representation &representation() const
  {
    return value_;
  }

  const Representation &to_rational() const
  {
    return value_;
  }

  template <class OS>
  void write(OS &out) const
  {
    if (is_inf_) {
      if (CGAL::sign(value_)==CGAL::POSITIVE) out << "inf";
      else out << "-inf";
    }
    else {
      out << value_;
    }
  }

  const std::pair<NT, NT>& isolating_interval() const
  {
    static std::pair<NT, NT> ret;
    ret= std::make_pair(value_, value_);
    return ret;
  }

  bool operator<(const This &o) const
  {
    return compare(*this, o)==-1;
  }
  bool operator>(const This &o) const
  {
    return compare(*this, o) ==1;
  }
  bool operator<=(const This &o) const
  {
    return compare(*this, o)!= 1;
  }
  bool operator>=(const This &o) const
  {
    return compare(*this, o) != -1;
  }
  bool operator==(const This &o) const
  {
    return compare(*this, o) ==0;
  }
  bool operator!=(const This &o) const
  {
    return compare(*this, o) !=0;
  }
  This operator-() const
  {
    return This(-value_, is_inf_);
  }

  int multiplicity() const
  {
    return mult_;
  }

  /*void write_type() const {
    std::cout << "general" << std::endl;
    }*/

 
  bool is_rational() const
  {
    // not really true
    return true;
  }
  static This infinity_rep() {
    return This(NT(1), true);
  }
protected:

  static int compare(const This &a, const This &b) {
    int ret=0;
    if (a.is_inf_ && b.is_inf_) {
      if (CGAL::sign(a.value_)== CGAL::sign(b.value_)) ret= 0;
      else if (CGAL::sign(a.value_) == CGAL::POSITIVE ) ret= 1;
      else ret= -1;
    } else if (b.is_inf_) {
      CGAL_assertion(!a.is_inf_);
      if (CGAL::sign(b.value_)== CGAL::NEGATIVE) ret= 1;
      else ret= -1;
    } else if (a.is_inf_) {
      CGAL_assertion(!b.is_inf_);
      if (CGAL::sign(a.value_)== CGAL::NEGATIVE) ret= -1;
      else ret= 1;
    }
    else {
      if (a.value_ > b.value_) ret= 1;
      else if (b.value_ == a.value_) ret= 0;
      else ret= -1;
    }
    return ret;
  }

  Explicit_root(const NT &nt, bool isinf): value_(nt), is_inf_(isinf), mult_(1) {
#ifndef NDEBUG
    approximation_=compute_double();
#endif
  }

  NT value_;
  bool is_inf_;
  int mult_;
#ifndef NDEBUG
  double approximation_;
#endif
};

template <class NT>
std::ostream &operator<<(std::ostream &out, const Explicit_root<NT> &r)
{
  r.write(out);
  return out;
}


/*
  template <class NT>
  double to_double(const Explicit_root<NT> &r){
  return r.to_double();
  }

  template <class NT>
  std::pair<double, double>  to_interval(const Explicit_root<NT> &r){
  return r.to_interval();
  }
*/

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

CGAL_BEGIN_NAMESPACE

template <class NT>
double to_double(const CGAL_POLYNOMIAL_NS::internal::Explicit_root<NT> &r)
{
  return r.compute_double();
}


template <class NT>
std::pair<double, double> to_interval(const CGAL_POLYNOMIAL_NS::internal::Explicit_root<NT> &r)
{
  return r.compute_interval();
}


CGAL_END_NAMESPACE

namespace std
{
  template <class Tr>
  class numeric_limits<CGAL_POLYNOMIAL_NS::internal::Explicit_root<Tr> >: public numeric_limits<Tr>
  {
  public:
    typedef numeric_limits<Tr> P;
    typedef CGAL_POLYNOMIAL_NS::internal::Explicit_root<Tr> T;
    static const bool is_specialized = true;
    static T min BOOST_PREVENT_MACRO_SUBSTITUTION () throw() {return T((P::min)());}
    static T max BOOST_PREVENT_MACRO_SUBSTITUTION () throw() {return T((P::max)());}
    static const bool has_infinity=true;
    static T infinity() throw() {return T::infinity_rep();}
  };
};
#endif
