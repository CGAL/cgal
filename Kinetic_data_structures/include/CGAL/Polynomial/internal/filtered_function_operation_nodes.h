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

#ifndef CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_OPS_H
#define CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_OPS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/filtered_function_node_bases.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Traits>
class Filtered_function_node_plus: public Filtered_function_node_binary_operation<Traits>
{
  typedef Filtered_function_node_plus<Traits> This;
  typedef Filtered_function_node_binary_operation<Traits> P;
public:
  Filtered_function_node_plus(const typename P::Handle &lc,
			      const typename P::Handle &rc): P(lc, rc) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(P::left_child()->interval_function() + P::right_child()->interval_function());
  }
  virtual ~Filtered_function_node_plus(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << "(";
      this->lc_->write(out);
      out << " + ";
      this->rc_->write(out);
      out << ")";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(P::left_child()->exact_function() + P::right_child()->exact_function());
  }
};

template <class Traits>
class Filtered_function_node_times: public Filtered_function_node_binary_operation<Traits>
{
  typedef Filtered_function_node_times<Traits> This;
  typedef Filtered_function_node_binary_operation<Traits> P;
public:
  Filtered_function_node_times(const typename P::Handle &lc,
			       const typename P::Handle &rc): P(lc, rc) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(P::left_child()->interval_function() * P::right_child()->interval_function());
  }
  virtual ~Filtered_function_node_times(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << "(";
      this->lc_->write(out);
      out << " * ";
      this->rc_->write(out);
      out << ")";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(P::left_child()->exact_function() * P::right_child()->exact_function());
  }
};

template <class Traits>
class Filtered_function_node_minus: public Filtered_function_node_binary_operation<Traits>
{
  typedef Filtered_function_node_minus<Traits> This;
  typedef Filtered_function_node_binary_operation<Traits> P;
public:
  Filtered_function_node_minus(const typename P::Handle &lc,
			       const typename P::Handle &rc): P(lc, rc) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(P::left_child()->interval_function() - P::right_child()->interval_function());
  }
  virtual ~Filtered_function_node_minus(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << "(";
      this->lc_->write(out);
      out << " - ";
      this->rc_->write(out);
      out << ")";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(P::left_child()->exact_function() - P::right_child()->exact_function());
  }
};

template <class Traits>
class Filtered_function_node_unary_minus: public Filtered_function_node_unary_operation<Traits>
{
  typedef Filtered_function_node_unary_minus<Traits> This;
  typedef Filtered_function_node_unary_operation<Traits> P;
public:
  Filtered_function_node_unary_minus(const typename P::Handle &c): P(c) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(-P::child()->interval_function());
  }
  virtual ~Filtered_function_node_unary_minus(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << "-";
      this->child()->write(out);
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(-P::child()->exact_function());
  }
};

template <class Traits>
class Filtered_function_node_times_constant: public Filtered_function_node_unary_operation<Traits>
{
  typedef Filtered_function_node_times_constant<Traits> This;
  typedef Filtered_function_node_unary_operation<Traits> P;
public:
  Filtered_function_node_times_constant(const typename P::Handle &c, const typename P::Exact_function::NT &cst): P(c), c_(cst) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(P::child()->interval_function()
			     *CGAL_POLYNOMIAL_NS::To_interval<typename P::Exact_function::NT>()(c_));
  }
  virtual ~Filtered_function_node_times_constant(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << c_ << " * ";
      this->child()->write(out);
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(c_*P::child()->exact_function());
  }
  typename P::Exact_function::NT c_;
};

template <class Traits>
class Filtered_function_node_plus_constant: public Filtered_function_node_unary_operation<Traits>
{
  typedef Filtered_function_node_plus_constant<Traits> This;
  typedef Filtered_function_node_unary_operation<Traits> P;
public:
  Filtered_function_node_plus_constant(const typename P::Handle &c, const typename P::Exact_function::NT &cst): P(c), c_(cst) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(typename P::Interval_function::NT(CGAL_POLYNOMIAL_NS::To_interval<typename P::Exact_function::NT>()(c_))
			     + P::child()->interval_function());
  }
  virtual ~Filtered_function_node_plus_constant(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << "(";
      out << c_ << " + ";
      this->child()->write(out);
      out << ")";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(c_+P::child()->exact_function());
  }
  typename P::Exact_function::NT c_;
};

template <class Traits>
class Filtered_function_node_times_double_constant: public Filtered_function_node_unary_operation<Traits>
{
  typedef Filtered_function_node_times_double_constant<Traits> This;
  typedef Filtered_function_node_unary_operation<Traits> P;
public:
  Filtered_function_node_times_double_constant(const typename P::Handle &c, double d): P(c), c_(d) {
    P::set_interval_function(typename P::Interval_function::NT(CGAL_POLYNOMIAL_NS::To_interval<double>()(c_))*P::child()->interval_function());
  }
  virtual ~Filtered_function_node_times_double_constant(){}
  virtual void write(std::ostream &out) const
  {
    if (P::has_exact_function()) {
      out << c_<< " * ";
      this->child()->write(out);
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(typename P::Exact_function::NT(c_)*P::child()->exact_function());
  }
  double c_;
};

template <class Traits>
class Filtered_function_node_plus_double_constant: public Filtered_function_node_unary_operation<Traits>
{
  typedef Filtered_function_node_plus_double_constant<Traits> This;
  typedef Filtered_function_node_unary_operation<Traits> P;
public:
  Filtered_function_node_plus_double_constant(const typename P::Handle &c, double cst): P(c), c_(cst) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function( P::child()->interval_function()+ CGAL_POLYNOMIAL_NS::To_interval<double>()(c_));
  }
  virtual ~Filtered_function_node_plus_double_constant(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      out << "(";
      out << c_<< " + ";
      this->child()->write(out);
      out << ")";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(P::child()->exact_function() + typename P::Exact_function::NT(c_));
  }
  double c_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
