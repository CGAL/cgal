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

#ifndef CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_BASES_H
#define CGAL_POLYNOMIAL_INTERNAL_VIRTUAL_FUNCTION_BASES_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <iostream>

namespace CGAL { namespace POLYNOMIAL { namespace internal {
#ifndef NDEBUG
int function_handles_created=0;
int function_handles_exact=0;
#endif

template <class Traits_t>
class Filtered_function_node:
  public CGAL::Kinetic::Ref_counted<Filtered_function_node<Traits_t > >
{
  typedef Filtered_function_node<Traits_t> This;
public:
  typedef Traits_t Traits;
  typedef typename Traits::Interval_function Interval_function;
  typedef typename Traits::Exact_function Exact_function;
  typedef typename Traits::Exact_to_interval_converter Interval_function_converter;

  Filtered_function_node(const Interval_function_converter &ifc): 
    has_exact_(false), ifc_(ifc) {
#ifndef NDEBUG
    ++function_handles_created;
#endif
  }
  virtual ~Filtered_function_node(){};

  const Interval_function& interval_function() const
  {
    return fi_;
  }
  const Exact_function& exact_function() const
  {
    if (!has_exact_) {
      generate_exact_function();
    }
    return fe_;
  }

  virtual void write(std::ostream &out) const =0;

  bool has_exact_function() const
  {
    return has_exact_;
  }

  Interval_function_converter interval_function_converter() const
  {
    return ifc_;
  }

  //virtual Filtered_function_node* clone() const =0;
protected:
  void set_interval_function(const Interval_function &fi) const
  {
    fi_= fi;
  }
  void set_exact_function(const Exact_function &fe) const
  {
    fe_= fe;
    has_exact_=true;
#ifndef NDEBUG
    ++function_handles_exact;
#endif
  }

  virtual void generate_exact_function() const =0;

private:
  mutable Interval_function fi_;
  mutable Exact_function fe_;
  mutable bool has_exact_;

  Interval_function_converter ifc_;
};

template <class Traits>
class Filtered_function_node_binary_operation: public Filtered_function_node<Traits>
{
  typedef Filtered_function_node_binary_operation<Traits> This;
  typedef Filtered_function_node<Traits> P;
public:
  Filtered_function_node_binary_operation(const typename P::Handle &lc,
					  const typename P::Handle &rc)
    : P(lc->interval_function_converter()), lc_(lc), rc_(rc){}
  virtual ~Filtered_function_node_binary_operation(){}

  const typename P::Handle& left_child() const
  {
    return lc_;
  }
  const typename P::Handle& right_child() const
  {
    return rc_;
  }
protected:
  void set_exact_function(const typename P::Exact_function &ef) const
  {
    P::set_exact_function(ef);
    // We don't need the children any more
    lc_= NULL;
    rc_= NULL;
    P::set_interval_function(P::interval_function_converter()(P::exact_function()));

  }
  bool has_children() const
  {
    return lc_!= NULL;
  }
  mutable typename P::Handle lc_, rc_;
};

template <class Traits>
class Filtered_function_node_unary_operation: public Filtered_function_node<Traits>
{
  typedef Filtered_function_node_unary_operation<Traits> This;
  typedef Filtered_function_node<Traits> P;
public:
  Filtered_function_node_unary_operation(const typename P::Handle &c):
    P(c->interval_function_converter()), c_(c){}
  virtual ~Filtered_function_node_unary_operation(){}
  const typename P::Handle& child() const
  {
    return c_;
  }
protected:
  bool has_child() const
  {
    return c_!=NULL;
  }
  void set_exact_function(const typename P::Exact_function &ef) const
  {
    // We don't need the children any more
    P::set_exact_function(ef);
    c_= NULL;
    P::set_interval_function(P::interval_function_converter()(P::exact_function()));

  }

  mutable typename P::Handle c_;
};

template <class Traits, class EF,class IF>
class Filtered_function_node_unary_transform: public Filtered_function_node_unary_operation<Traits>
{
  typedef Filtered_function_node_unary_transform<Traits, EF, IF> This;
  typedef Filtered_function_node_unary_operation<Traits> P;
public:
  Filtered_function_node_unary_transform(const typename P::Handle &c,
					 const EF &fe, const IF &fi): P(c), fe_(fe) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(fi(P::child()->interval_function()));
  }
  virtual ~Filtered_function_node_unary_transform(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      fe_.write(out);
      out << "(";
      this->child()->write(out);
      out << ", t)";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(fe_(P::child()->exact_function()));
  }
  EF fe_;
};

template <class Traits, class EF,class IF>
class Filtered_function_node_binary_transform: public Filtered_function_node_binary_operation<Traits>
{
  typedef Filtered_function_node_binary_transform<Traits, EF, IF> This;
  typedef Filtered_function_node_binary_operation<Traits> P;
public:
  Filtered_function_node_binary_transform(const typename P::Handle &lc,
					  const typename P::Handle &rc,
					  const EF &fe, const IF &fi): P(lc, rc), fe_(fe) {
    CGAL_POLYNOMIAL_NS::Interval_arithmetic_guard ig;
    P::set_interval_function(fi(P::left_child()->interval_function(),
				P::right_child()->interval_function()));
  }
  virtual ~Filtered_function_node_binary_transform(){}
  virtual void write(std::ostream &out) const
  {
    if (!P::has_exact_function()) {
      fe_.write(out);
      out << "(";
      this->left_child()->write(out);
      out << ", ";
      this->right_child()->write(out);
      out << ", t)";
    }
    else {
      out << P::exact_function();
    }
  }
protected:
  virtual void generate_exact_function() const
  {
    P::set_exact_function(fe_(P::left_child()->exact_function(),
			      P::right_child()->exact_function()));
  }
  EF fe_;
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
