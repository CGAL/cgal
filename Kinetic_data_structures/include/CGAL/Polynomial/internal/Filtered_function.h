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

#ifndef CGAL_POLYNOMIAL_FILTERED_FUNCTION_H
#define CGAL_POLYNOMIAL_FILTERED_FUNCTION_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/filtered_function_node_bases.h>
#include <CGAL/Polynomial/internal/filtered_function_leaf_nodes.h>
#include <CGAL/Polynomial/internal/filtered_function_operation_nodes.h>
#include <iterator>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! A function which supports filtering.
/*!  This class wraps a tree for generating functions, allowing the
  exact version of the function to be generated in a lazy manner
*/
template <class Traits >
class Filtered_function
{
protected:
  typedef typename Traits::Exact_function EF;
  typedef typename Traits::Interval_function IF;
  typedef typename Traits::Exact_to_interval_converter IFC;
  typedef Filtered_function<Traits> This;
  typedef typename internal::Filtered_function_node<Traits> VF;

  //! \todo should be hidden
  typedef typename VF::Handle VFP;
  typedef VF Node_base;
public:
  typedef Node_base* Raw_pointer;
  typedef IFC Interval_function_converter;
  typedef EF Exact_function;
  typedef IF Interval_function;
  typedef typename Exact_function::NT NT;

  Filtered_function(){}

  Filtered_function(const NT &bt) {
    ptr_= new internal::Filtered_function_node_constant< Traits>(bt, IFC());
  }

  Filtered_function(double d) {
    ptr_= new internal::Filtered_function_node_double_constant<Traits>(d, IFC());
  }

  Filtered_function(int d) {
    ptr_= new internal::Filtered_function_node_double_constant<Traits>(d, IFC());
  }

  template <class It>
  Filtered_function(It b, It e) {
    typename std::iterator_traits<It>::value_type t;
    ptr_= make_node(b, e,  t);
  }

  Filtered_function(const EF &ef) {
    ptr_= new internal::Filtered_function_node_explicit<Traits>(ef, IFC());
  }

  Filtered_function(Raw_pointer vfp): ptr_(vfp){}

  //! Return the function computed using an interval type.
  const Interval_function& interval_function() const
  {
    return ptr_->interval_function();
  }

  //! Return the function computed using an exact type.
  const Exact_function& exact_function() const
  {
    return ptr_->exact_function();
  }

  //! If you don't want to rewrite some code.
  NT operator[](unsigned int i) const
  {
    return exact_function()[i];
  }

  //! Note this is polynomial equality, not equality of roots or anything like that
  /*!
    I think this is optimal.
  */
  bool operator==(const This &o) const
  {
    if (ptr_ == o.ptr_) return true;      // want pointer comparison
	   int md= (std::min)(interval_function().degree(), o.interval_function().degree());
    for (int i=interval_function().degree(); i>= md+1; --i) {
      if (!interval_function()[i].do_overlap(0)) return false;
    }
    for (int i=o.interval_function().degree(); i>= md+1; --i) {
      if (!o.interval_function()[i].do_overlap(0)) return false;
    }
    for (int i=md; i>=0; --i) {
      if (!interval_function()[i].do_overlap(o.interval_function()[i])) return false;
    }
    if (exact_function().degree() != o.exact_function().degree()) return false;
    else for (int i=0; i<= o.exact_function().degree(); ++i) {
      if (exact_function()[i] != o.exact_function()[i]) return false;
    }
    return true;
  }

  //! neq
  bool operator!=(const This &o) const
  {
    return !operator==(o);
  }

  void write(std::ostream &out) const
  {
    if (ptr_ == NULL) out << "NULL";
    else ptr_->write(out);
  }

  void print() const
  {
    write(std::cout);
  }

  //! Compute the degree.
  /*!
    To do so, first look at the interval function and see if we can
    determine the degree from it.  If not, fall back on exact.
  */
  int degree() const
  {
    if (ptr_->has_exact_function()) return exact_function().degree();
    if (interval_function().degree()<0) return interval_function().degree();
    else if (extended_sign(interval_function()[interval_function().degree()]) ==EXTENDED_POSITIVE
	     || extended_sign(interval_function()[interval_function().degree()]) ==EXTENDED_NEGATIVE) {
      return interval_function().degree();
    } else return exact_function().degree();
  }

  //! Check if it is constant
  /*!
    \todo check if the exact is know and use it if it is
  */
  bool is_zero() const
  {
    if (interval_function().is_zero()) return true;
    //else if (!interval_function().is_zero()) return false;
    else return exact_function().is_zero();
  }

  //! Check if it is constant
  /*!
    \todo check if the exact is know and use it if it is
  */
  bool is_constant() const
  {
    if (interval_function().is_constant()) return true;
    //else if (!interval_function().is_constant()) return false;
    else return exact_function().is_constant();
  }

  This operator+(const This &o) const
  {
    return This(new internal::Filtered_function_node_plus<Traits>(ptr_, o.ptr_));
  }
  This operator+(const NT &o) const
  {
    return This(new internal::Filtered_function_node_plus_constant<Traits>(ptr_, o));
  }
  This operator+(double o) const
  {
    return This(new internal::Filtered_function_node_plus_double_constant<Traits>(ptr_, o));
  }
  This operator-() const
  {
    return This(new internal::Filtered_function_node_unary_minus<Traits>(ptr_));
  }
  This operator-(const This &o) const
  {
    return This(new internal::Filtered_function_node_minus<Traits>(ptr_, o.ptr_));
  }
  This operator-(const NT &o) const
  {
    return This(new internal::Filtered_function_node_plus_constant<Traits>(ptr_, -o));
  }
  This operator-(double o) const
  {
    return This(new internal::Filtered_function_node_plus_double_constant<Traits>(ptr_, -o));
  }
  This operator*(const This &o) const
  {
    return This(new internal::Filtered_function_node_times<Traits>(ptr_, o.ptr_));
  }
  This operator*(const NT &o) const
  {
    return This(new internal::Filtered_function_node_times_constant<Traits>(ptr_, o));
  }
  This operator*(double o) const
  {
    return This(new internal::Filtered_function_node_times_double_constant<Traits>(ptr_, o));
  }
  This operator/(const NT &o) const
  {
    return This(new internal::Filtered_function_node_times_constant<Traits>(ptr_, 1/o));
  }
  This operator/(double o) const
  {
    return This(new internal::Filtered_function_node_times_double_constant<Traits>(ptr_, 1/o));
  }

  Raw_pointer tree() const
  {
    return ptr_.get();
  }

protected:
  template <class It, class NTT>
  VF* make_node(It b, It e, const NTT &) {
    return new internal::Filtered_function_node_explicit<Traits>(Exact_function(b,e), IFC());
  }
  template <class It>
  VF* make_node(It b, It e, double) {
    return new internal::Filtered_function_node_explicit<Traits>(Interval_function(b,e), IFC());
  }

  VFP ptr_;
};

template <class Traits>
std::ostream &operator<<(std::ostream &out, const Filtered_function<Traits> &ff)
{
  ff.write(out);
  return out;
}


template <class Traits>
Filtered_function<Traits> operator*(const typename Traits::Exact_function::NT &nt, const Filtered_function<Traits> &ff)
{
  return ff.operator*(nt);
}


template <class Traits>
Filtered_function<Traits> operator*(double nt, const Filtered_function<Traits> &ff)
{
  return ff.operator*(nt);
}


template <class Traits>
Filtered_function<Traits> operator+(double nt, const Filtered_function<Traits> &ff)
{
  return ff.operator+(nt);
}


template <class Traits>
Filtered_function<Traits> operator+(const typename Traits::Exact_function::NT &nt, const Filtered_function<Traits> &ff)
{
  return ff.operator*(nt);
}


} } } //namespace CGAL::POLYNOMIAL::internal
#endif
