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

#ifndef CGAL_POLYNOMIAL_STURM_ISOLATING_INTERVAL_H
#define CGAL_POLYNOMIAL_STURM_ISOLATING_INTERVAL_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Filtered_number.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE


//! Define the interface for a bounding interval By convention the intervals are closed.
/*!
  Just remember the bounds.  
*/

template <class FT_t>
class Sturm_isolating_interval {
  typedef Sturm_isolating_interval<FT_t> This;
public:
  typedef Filtered_number<FT_t> NT;
  Sturm_isolating_interval(){}
  Sturm_isolating_interval(const NT &lbi)
    : lb_(lbi), ub_(lbi) {}
  Sturm_isolating_interval(const NT &lbi, const NT &ubi)
    : lb_(lbi), ub_(ubi)
  {
    if(0) print(); // make sure it is instantiated
    if (lb() > ub()){
      std::cerr << "lb_=" << lb() << " and ub_=" << ub() << std::endl;
    }
    CGAL_assertion(lb() <= ub());
  }

  typedef enum Bound {UPPER, LOWER} Bound;

  template <class Function>
  typename Function::result_type
  apply_to_bound(Function &f, Bound b) const
  {
    if (b == UPPER){
      return apply(f, ub());
    } else {
      return apply(f, lb());
    }
  }

  template <class Function>
  typename Function::result_type
  apply_to_interval(const Function &f) const
  {
    return apply(f, lb(), ub());
  }


  This collapse_to_bound(Bound b) const
  {
    if (b == UPPER) {
      return This(ub());
    } else {
      return This(lb());
    }
  }


  std::pair<This,This> split() const
  {
    NT mid = middle();
    return std::pair<This,This>(This(lb(), middle()),
				This(middle(), ub()));
  }

  This middle_half() const
  {
    NT mid = middle();
    NT mid_left = midpoint(lb(), middle());
    NT mid_right = midpoint(middle(), ub());

    return This(mid_left, mid_right);
  }


  This first_half() const {
    return This(lb(), middle());
  }

  This second_half() const {
    return This(middle(), ub());
  }

  bool is_singular() const { return lb() == ub(); }

#if 0
  bool overlaps(const This &o) const
  {
    if (lb_ > o.lb_ && lb_ < o.ub_) return true;
    else if (ub_ > o.lb_ && lb_ < o.ub_) return true;
    else if (o.ub_ > lb_ && lb_ < ub_) return true;
    else if (o.lb_ > lb_ && lb_ < ub_) return true;
    else if (*this==o) return true;
    else return false;
  }

  typedef enum Relation {LESS, BELOW, CONTAINED, CONTAINS, EQUAL, ABOVE, GREATER} Relation;
  
  Relation relation(const This &o) const {
    if (*this < o) return LESS;
    else if (o < *this) return GREATER;
    else if (lb() <= o.lb() && ub() >= o.ub()) return CONTAINS;
    else if (o.lb() <= lb() && o.ub() >= ub()) return CONTAINED;
    else if (*this == o) return EQUAL;
    else if (lb() < o.lb()) return BELOW;
    else return ABOVE;
  }

  This chop(const This &o, Bound bd) const {
    if (bd== UPPER) {
      CGAL_assertion(lb()== o.lb() && ub() > o.ub());
      return This(o.ub(), ub());
    } else {
      CGAL_assertion(ub()== o.ub() && lb() < o.lb());
      return This(lb(),o.lb());
    }
  }
  
  void write() const {
    write(std::cout);
  }

  void split_on(const This &o, This &a, This &b, This &c) const {
    CGAL_precondition_code(Relation rel= relation(o));
    CGAL_precondition(rel== BELOW || rel == CONTAINS);
    a= This(lb(), o.lb());
    if (o.ub() <= ub()){
      b= This(o);
      c= This(o.ub(), ub());
    } else {
      b= This(o.lb(), ub());
      c= This(ub(), o.ub());
    }
  }
#endif

  std::pair<This,This> split_at(const NT& x) const
  {
    CGAL_precondition( !is_singular() &&
		       lower_bound() < x && x < upper_bound() );
    return std::pair<This,This>( This(lower_bound(), x),
				 This(x, upper_bound()) );
  }

  This interval_around_bound(Bound b) const
  {
    NT mid = middle();

    if ( b == LOWER ) {
      return This(lower_bound().inf(), mid);
    } // end-if


    // b == UPPER
    return This(mid, upper_bound().sup());
  }

  /**/
  bool operator<(const This &o) const {
    if (ub() < o.lb()) return true;
    if (ub() == o.lb() && (!is_singular() || !o.is_singular())) return true; 
    else return false;
  }
  bool operator>(const This &o) const {
    return o < *this;
  }
  bool operator>=(const This &o) const {
    return *this >0 || *this==o;
  }
  bool operator<=(const This &o) const {
    return *this < 0 || *this==o;
  }
  bool operator==(const This &o) const {
    return lb()==o.lb() && ub()==o.ub();
  }
  
  double approximate_width() const {
    return CGAL::to_double(ub()) - CGAL::to_double(lb());
  }

  NT width() const {
    return upper_bound().nt() - lower_bound().nt();
  }

#if 0
  bool contains(const This &o){
    return lb() <= o.lb() && ub() >= o.ub();
  }
#endif

  template <class OStream>
  void write(OStream &out) const {
    if (is_singular()){
      out << lb();
    } else {
      out << "(" << lb() << "..." << ub() << ")";
    }
  }

  void print() const {
    write(std::cout);
  }

  This operator-() const {
    return This(-ub(), -lb());
  }

  std::pair<double, double> to_interval() const {
    std::pair<double, double> lbi =
      CGAL::to_interval(lb()), ubi= CGAL::to_interval(ub());

    return std::pair<double, double>(lbi.first, ubi.second);
  }

#if 0
  //\todo to exact interval somehow

  This operator||(const This &o){
    return This(std::min(lb(), o.lb()), std::max(ub(), o.ub()));
  }
#endif

  const NT &lower_bound() const {
    return lb_;
  }
  const NT &upper_bound() const {
    return ub_;
  }

  void set_upper(const NT& u) {
    ub_ = u;
  }

  void set_lower(const NT& l) {
    lb_ = l;
  }

  bool is_double() const {
    return lower_bound().is_double() && upper_bound().is_double();
  }

  NT middle() const {
    return midpoint(ub(),lb());
  }

protected:
  NT lb_, ub_;

  const NT &lb() const {
    return lb_;
  }
  const NT &ub() const {
    return ub_;
  }

 
};

template <class OStream, class NT>
OStream &operator<<(OStream &out, const Sturm_isolating_interval<NT> &ii){
  ii.write(out);
  return out;
}

template <class NT>
std::pair<double, double> to_interval(const Sturm_isolating_interval<NT> &ii){
  return ii.to_interval();
}

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

CGAL_BEGIN_NAMESPACE
template <class NT>
std::pair<double, double>
to_interval(const CGAL_POLYNOMIAL_NS::internal::Sturm_isolating_interval<NT> &ii)
{
  return ii.to_interval();
}
CGAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_STURM_ISOLATING_INTERVAL_H
