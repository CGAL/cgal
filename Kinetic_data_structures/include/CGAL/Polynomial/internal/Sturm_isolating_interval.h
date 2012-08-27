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

#ifndef CGAL_POLYNOMIAL_STURM_ISOLATING_INTERVAL_H
#define CGAL_POLYNOMIAL_STURM_ISOLATING_INTERVAL_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Filtered_number.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! Define the interface for a bounding interval By convention the intervals are closed.
/*!
  Just remember the bounds.
*/

template <class FT_t>
class Sturm_isolating_interval
{
  typedef Sturm_isolating_interval<FT_t> This;
public:
  typedef Filtered_number<FT_t> NT;
  Sturm_isolating_interval():b_(1, -1){}
  Sturm_isolating_interval(const NT &lbi)
    : b_(lbi, lbi) {}
  Sturm_isolating_interval(const NT &lbi, const NT &ubi)
    : b_(lbi, ubi) {
    if(0) print();                        // make sure it is instantiated
    if (lb() > ub()) {
      std::cerr << "b_.first=" << lb() << " and b_.second=" << ub() << std::endl;
    }
    CGAL_assertion(lb() <= ub());
  }
  bool is_valid() const {
    return b_.first <= b_.second;
  }

  typedef enum Bound {UPPER, LOWER}
    Bound;

  template <class Function>
  typename Function::result_type
  apply_to_bound(Function &f, Bound b) const
  {
    CGAL_precondition(is_valid());
    if (b == UPPER) {
      return apply(f, ub());
    }
    else {
      return apply(f, lb());
    }
  }

  template <class Function>
  typename Function::result_type
  apply_to_interval(const Function &f) const
  {
    CGAL_precondition(is_valid());
    return apply(f, lb(), ub());
  }

  This collapse_to_bound(Bound b) const
  {
    CGAL_precondition(is_valid());
    if (b == UPPER) {
      return This(ub());
    }
    else {
      return This(lb());
    }
  }

  std::pair<This,This> split() const
  {
    CGAL_precondition(is_valid());
    NT mid = middle();
    return std::pair<This,This>(This(lb(), middle()),
				This(middle(), ub()));
  }

  This middle_half() const
  {
    CGAL_precondition(is_valid());
    NT mid = middle();
    NT mid_left = midpoint(lb(), middle());
    NT mid_right = midpoint(middle(), ub());

    return This(mid_left, mid_right);
  }

  This first_half() const
  {
    CGAL_precondition(is_valid());
    return This(lb(), middle());
  }

  This second_half() const
  {
    CGAL_precondition(is_valid());
    return This(middle(), ub());
  }

  bool is_singular() const { return lb() == ub(); }

#if 0
  bool overlaps(const This &o) const
  {
    if (b_.first > o.b_.first && b_.first < o.b_.second) return true;
    else if (b_.second > o.b_.first && b_.first < o.b_.second) return true;
    else if (o.b_.second > b_.first && b_.first < b_.second) return true;
    else if (o.b_.first > b_.first && b_.first < b_.second) return true;
    else if (*this==o) return true;
    else return false;
  }

  typedef enum Relation {LESS, BELOW, CONTAINED, CONTAINS, EQUAL, ABOVE, GREATER}
    Relation;

  Relation relation(const This &o) const
  {
    CGAL_precondition(is_valid());
    if (*this < o) return LESS;
    else if (o < *this) return GREATER;
    else if (lb() <= o.lb() && ub() >= o.ub()) return CONTAINS;
    else if (o.lb() <= lb() && o.ub() >= ub()) return CONTAINED;
    else if (*this == o) return EQUAL;
    else if (lb() < o.lb()) return BELOW;
    else return ABOVE;
  }

  This chop(const This &o, Bound bd) const
  {
    CGAL_precondition(is_valid());
    if (bd== UPPER) {
      CGAL_assertion(lb()== o.lb() && ub() > o.ub());
      return This(o.ub(), ub());
    }
    else {
      CGAL_assertion(ub()== o.ub() && lb() < o.lb());
      return This(lb(),o.lb());
    }
  }

  void write() const
  {
    write(std::cout);
  }

  void split_on(const This &o, This &a, This &b, This &c) const
  {
    CGAL_precondition(is_valid());
    CGAL_precondition_code(Relation rel= relation(o));
    CGAL_precondition(rel== BELOW || rel == CONTAINS);
    a= This(lb(), o.lb());
    if (o.ub() <= ub()) {
      b= This(o);
      c= This(o.ub(), ub());
    }
    else {
      b= This(o.lb(), ub());
      c= This(ub(), o.ub());
    }
  }
#endif

  std::pair<This,This> split_at(const NT& x) const
  {
    CGAL_precondition(is_valid());
    CGAL_precondition( !is_singular() &&
		       lower_bound() < x && x < upper_bound() );
    return std::pair<This,This>( This(lower_bound(), x),
				 This(x, upper_bound()) );
  }

  This interval_around_bound(Bound b) const
  {
    CGAL_precondition(is_valid());
    NT mid = middle();

    if ( b == LOWER ) {
      return This(lower_bound().inf(), mid);
    }                                     // end-if

    // b == UPPER
    return This(mid, upper_bound().sup());
  }

  /**/
  bool operator<(const This &o) const
  {
    if (ub() < o.lb()) return true;
    if (ub() == o.lb() && (!is_singular() || !o.is_singular())) return true;
    else return false;
  }
  bool operator>(const This &o) const
  {
    return o < *this;
  }
  bool operator>=(const This &o) const
  {
    return *this >0 || *this==o;
  }
  bool operator<=(const This &o) const
  {
    return *this < 0 || *this==o;
  }
  bool operator==(const This &o) const
  {
    return lb()==o.lb() && ub()==o.ub();
  }

  double approximate_width() const
  {
    return CGAL::to_double(ub()) - CGAL::to_double(lb());
  }

  NT width() const
  {
    CGAL_precondition(is_valid());
    return upper_bound().nt() - lower_bound().nt();
  }

#if 0
  bool contains(const This &o) {
    return lb() <= o.lb() && ub() >= o.ub();
  }
#endif

  template <class OStream>
  void write(OStream &out) const
  {
    if (is_singular()) {
      out << lb();
    }
    else {
      out << "(" << lb() << "..." << ub() << ")";
    }
  }

  void print() const
  {
    write(std::cout);
  }

  This operator-() const
  {
    return This(-ub(), -lb());
  }

  std::pair<double, double> compute_interval() const
  {
    CGAL_precondition(is_valid());
    std::pair<double, double> lbi =
      CGAL_POLYNOMIAL_TO_INTERVAL(lb()), ubi= CGAL_POLYNOMIAL_TO_INTERVAL(ub());

    return std::pair<double, double>(lbi.first, ubi.second);
  }

#if 0
  //\todo to exact interval somehow

  This operator||(const This &o) {
    return This((std::min)(lb(), o.lb()), (std::max)(ub(), o.ub()));
  }
#endif

  const NT &lower_bound() const
  {
    return b_.first;
  }
  const NT &upper_bound() const
  {
    return b_.second;
  }

  void set_upper(const NT& u) {
    b_.second = u;
  }

  void set_lower(const NT& l) {
    b_.first = l;
  }

  bool is_double() const
  {
    return lower_bound().is_double() && upper_bound().is_double();
  }

  NT middle() const
  {
    return midpoint(ub(),lb());
  }

  const std::pair<NT, NT>& to_pair() const {
    return b_;
  }

protected:
  std::pair<NT, NT> b_;

  const NT &lb() const
  {
    return b_.first;
  }
  const NT &ub() const
  {
    return b_.second;
  }

};

template <class OStream, class NT>
OStream &operator<<(OStream &out, const Sturm_isolating_interval<NT> &ii)
{
  ii.write(out);
  return out;
}


template <class NT>
std::pair<double, double> to_interval(const Sturm_isolating_interval<NT> &ii)
{
  return ii.compute_interval();
}


} } } //namespace CGAL::POLYNOMIAL::internal

namespace CGAL {
template <class NT>
std::pair<double, double>
to_interval(const CGAL_POLYNOMIAL_NS::internal::Sturm_isolating_interval<NT> &ii)
{
  return ii.compute_interval();
}


} //namespace CGAL
#endif                                            // CGAL_POLYNOMIAL_STURM_ISOLATING_INTERVAL_H
