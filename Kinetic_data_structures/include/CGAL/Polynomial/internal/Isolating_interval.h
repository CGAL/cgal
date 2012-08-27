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

#ifndef CGAL_POLYNOMIAL_INTERNAL_ISOLATING_INTERVAL_H
#define CGAL_POLYNOMIAL_INTERNAL_ISOLATING_INTERVAL_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <iostream>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

//! necessary to support filtered_numbers which I don't like
template <class NT, class Functor>
typename Functor::result_type apply(const Functor &f,  const NT &a)
{
  return f(a);
}


//! necessary to support filtered_numbers which I don't like
template <class NT, class Functor>
typename Functor::result_type apply(const Functor &f,  const NT &a,
				    const NT &b)
{
  return f( a, b);
}


//! necessary to support filtered_numbers which I don't like
template <class NT, class Functor, class Data>
typename Functor::result_type apply(const Functor &f,  const NT &a,
				    const NT &b, const Data &da, const Data &db)
{
  return f( a, b, da, db);
}


//! Define the interface for a bounding interval By convention the intervals are closed.
/*!
  Just remember the bounds.

*/
template <class FT_t>
class Isolating_interval
{
  typedef Isolating_interval<FT_t> This;
public:
  typedef FT_t NT;
  Isolating_interval(): b_(1, -1){}
  Isolating_interval(const NT &lbi): b_(lbi, lbi) {
  }
  Isolating_interval(const NT &lbi, const NT &ubi): b_(lbi, ubi) {
    if(0) print();                        // make sure it is instantiated
    if (lb() > ub()) {
      std::cerr << "b_.first=" << lb() << " and b_.second=" << ub() << std::endl;
    }
    CGAL_Polynomial_assertion(lb() <= ub());
  }

  typedef enum Endpoint {UPPER, LOWER}
    Endpoint;

  //! Damn this is awful, but I can't think of another way for Filtered interval
  template <class Function>
  typename Function::result_type apply_to_endpoint(const Function &f, Endpoint b) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    if (b==UPPER) {
      return apply(f,ub());
    }
    else {
      return apply(f, lb());
    }
  }

  //! Damn this is awful, but I can't think of another way for Filtered interval
  template <class Function>
  typename Function::result_type apply_to_midpoint(const Function &f) {
    CGAL_Polynomial_assertion(lb() <= ub());
    return apply(f,middle());
  }

  //!
  template <class Function>
  typename Function::result_type apply_to_interval(const Function &f) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return apply(f, lb(), ub());
  }

  //! Allow some extra data. This is kind of an evil hack.
  template <class Function, class Data>
  typename Function::result_type apply_to_interval(const Function &f,
						   const Data &d0, const Data &d1) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return apply(f, lb(), ub(), d0, d1);
  }

  This upper_endpoint_interval() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This(ub());
  }

  This lower_endpoint_interval() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This(lb());
  }

  //! Make the interval singular on bound b.
  This endpoint_interval(Endpoint b) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    if (b==UPPER) {
      return This(ub());
    }
    else {
      return This(lb());
    }
  }

  This midpoint_interval() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This(middle());
  }

  //! Split in half
  This first_half() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This(lb(), middle());
  }
  //! Split in half
  This second_half() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This(middle(), ub());
  }

  //! Return both the first and second half; the midpoint is computed once
  std::pair<This,This> split() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    NT mid = middle();
    return std::pair<This,This>(This(lb(), mid), This(mid, ub()));
  }

  //! Returns the interval [(a+(a+b)/2)/2,((a+b)/2+b)/2]
  This middle_half() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    NT mid = middle();
    NT mid_left = midpoint(lb(), middle());
    NT mid_right = midpoint(middle(), ub());

    return This(mid_left, mid_right);
  }

  std::pair<This,This> split_at(const NT& x) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_precondition( !is_singular() &&
		       lower_bound() < x && x < upper_bound() );
    return std::pair<This,This>( This(lower_bound(), x),
				 This(x, upper_bound()) );
  }

  This split_at_first_half(const NT& x) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_precondition( !is_singular() &&
		       lower_bound() < x && x < upper_bound() );
    return This(lower_bound(), x);
  }
  This split_at_second_half(const NT& x) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_precondition( !is_singular() &&
		       lower_bound() < x && x < upper_bound() );
    return This(x, upper_bound());
  }

  This interval_around_endpoint(Endpoint b) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    NT mid = middle();

    if ( b == LOWER ) {
      return This(lower_bound()-NT(1), mid);
    }                                     // end-if

    // b == UPPER
    return This(mid, upper_bound()+NT(1));
  }

  //! See if the interval only contains a point.
  bool is_singular() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return lb()==ub();
  }

  const NT &to_nt() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_precondition(is_singular());
    return lb();
  }

  /*bool overlaps(const This &o) const {
    if (b_.first > o.b_.first && b_.first < o.b_.second) return true;
    else if (b_.second > o.b_.first && b_.first < o.b_.second) return true;
    else if (o.b_.second > b_.first && b_.first < b_.second) return true;
    else if (o.b_.first > b_.first && b_.first < b_.second) return true;
    else if (*this==o) return true;
    else return false;
    }*/

  //! Compare this interval to another.
  Order order(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    if (*this < o) return STRICTLY_BELOW;
    else if (o < *this) return STRICTLY_ABOVE;
    else if (lb() <= o.lb() && ub() >= o.ub()) return CONTAINS;
    else if (o.lb() <= lb() && o.ub() >= ub()) return CONTAINED;
    else if (*this == o) return EQUAL;
    else if (lb() < o.lb()) return BELOW;
    else return ABOVE;
  }

  //! Cut this at bound b on interval o
  /*!  i.e. If bd is UPPER, then o is assumed to be below this (but
    overlaping) and the new interval goes from o.upper_bound() to
    this->upper_bound()
  */
  /*This chop(const This &o, Endpoint bd) const {
    if (bd== UPPER) {
    Polynomial_assertion(lb()== o.lb() && ub() > o.ub());
    return This(o.ub(), ub());
    } else {
    Polynomial_assertion(ub()== o.ub() && lb() < o.lb());
    return This(lb(),o.lb());
    }
    }*/

  int number_overlap_intervals(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    Order ord= order(o);
    switch(ord) {
    case CONTAINS:
      return 3;
    case CONTAINED:
      return 1;
    case BELOW:
    case ABOVE:
      return 2;
    case STRICTLY_ABOVE:
    case STRICTLY_BELOW:
      return 1;
    default:
      return -1;
    }
  }

  // If we have two partially overllaping intervals, then
  // *_overlap_interval  return the three possible intervals of *this:
  // the common interval and the non-common remainders of the two
  // intervals.
  This first_overlap_interval(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    Order ord= order(o);
    switch(ord) {
    case CONTAINS:
      return This(lb(), o.lb());
    case CONTAINED:
    case STRICTLY_ABOVE:
    case STRICTLY_BELOW:
      return *this;
    case BELOW:
      return This(lb(), o.lb());
    case ABOVE:
      return This(lb(), o.ub());
    default:
      return This();
    };
  }

  This second_overlap_interval(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    Order ord= order(o);
    switch(ord) {
    case CONTAINS:
      return This(o.lb(), o.ub());
    case BELOW:
      return This(o.lb(), ub());
    case ABOVE:
      return This(o.ub(), ub());
    default:
      return This();
    };
  }

  This third_overlap_interval(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    Order ord= order(o);
    switch(ord) {
    case CONTAINS:
      return This(o.ub(), ub());
    default:
      return This();
    };
  }

  /*std::pair<This, This> split_with(const This &o) const {
    bool I_would_like_to_get_rid_of_this;
    CGAL_precondition(0);
    Order ord= order(o);
    typedef std::pair<This, This> IP;
    switch(ord){
    case STRICTLY_BELOW:
    return IP(*this, endpoint_interval(UPPER));
    case STRICTLY_ABOVE:
    return IP(endpoint_interval(LOWER, *this));
    case BELOW:
    return IP(This(lb(), o.lb()), This(ub(), o.ub()));
    case CONTAINS:
    return IP(This(lb(), o.lb()), This();
    }
    }*/

  //! Split this and o into 3 parts
  /*!
    this and o must overlap and this must be below or containing o.
    Then the three regions created by the 4 endpoints are returned.
  */
  /*
    void subintervals(const This &o, This &a, This &b, This &c) const {
    bool I_would_like_to_get_rid_of_this;

    Order ord= order(o);
    switch (ord){
    case BELOW:
    a=This(lb(), o.lb());
    b=This(o.lb(), ub());
    c=This(ub(), o.ub());
    break;
    case CONTAINS:
    a=This(lb, o.lb());
    b=o;
    c=This(o.ub(), ub());
    break;
    case CONTAINED:
    a=This(o.lb, lb());
    b=*this;
    c=This(ub(), o.ub());
    break;
    case ABOVE:
    a=This(o.lb(), lb());
    b=This(lb(), o.ub());
    c=This(o.ub(), ub());
    break;
    case STRICTLY_BELOW:
    a=*this;
    b=This(ub(), o.lb());
    c=o;
    break;
    case STRICTLY_ABOVE:
    a=o;
    b=This(o.ub(), lb());
    c=*this;
    break;
    default:
    CGAL_error();
    }
    CGAL_postcondition(a.ub()==b.lb());
    CGAL_postcondition(b.ub()==c.lb());
    }*/

  /**/
  bool operator<(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    if (ub() < o.lb()) return true;
    return (ub() == o.lb() && (!is_singular() || !o.is_singular()));
    /*
      if (ub() == o.lb() && (!is_singular() || !o.is_singular())) return true;
      else return false;
    */
  }
  bool operator>(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    return o < *this;
  }
  bool operator>=(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    return *this >0 || *this==o;
  }
  bool operator<=(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    return *this < 0 || *this==o;
  }
  bool operator==(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    return lb()==o.lb() && ub()==o.ub();
  }
  bool operator!=(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    CGAL_Polynomial_assertion(o.lb() <= o.ub());
    return lb()!=o.lb() || ub()!=o.ub();
  }

  //! Approximate the width with doubles.
  double approximate_width() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return (std::max)(to_double(ub()) - to_double(lb()), 0.0);
  }

  double approximate_relative_width() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return approximate_width()/(std::max)(to_double(abs(ub())), to_double(abs(lb())));
  }

  bool contains(const This &o) {
    CGAL_Polynomial_assertion(lb() <= ub());
    return lb() <= o.lb() && ub() >= o.ub();
  }

  bool contains(const NT &o) {
    return lb() <= o && ub() >= o;
  }

  template <class OStream>
  void write(OStream &out) const
  {
    if (is_singular()) {
      out << lb();
    }
    else {
      out << "(" << lb() << "..." << ub() << ")";
    }
    if (0) print();
  }
  void print() const ;

  This operator-() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This(-ub(), -lb());
  }

  //! return an interval
  std::pair<double, double> double_interval() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    std::pair<double, double>
      lbi= CGAL_POLYNOMIAL_TO_INTERVAL(lb()),
      ubi= CGAL_POLYNOMIAL_TO_INTERVAL(ub());
    return std::pair<double, double>(lbi.first, ubi.second);
  }

  const NT& lower_bound() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    //    bool not_recommended;
    return lb();
  }
  const NT& upper_bound() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    //    bool not_recommended;
    return ub();
  }

  void set_upper(const NT& u) {
    CGAL_Polynomial_assertion(lb() <= ub());
    b_.second = u;
  }

  void set_lower(const NT& l) {
    CGAL_Polynomial_assertion(lb() <= ub());
    b_.first = l;
  }

  /*std::pair<NT, NT> () const {
    return std::pair<NT,NT>(b_.first, b_.second);
    }*/

  //! find the min interval containing both.
  /*This operator||(const This &o){
    return This((std::min)(lb(), o.lb()), (std::max)(ub(), o.ub()));
    }*/

  //! Union
  This operator||(const This &o) const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return This((std::min)(lb(), o.lb()), (std::max)(ub(), o.ub()));
  }

  NT middle() const
  {
    CGAL_Polynomial_assertion(lb() <= ub());
    return NT(0.5)*(ub()+lb());
  }

  const std::pair<NT,NT>& to_pair() const {
    return b_;
  }

protected:
  std::pair<NT,NT> b_;

  NT midpoint(const NT& a, const NT& b) const
  {
    return NT(0.5)*(a+b);
  }

  const NT &lb() const
  {
    return b_.first;
  }
  const NT &ub() const
  {
    return b_.second;
  }

  /*static NT infinity() {
    if (std::numeric_limits<NT>::has_infinity){
    return std::numeric_limits<NT>::infinity();
    } else if (std::numeric_limits<NT>::is_bounded){
    return (std::numeric_limits<NT>::max)();
    } else {
    return NT(1000000000);
    }
    }*/
};

template <class NT>
void Isolating_interval<NT>::print() const
{
  write(std::cout);
  std::cout << std::endl;
}


template <class OStream, class NT>
OStream &operator<<(OStream &out, const Isolating_interval<NT> &ii)
{
  ii.write(out);
  return out;
}


/*template <class NT>
  std::pair<double, double> to_interval(const Isolating_interval<NT> &ii){
  return ii.to_interval();
  }*/

} } } //namespace CGAL::POLYNOMIAL::internal

#ifdef CGAL_POLYNOMIAL_USE_CGAL

namespace CGAL {
template <class NT>
std::pair<double, double> to_interval(const typename CGAL_POLYNOMIAL_NS::internal::Isolating_interval<NT> &ii)
{
  return ii.double_interval();
}


} //namespace CGAL
#endif
#endif
