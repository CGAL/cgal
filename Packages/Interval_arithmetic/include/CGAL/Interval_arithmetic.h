// Copyright (c) 1998-2005  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERVAL_ARITHMETIC_H
#define CGAL_INTERVAL_ARITHMETIC_H

// This file contains the description of the following classes:
// - Interval_nt<false>  It's a number type that needs the FPU rounding mode
//                       to be set to +inf.  It is also typedef'd to
//                       Interval_nt_advanced for backward compatibility.
// - Interval_nt<true>   Same but it does the rounding mode itself so you
//                       don't have to worry about it.  But it's slower.
//
// Note: When rounding is towards +infinity, to make an operation rounded
// towards -infinity, it's enough to take the opposite of some of the operand,
// and the opposite of the result (see operator+, operator*,...).

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/FPU.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Interval_nt_fwd.h>

CGAL_BEGIN_NAMESPACE

template <bool Protected = true>
class Interval_nt
{
  typedef Interval_nt<Protected>     IA;
  typedef std::pair<double, double>  Pair;

public:

  typedef double      value_type;

  typedef Tag_false   Has_gcd;
  typedef Tag_true    Has_division;
  typedef Tag_true    Has_sqrt;

  // We may have to look back at these...
  typedef Tag_false   Has_exact_ring_operations;
  typedef Tag_false   Has_exact_division;
  typedef Tag_false   Has_exact_sqrt;

  typedef std::exception unsafe_comparison;

  Interval_nt() {}

  Interval_nt(int i)
    : _inf(i), _sup(i) {}

  Interval_nt(double d)
    : _inf(d), _sup(d) {}

  Interval_nt(double i, double s)
    : _inf(i), _sup(s)
  {
      // VC++ should use instead : (i<=s) || !is_valid(i) || !is_valid(s)
      // Or should I use is_valid() ? or is_valid_or_nan() ?
    CGAL_assertion_msg(!(i>s),
	      " Variable used before being initialized (or CGAL bug)");
  }

  // This copy ctor, normally equivalent to the one created by the compiler,
  // appears to fix a code generation problem with GCC 3.0.4...
  // (see test/IA/gcc_3.0.bug.C).
  Interval_nt(const Interval_nt & i)
    : _inf(i._inf), _sup(i._sup) {}

  Interval_nt(const Pair & p)
    : _inf(p.first), _sup(p.second) {}

  static unsigned number_of_failures()
  { return Uncertain<bool>::number_of_failures(); }

  IA operator-() const { return IA (-sup(), -inf()); }

  IA & operator+= (const IA &d) { return *this = *this + d; }
  IA & operator-= (const IA &d) { return *this = *this - d; }
  IA & operator*= (const IA &d) { return *this = *this * d; }
  IA & operator/= (const IA &d) { return *this = *this / d; }
  
  bool is_point() const
  {
    return sup() == inf();
  }

  bool is_same (const IA & d) const
  {
    return inf() == d.inf() && sup() == d.sup();
  }

  bool do_overlap (const IA & d) const
  {
    return !(d.inf() > sup() || d.sup() < inf());
  }

  const double & inf() const { return _inf; }
  const double & sup() const { return _sup; }

  std::pair<double, double> pair() const
  {
    return std::pair<double, double>(_inf, _sup);
  }

  static IA largest()
  {
    return IA(-CGALi::infinity, CGALi::infinity);
  } 

  static IA smallest()
  {
    return IA(-CGAL_IA_MIN_DOUBLE, CGAL_IA_MIN_DOUBLE);
  }

private:
  // Pair inf_sup;
  double _inf, _sup;
};

template <bool Protected>
inline
Uncertain<bool>
operator<(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{
  if (a.sup()  < b.inf()) return true;
  if (a.inf() >= b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{
  if (a.sup() <= b.inf()) return true;
  if (a.inf() >  b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{
  if (b.inf() >  a.sup() || b.sup() <  a.inf()) return false;
  if (b.inf() == a.sup() && b.sup() == a.inf()) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(const Interval_nt<Protected> &a, const Interval_nt<Protected> &b)
{ return ! (a == b); }

// Mixed operators.

template <bool Protected>
inline
Uncertain<bool>
operator<(int a, const Interval_nt<Protected> &b)
{
  if (a < b.inf()) return true;
  if (a >= b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>(int a, const Interval_nt<Protected> &b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(int a, const Interval_nt<Protected> &b)
{
  if (a <= b.inf()) return true;
  if (a >  b.sup()) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(int a, const Interval_nt<Protected> &b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(int a, const Interval_nt<Protected> &b)
{
  if (b.inf() >  a || b.sup() <  a) return false;
  if (b.inf() == a && b.sup() == a) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(int a, const Interval_nt<Protected> &b)
{ return ! (a == b); }

template <bool Protected>
inline
Uncertain<bool>
operator<(const Interval_nt<Protected> &a, int b)
{
  if (a.sup()  < b) return true;
  if (a.inf() >= b) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>(const Interval_nt<Protected> &a, int b)
{ return b < a; }

template <bool Protected>
inline
Uncertain<bool>
operator<=(const Interval_nt<Protected> &a, int b)
{
  if (a.sup() <= b) return true;
  if (a.inf() >  b) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator>=(const Interval_nt<Protected> &a, int b)
{ return b <= a; }

template <bool Protected>
inline
Uncertain<bool>
operator==(const Interval_nt<Protected> &a, int b)
{
  if (b >  a.sup() || b <  a.inf()) return false;
  if (b == a.sup() && b == a.inf()) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
operator!=(const Interval_nt<Protected> &a, int b)
{ return ! (a == b); }


template <bool Protected>
inline
double
to_double (const Interval_nt<Protected> & d)
{
  return (d.sup() + d.inf()) * 0.5;
  // This may overflow...
}

// Returns true if the interval is a unique representable double.
template <bool Protected>
inline
bool
fit_in_double (const Interval_nt<Protected> & d, double &r)
{
  r = d.inf();
  return d.is_point();
}

template <bool Protected>
inline
std::pair<double, double>
to_interval (const Interval_nt<Protected> & d)
{
  return d.pair();
}


template <bool Protected>
inline
bool
is_valid (const Interval_nt<Protected> & d)
{
  return CGAL::is_valid(d.inf()) && 
         CGAL::is_valid(d.sup()) && 
         d.inf() <= d.sup();
}

template <bool Protected>
inline
bool
is_finite (const Interval_nt<Protected> & d)
{
  return CGAL::is_finite(d.inf()) && CGAL::is_finite(d.sup());
}

template <bool Protected>
inline
io_Operator
io_tag (const Interval_nt<Protected> &)
{
  return io_Operator();
}

template <bool Protected>
std::ostream & operator<< (std::ostream &os, const Interval_nt<Protected> & I )
{
  return os << "[" << I.inf() << ";" << I.sup() << "]";
}

template <bool Protected>
std::istream & operator>> (std::istream &is, Interval_nt<Protected> & I)
{
    double d;
    is >> d;
    I = d;
    return is;
}


typedef Interval_nt<false> Interval_nt_advanced;  // for back-compatibility


template <bool Protected>
inline
Interval_nt<Protected>
operator+ (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  Protect_FPU_rounding<Protected> P;
  return Interval_nt<Protected> (-CGAL_IA_SUB(-a.inf(), b.inf()),
                                  CGAL_IA_ADD(a.sup(), b.sup()));
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (double a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)+b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (const Interval_nt<Protected> & a, double b)
{
  return a+Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (int a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)+b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator+ (const Interval_nt<Protected> & a, int b)
{
  return a+Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  Protect_FPU_rounding<Protected> P;
  return Interval_nt<Protected>(-CGAL_IA_SUB(b.sup(), a.inf()),
                                 CGAL_IA_SUB(a.sup(), b.inf()));
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (double a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)-b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (const Interval_nt<Protected> & a, double b)
{
  return a-Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (int a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)-b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator- (const Interval_nt<Protected> & a, int b)
{
  return a-Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  typedef Interval_nt<Protected> IA;
  Protect_FPU_rounding<Protected> P;
  if (a.inf() >= 0.0)					// e>=0
  {
    // b>=0     [a.inf()*b.inf(); a.sup()*b.sup()]
    // b<=0     [a.sup()*b.inf(); a.inf()*b.sup()]
    // b~=0     [a.sup()*b.inf(); a.sup()*b.sup()]
    double aa = a.inf(), bb = a.sup();
    if (b.inf() < 0.0)
    {
	aa = bb;
	if (b.sup() < 0.0)
	    bb = a.inf();
    }
    return IA(-CGAL_IA_MUL(aa, -b.inf()), CGAL_IA_MUL(bb, b.sup()));
  }
  else if (a.sup()<=0.0)				// e<=0
  {
    // b>=0     [a.inf()*b.sup(); a.sup()*b.inf()]
    // b<=0     [a.sup()*b.sup(); a.inf()*b.inf()]
    // b~=0     [a.inf()*b.sup(); a.inf()*b.inf()]
    double aa = a.sup(), bb = a.inf();
    if (b.inf() < 0.0)
    {
	aa=bb;
	if (b.sup() < 0.0)
	    bb=a.sup();
    }
    return IA(-CGAL_IA_MUL(bb, -b.sup()), CGAL_IA_MUL(aa, b.inf()));
  }
  else						// 0 \in [inf();sup()]
  {
    if (b.inf()>=0.0)				// d>=0
      return IA(-CGAL_IA_MUL(-a.inf(), b.sup()),
                 CGAL_IA_MUL(a.sup(), b.sup()));
    if (b.sup()<=0.0)				// d<=0
      return IA(-CGAL_IA_MUL(a.sup(), -b.inf()),
                 CGAL_IA_MUL(a.inf(), b.inf()));
        					// 0 \in d
    double tmp1 = CGAL_IA_MUL(-a.inf(),  b.sup());
    double tmp2 = CGAL_IA_MUL( a.sup(), -b.inf());
    double tmp3 = CGAL_IA_MUL( a.inf(),  b.inf());
    double tmp4 = CGAL_IA_MUL( a.sup(),  b.sup());
    return IA(-std::max(tmp1,tmp2), std::max(tmp3,tmp4));
  }
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (double a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)*b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (const Interval_nt<Protected> & a, double b)
{
  return a*Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (int a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)*b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator* (const Interval_nt<Protected> & a, int b)
{
  return a*Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (const Interval_nt<Protected> &a, const Interval_nt<Protected> & b)
{
  typedef Interval_nt<Protected> IA;
  Protect_FPU_rounding<Protected> P;
  if (b.inf() > 0.0)				// b>0
  {
    // e>=0	[a.inf()/b.sup(); a.sup()/b.inf()]
    // e<=0	[a.inf()/b.inf(); a.sup()/b.sup()]
    // e~=0	[a.inf()/b.inf(); a.sup()/b.inf()]
    double aa = b.sup(), bb = b.inf();
    if (a.inf() < 0.0)
    {
	aa = bb;
	if (a.sup() < 0.0)
	    bb = b.sup();
    };
    return IA(-CGAL_IA_DIV(-a.inf(), aa), CGAL_IA_DIV(a.sup(), bb));
  }
  else if (b.sup()<0.0)			// b<0
  {
    // e>=0	[a.sup()/b.sup(); a.inf()/b.inf()]
    // e<=0	[a.sup()/b.inf(); a.inf()/b.sup()]
    // e~=0	[a.sup()/b.sup(); a.inf()/b.sup()]
    double aa = b.sup(), bb = b.inf();
    if (a.inf() < 0.0)
    {
	bb = aa;
	if (a.sup() < 0.0)
	    aa = b.inf();
    };
    return IA(-CGAL_IA_DIV(-a.sup(), aa), CGAL_IA_DIV(a.inf(), bb));
  }
  else					// b~0
    return IA::largest();
	   // We could do slightly better -> [0;infinity] when b.sup()==0,
	   // but is this worth ?
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (double a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)/b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (const Interval_nt<Protected> & a, double b)
{
  return a/Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (int a, const Interval_nt<Protected> & b)
{
  return Interval_nt<Protected>(a)/b;
}

template <bool Protected>
inline
Interval_nt<Protected>
operator/ (const Interval_nt<Protected> & a, int b)
{
  return a/Interval_nt<Protected>(b);
}

template <bool Protected>
inline
Interval_nt<Protected>
sqrt (const Interval_nt<Protected> & d)
{
  Protect_FPU_rounding<Protected> P;  // not optimal here.
  // sqrt([+a,+b]) => [sqrt(+a);sqrt(+b)]
  // sqrt([-a,+b]) => [0;sqrt(+b)] => assumes roundoff error.
  // sqrt([-a,-b]) => [0;sqrt(-b)] => assumes user bug (unspecified result).
  FPU_set_cw(CGAL_FE_DOWNWARD);
  double i = (d.inf() > 0.0) ? CGAL_IA_SQRT(d.inf()) : 0.0;
  FPU_set_cw(CGAL_FE_UPWARD);
  return Interval_nt<Protected>(i, CGAL_IA_SQRT(d.sup()));
}

template <bool Protected>
inline
Interval_nt<Protected>
min (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  return Interval_nt<Protected>(std::min(d.inf(), e.inf()),
		                std::min(d.sup(), e.sup()));
}

template <bool Protected>
inline
Interval_nt<Protected>
max (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  return Interval_nt<Protected>(std::max(d.inf(), e.inf()),
		                std::max(d.sup(), e.sup()));
}

template <bool Protected>
inline
Interval_nt<Protected>
square (const Interval_nt<Protected> & d)
{
  Protect_FPU_rounding<Protected> P;
  if (d.inf()>=0.0)
      return Interval_nt<Protected>(-CGAL_IA_MUL(d.inf(), -d.inf()),
	                             CGAL_IA_MUL(d.sup(), d.sup()));
  if (d.sup()<=0.0)
      return Interval_nt<Protected>(-CGAL_IA_MUL(d.sup(), -d.sup()),
	     	                     CGAL_IA_MUL(d.inf(), d.inf()));
  return Interval_nt<Protected>(0.0, CGAL_IA_SQUARE(std::max(-d.inf(),
							     d.sup())));
}

template <bool Protected>
inline
Interval_nt<Protected>
abs (const Interval_nt<Protected> & d)
{
  if (d.inf() >= 0.0) return d;
  if (d.sup() <= 0.0) return -d;
  return Interval_nt<Protected>(0.0, std::max(-d.inf(), d.sup()));
}

template <bool Protected>
inline
Uncertain<Sign>
sign (const Interval_nt<Protected> & d)
{
  if (d.inf() > 0.0) return POSITIVE;
  if (d.sup() < 0.0) return NEGATIVE;
  if (d.inf() == d.sup()) return ZERO;
  return Uncertain<Sign>::indeterminate();
}

template <bool Protected>
inline
Uncertain<Comparison_result>
compare (const Interval_nt<Protected> & d, const Interval_nt<Protected> & e)
{
  if (d.inf() > e.sup()) return LARGER;
  if (e.inf() > d.sup()) return SMALLER;
  if (e.inf() == d.sup() && d.inf() == e.sup()) return EQUAL;
  return Uncertain<Comparison_result>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
is_zero (const Interval_nt<Protected> & d)
{
  if (d.inf() > 0.0) return false;
  if (d.sup() < 0.0) return false;
  if (d.inf() == d.sup()) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
is_one (const Interval_nt<Protected> & d)
{
  if (d.inf() > 1) return false;
  if (d.sup() < 1) return false;
  if (d.inf() == d.sup()) return true;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
is_positive (const Interval_nt<Protected> & d)
{
  if (d.inf() > 0.0) return true;
  if (d.sup() <= 0.0) return false;
  return Uncertain<bool>::indeterminate();
}

template <bool Protected>
inline
Uncertain<bool>
is_negative (const Interval_nt<Protected> & d)
{
  if (d.inf() >= 0.0) return false;
  if (d.sup() < 0.0) return true;
  return Uncertain<bool>::indeterminate();
}

inline
std::pair<double, double>
to_interval (const long & l)
{
#ifndef __BORLANDC__ // The stupid Borland compiler generates warnings...
  if (sizeof(double) > sizeof(long)) {
    // On 64bit platforms, a long doesn't fit exactly in a double.
    // Well, a perfect fix would be to use std::numeric_limits<>, but...
    Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
    Interval_nt<false> approx ((double) l);
    FPU_set_cw(CGAL_FE_UPWARD);
    approx += Interval_nt<false>::smallest();
    return approx.pair();
  }
  else
#endif
    return std::pair<double,double>(l,l);
}

CGAL_END_NAMESPACE

#endif // CGAL_INTERVAL_ARITHMETIC_H
