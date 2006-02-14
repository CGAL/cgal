// Copyright (c) 2005  Utrecht University (The Netherlands),
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

#ifndef CGAL_NUMBER_TYPE_CHECKER_H
#define CGAL_NUMBER_TYPE_CHECKER_H

#include <CGAL/basic.h>
#include <CGAL/Number_type_checker_fwd.h>

// A number type class, parameterized by 2 number types NT1 and NT2.
// It runs all operations on parallel over NT1 and NT2.
// It is also parameterized by a comparator which compares the values
// of NT1 and NT2 after all arithmetic operations.

// #define CGAL_NT_CHECK_DEBUG(s) std::cerr << s << std::endl
#define CGAL_NT_CHECK_DEBUG(s)

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// std::equal_to<T> sticks things to one type, so I define my own.
struct Equal_to
{
  template < typename NT1, typename NT2 >
  bool operator()(const NT1& a, const NT2& b) const
  { return a == b; }
};

} // namespace CGALi


template < typename NT1, typename NT2, typename Cmp = CGALi::Equal_to >
class Number_type_checker
{
  NT1 _n1;
  NT2 _n2;
  Cmp cmp;
public:

  typedef Number_type_checker<NT1, NT2, Cmp> Self;

  // What if they differ between NT1 and NT2 ?
  typedef typename Number_type_traits<NT1>::Has_gcd        Has_gcd;
  typedef typename Number_type_traits<NT1>::Has_division   Has_division;
  typedef typename Number_type_traits<NT1>::Has_sqrt       Has_sqrt;
  typedef typename Number_type_traits<NT1>::Has_exact_ring_operations
                                            Has_exact_ring_operations;
  typedef typename Number_type_traits<NT1>::Has_exact_division
                                                           Has_exact_division;
  typedef typename Number_type_traits<NT1>::Has_exact_sqrt Has_exact_sqrt;

  Number_type_checker() {}
  Number_type_checker(int i)
    : _n1(i), _n2(i) { CGAL_assertion(is_valid()); }
  Number_type_checker(double d)
    : _n1(d), _n2(d) { CGAL_assertion(is_valid()); }
  Number_type_checker(const NT1 &n1, const NT2 &n2)
    : _n1(n1), _n2(n2) { CGAL_assertion(is_valid()); }

  // The following need to be dependant on NT1 != {NT2,int,double} ...
  //Number_type_checker(const NT1 &n1) : _n1(n1), _n2(n1) {}
  //Number_type_checker(const NT2 &n2) : _n1(n2), _n2(n2) {}

  Self& operator+=(const Self &a)
  { n1() += a.n1(); n2() += a.n2(); CGAL_assertion(is_valid()); return *this; }

  Self& operator-=(const Self &a)
  { n1() -= a.n1(); n2() -= a.n2(); CGAL_assertion(is_valid()); return *this; }

  Self& operator*=(const Self &a)
  { n1() *= a.n1(); n2() *= a.n2(); CGAL_assertion(is_valid()); return *this; }

  Self& operator/=(const Self &a)
  { n1() /= a.n1(); n2() /= a.n2(); CGAL_assertion(is_valid()); return *this; }

  // Accessors and setters.
  const NT1& n1() const { return _n1; }
  const NT2& n2() const { return _n2; }

  NT1& n1() { return _n1; }
  NT2& n2() { return _n2; }

  // Validity checking.
  bool is_valid() const
  {
    CGAL_NT_CHECK_DEBUG("Checking...");
    bool b = cmp(_n1, _n2);
    if (!b)
    {
      CGAL_NT_CHECK_DEBUG("Different values :");
      CGAL_NT_CHECK_DEBUG("n1 = " << _n1);
      CGAL_NT_CHECK_DEBUG("n2 = " << _n2);
    }
    return b;
  }
};

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() + b.n1(), a.n2() + b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(a + b.n1(), a + b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator+(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator+");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() + b, a.n2() + b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator-");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() - b.n1(), a.n2() - b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator-");
   return Number_type_checker<NT1, NT2, Cmp>(a - b.n1(), a - b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator-");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() - b, a.n2() - b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator-(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   CGAL_NT_CHECK_DEBUG("unary operator-");
   return Number_type_checker<NT1, NT2, Cmp>(-a.n1(), -a.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator*(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator*");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() * b.n1(), a.n2() * b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator*(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator*");
   return Number_type_checker<NT1, NT2, Cmp>(a * b.n1(), a * b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator*(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator*");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() * b, a.n2() * b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator/(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() / b.n1(), a.n2() / b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator/(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(a / b.n1(), a / b.n2());
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
operator/(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(a.n1() / b, a.n2() / b);
}

template < typename NT1, typename NT2, typename Cmp >
Number_type_checker<NT1, NT2, Cmp>
sqrt(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   CGAL_NT_CHECK_DEBUG("operator/");
   return Number_type_checker<NT1, NT2, Cmp>(CGAL::sqrt(a.n1()),
                                             CGAL::sqrt(a.n2()));
}


template < typename NT1, typename NT2, typename Cmp >
bool
operator==(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() == b.n1();
   bool b2 = a.n2() == b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator!=(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() != b.n1();
   bool b2 = a.n2() != b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() < b.n1();
   bool b2 = a.n2() < b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>(const Number_type_checker<NT1, NT2, Cmp> &a,
          const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() > b.n1();
   bool b2 = a.n2() > b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<=(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() <= b.n1();
   bool b2 = a.n2() <= b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>=(const Number_type_checker<NT1, NT2, Cmp> &a,
           const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = a.n1() >= b.n1();
   bool b2 = a.n2() >= b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator==(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() == i;
   bool b2 = a.n2() == i;
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator!=(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() != i;
   bool b2 = a.n2() != i;
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() < i;
   bool b2 = a.n2() < i;
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() > i;
   bool b2 = a.n2() > i;
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<=(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() <= i;
   bool b2 = a.n2() <= i;
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>=(const Number_type_checker<NT1, NT2, Cmp> &a, int i)
{
   bool b1 = a.n1() >= i;
   bool b2 = a.n2() >= i;
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator==(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i == b.n1();
   bool b2 = i == b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator!=(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i != b.n1();
   bool b2 = i != b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i < b.n1();
   bool b2 = i < b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i > b.n1();
   bool b2 = i > b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator<=(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i <= b.n1();
   bool b2 = i <= b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
operator>=(int i, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   bool b1 = i >= b.n1();
   bool b2 = i >= b.n2();
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
Sign
sign(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   Sign s1 = CGAL::sign(a.n1());
   Sign s2 = CGAL::sign(a.n2());
   CGAL_assertion(s1 == s2);
   return s1;
}

template < typename NT1, typename NT2, typename Cmp >
Comparison_result
compare(const Number_type_checker<NT1, NT2, Cmp> &a,
        const Number_type_checker<NT1, NT2, Cmp> &b)
{
   Comparison_result c1 = CGAL::compare(a.n1(), b.n1());
   Comparison_result c2 = CGAL::compare(a.n2(), b.n2());
   CGAL_assertion(c1 == c2);
   return c1;
}

template < typename NT1, typename NT2, typename Cmp >
Comparison_result
compare(int a, const Number_type_checker<NT1, NT2, Cmp> &b)
{
   Comparison_result c1 = CGAL::compare(a, b.n1());
   Comparison_result c2 = CGAL::compare(a, b.n2());
   CGAL_assertion(c1 == c2);
   return c1;
}

template < typename NT1, typename NT2, typename Cmp >
Comparison_result
compare(const Number_type_checker<NT1, NT2, Cmp> &a, int b)
{
   Comparison_result c1 = CGAL::compare(a.n1(), b);
   Comparison_result c2 = CGAL::compare(a.n2(), b);
   CGAL_assertion(c1 == c2);
   return c1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
is_finite(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   bool b1 = CGAL::is_finite(a.n1());
   bool b2 = CGAL::is_finite(a.n2());
   CGAL_assertion(b1 == b2);
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
bool
is_valid(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   bool b1 = CGAL::is_valid(a.n1());
   bool b2 = CGAL::is_valid(a.n2());
   CGAL_assertion(b1 == b2);
   // Should we also call a.is_valid() ?
   return b1;
}

template < typename NT1, typename NT2, typename Cmp >
std::pair<double, double>
to_interval(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   std::pair<double, double> i1 = CGAL::to_interval(a.n1());
   std::pair<double, double> i2 = CGAL::to_interval(a.n2());
   // Here we could check that there is a common point.
   // CGAL_assertion( ??? );

   // We return one of the two.
   return i1;
}

template < typename NT1, typename NT2, typename Cmp >
double
to_double(const Number_type_checker<NT1, NT2, Cmp> &a)
{
   double d1 = CGAL::to_double(a.n1());
   double d2 = CGAL::to_double(a.n2());
   // What can we check ?
   // We return one of the two.
   return d1;
}


// What to do with the IO ?
// - either pick one as the main, and convert to the second.
// - output/input pairs, and read both (=> problems with FP values).

template < typename NT1, typename NT2, typename Cmp >
std::ostream &
operator<< (std::ostream & os, const Number_type_checker<NT1, NT2, Cmp> &b)
{
  return os << b.n1();
}

template < typename NT1, typename NT2, typename Cmp >
std::istream &
operator>> (std::istream & is, Number_type_checker<NT1, NT2, Cmp> &b)
{
  is >> b.n1();
  b.n2() = b.n1(); // We hope that there is a conversion.
  return is;
}

template < typename NT1, typename NT2, typename Cmp >
io_Operator
io_tag(const Number_type_checker<NT1, NT2, Cmp> &)
{
  // Not sure what do to here.
  return io_Operator();
}

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_TYPE_CHECKER_H
