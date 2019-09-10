/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * File: Expr.h
 * Synopsis: a class of Expression in Level 3
 * 
 * Written by 
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu> 
 *       Vikram Sharma<sharma@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0+
 ***************************************************************************/

// We need to include BigFloat.h here because there is a circular dependency
// between Expr and BigFloat.
#include <CGAL/CORE/BigFloat.h>

#ifndef _CORE_EXPR_H_
#define _CORE_EXPR_H_

#include <CGAL/CORE/ExprRep.h>
#include <CGAL/assertions.h>

namespace CORE { 

/// \class Expr Expr.h
/// \brief Expr is a class of Expression in Level 3
typedef RCImpl<ExprRep> RCExpr;

class Expr : public RCExpr {
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  Expr() : RCExpr(new ConstDoubleRep()) {}
  /// constructor for <tt>int</tt>
  Expr(int i) : RCExpr(new ConstDoubleRep(i)) {}
  /// constructor for <tt>short</tt>
  Expr(short i) : RCExpr(new ConstDoubleRep(i)) {}
  /// constructor for <tt>unsigned int</tt>
  Expr(unsigned int ui) : RCExpr(new ConstDoubleRep(ui)) {}

  /// constructor for <tt>long</tt>
  Expr(long l) : RCExpr(new ConstRealRep(Real(l))) {}
  /// constructor for <tt>unsigned long</tt>
  Expr(unsigned long ul) : RCExpr(new ConstRealRep(Real(ul))) {}

  /// constructor for <tt>float</tt>
  /** \note the results of this constructor may appear unpredictable to the 
   *  user.  E.g.,  one may assume that new Expr(.1) is exactly equal to .1,
   *  but it will be print as
   *      .1000000000000000055511151231257827021181583404541015625.
   *  This is so because .1 cannot be represented exactly as a double
   *  (or, for that matter, as a binary fraction of any finite length). 
   *  The value is the closest double value determined by the compiler.
   */
  Expr(float f) : RCExpr(NULL) { // check for valid numbers
    // (i.e., not infinite and not NaN)
    if (! CGAL_CORE_finite(f)) {
      core_error(" ERROR : constructed an invalid float! ", __FILE__, __LINE__, false);
      if (get_static_AbortFlag())
        abort();
      get_static_InvalidFlag() = -1;
    }
    rep = new ConstDoubleRep(f);
  }
  /// constructor for <tt>double</tt>
  Expr(double d) : RCExpr(NULL) { // check for valid numbers
    // (i.e., not infinite and not NaN)
    if (! CGAL_CORE_finite(d)) {
      core_error(" ERROR : constructed an invalid double! ", __FILE__, __LINE__, false);
      if (get_static_AbortFlag())
        abort();
      get_static_InvalidFlag() = -2;
    }
    rep = new ConstDoubleRep(d);
  }

  /// constructor for <tt>BigInt</tt>
  Expr(const BigInt& I) : RCExpr(new ConstRealRep(Real(I))) {}
  /// constructor for <tt>BigRat</tt>
  Expr(const BigRat& R) : RCExpr(new ConstRealRep(Real(R))) {}

  /// constructor for <tt>BigFloat</tt>
  Expr(const BigFloat& F) : RCExpr(new ConstRealRep(Real(F))) {}

  /// constructor for <tt>const char*</tt>
  /** construct Expr from a string representation \a s
   * with precision \a prec. It is perfectly predictable:
   * new Expr(".1") is exactly equal to .1, as one would expect. Therefore,
   * it is generally recommended that the (String) constructor be used in
   * preference to the (double) constructor.
   */
  Expr(const char *s, const extLong& p = get_static_defInputDigits())
      : RCExpr(new ConstRealRep(Real(s, p))) {}

  /// constructor for <tt>std::string</tt>
  Expr(const std::string& s, const extLong& p = get_static_defInputDigits())
      : RCExpr(new ConstRealRep(Real(s, p))) {}

  /// constructor for <tt>Real</tt>
  Expr(const Real &r) : RCExpr(new ConstRealRep(r)) {}

  /// constructor for Polynomial node (n-th root)
  /** default value n=0 means the first positive root */
  template <class NT>
  Expr(const Polynomial<NT>& p, int n = 0)
      : RCExpr(new ConstPolyRep<NT>(p, n)) {}

  /// constructor for Polynomial node (root in Interval <tt>I</tt>)
  template <class NT>
  Expr(const Polynomial<NT>& p, const BFInterval& I)
      : RCExpr(new ConstPolyRep<NT>(p, I)) {}

  /// constructor for ExprRep
  Expr(ExprRep* p) : RCExpr(p) {}
  //@}

  /// \name Copy-Assignment-Destructors
  //@{
  /// copy constructor
  Expr(const Expr& rhs) : RCExpr(rhs) {
    rep->incRef();
  }

  /// = operator
  Expr& operator=(const Expr& rhs) {
    if (this != &rhs) {
      rep->decRef();
      rep = rhs.rep;
      rep->incRef();
    }
    return *this;
  }
  /// destructor
  ~Expr() {
    rep->decRef();
  }
  //@}

  /// \name Compound Assignment Operators
  //@{
  /// += operator
  Expr& operator+=(const Expr& e) {
    *this = new AddRep(rep, e.rep);
    return *this;
  }
  /// -= operator
  Expr& operator-=(const Expr& e) {
    *this = new SubRep(rep, e.rep);
    return *this;
  }
  /// *= operator
  Expr& operator*=(const Expr& e) {
    *this = new MultRep(rep, e.rep);
    return *this;
  }
  /// /= operator
  Expr& operator/=(const Expr& e) {
    if ((e.rep)->getSign() == 0) {
      core_error(" ERROR : division by zero ! ",__FILE__, __LINE__, false);
      if (get_static_AbortFlag())
        abort();
      get_static_InvalidFlag() = -3;
    }
    *this = new DivRep(rep, e.rep);
    return *this;
  }
  //@}

  /// \name Unary Minus, Increment and Decrement Operators
  //@{
  /// unary plus
  Expr operator+() const {
    return Expr(*this);
  }
  /// unary minus
  Expr operator-() const {
    return Expr(new NegRep(rep));
  }
  /// left increment operator (++i)
  Expr& operator++() {
    *this += 1;
    return *this;
  }
  /// right increment operator (i++)
  Expr operator++(int) {
    Expr t(*this);
    *this += 1;
    return t;
  }
  /// left decrement operator (--i)
  Expr& operator--() {
    *this -= 1;
    return *this;
  }
  /// right deccrement operator (i--)
  Expr operator--(int) {
    Expr t(*this);
    *this -= 1;
    return t;
  }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  void fromString(const char* s, const extLong& prec = get_static_defInputDigits()) {
    *this = Expr(s, prec);
  }
  /// convert to <tt>std::string</tt>
  /** give decimal string representation */
  std::string toString(long prec=get_static_defOutputDigits(), bool sci=false) const {
    return rep->toString(prec, sci);
  }
  //@}
  //

  /// \name Conversion Functions
  //@{
  /// convert to \c int
  int intValue() const {
    return approx(64, 1024).intValue();
  }
  /// convert to \c long
  long longValue() const {
    return approx(64, 1024).longValue();
  }
  /// convert to \c float
  float floatValue() const {
    return approx(53, 1024).floatValue();
  }
  /// convert to \c double
  /** chen: - use equivalent precision (rel:53, abs: 1024)
    as in IEEE double. enforce an evaluation in case
    before this has been done before casting. */
  double doubleValue() const {
    return approx(53, 1024).doubleValue();
  }
  /// convert to an interval defined by a pair of \c double
  /** If value is exact, the two \c double will coincide
   */
  CGAL_CORE_EXPORT void doubleInterval(double & lb, double & ub) const;
  /// convert to \c BigInt (approximate it first!)
  BigInt BigIntValue() const {
    return rep->BigIntValue();
  }
  /// convert to \c BigRat (approximate it first!)
  BigRat BigRatValue() const {
    return rep->BigRatValue();
  }
  /// convert to \c BigFloat (approximate it first!)
  /** Ought to allow BigFloatValue() take an optional precision argument */
  BigFloat BigFloatValue() const {
    return rep->BigFloatValue();
  }
  //@}

  /// \name Approximation Function
  //@{
  /// Compute approximation to combined precision [\a r, \a a].
  /** Here is the definition of what this means:
       If e is the exact value and ee is the approximate value,
       then  |e - ee| <= 2^{-a} or  |e - ee| <= 2^{-r} |e|. */
  const Real & approx(const extLong& relPrec = get_static_defRelPrec(),
                      const extLong& absPrec = get_static_defAbsPrec()) const {
    return rep->getAppValue(relPrec, absPrec);
  }
  //@}

  /// \name Helper Functions
  //@{
  //CONSTANTS:
  /// return Expr(0)
  CGAL_CORE_EXPORT static const Expr& getZero();

  /// return Expr(1)
  CGAL_CORE_EXPORT static const Expr& getOne();

  /// Has Exact Division
  static bool hasExactDivision() {
    return true;
  }

  /// get the sign
  int sign() const {
    return rep->getSign();
  }
  /// is zero?
  bool isZero() const {
    return sign() == 0;
  }
  /// absolute value
  Expr abs() const {
    return (sign() >= 0) ? +(*this) : -(*this);
  }

  /// compare function
  int cmp(const Expr& e) const {
    return rep == e.rep ? 0 : SubRep(rep, e.rep).getSign();
  }

  /// return the internal representation
  ExprRep* Rep() const {
    return rep;
  }
  /// get exponent of current approximate value
  long getExponent() const {
    return BigFloatValue().exp();
  }
  /// get mantissa of current approximate value
  BigInt getMantissa() const {
    return BigFloatValue().m();
  }
  //@}

public:
  /// \name Debug Helper Function
  //@{
  /// debug function
  void  debug(int mode = TREE_MODE, int level = DETAIL_LEVEL,
              int depthLimit = INT_MAX) const;
  //@}
  /// debug information levels
  enum {LIST_MODE, TREE_MODE, SIMPLE_LEVEL, DETAIL_LEVEL};
};// class Expr

#define CORE_EXPR_ZERO Expr::getZero()

/// I/O Stream operator<<
inline std::ostream& operator<<(std::ostream& o, const Expr& e) {
  o << *(const_cast<ExprRep*>(&e.getRep()));
  return o;
}
/// I/O Stream operator>>
inline std::istream& operator>>(std::istream& i, Expr& e) {
  Real rVal;
  i >> rVal; // precision is = get_static_defInputDigits()
  if (i)
    e = rVal;		// only assign when reading is successful.
  return i;
}

/// floor function
CGAL_CORE_EXPORT BigInt floor(const Expr&, Expr&);
/// power function
CGAL_CORE_EXPORT Expr pow(const Expr&, unsigned long);

/// addition
inline Expr operator+(const Expr& e1, const Expr& e2) {
  return Expr(new AddRep(e1.Rep(), e2.Rep()));
}
/// substraction
inline Expr operator-(const Expr& e1, const Expr& e2) {
  return Expr(new SubRep(e1.Rep(), e2.Rep()));
}
/// multiplication
inline Expr operator*(const Expr& e1, const Expr& e2) {
  return Expr(new MultRep(e1.Rep(), e2.Rep()));
}
/// division
inline Expr operator/(const Expr& e1, const Expr& e2) {
  if (e2.sign() == 0) {
    core_error(" ERROR : division by zero ! ", __FILE__, __LINE__, false);
    if (get_static_AbortFlag())
      abort();
    get_static_InvalidFlag() = -4;
  }
  return Expr(new DivRep(e1.Rep(), e2.Rep()));
}
/// modulo operator
inline Expr operator%(const Expr& e1, const Expr& e2) {
  Expr result;
  floor(e1/e2, result);
  return result;
}

/// operator ==
/** this is inefficient if you compare to zero:
 *  e.g., if (e != 0) {...} use e.isZero() instead */
inline bool operator==(const Expr& e1, const Expr& e2) {
  return e1.cmp(e2) == 0;
}
/// operator !=
inline bool operator!=(const Expr& e1, const Expr& e2) {
  return e1.cmp(e2) != 0;
}
/// operator <
inline bool operator< (const Expr& e1, const Expr& e2) {
  return e1.cmp(e2) < 0;
}
/// operator <=
inline bool operator<=(const Expr& e1, const Expr& e2) {
  return e1.cmp(e2) <= 0;
}
/// operator <
inline bool operator> (const Expr& e1, const Expr& e2) {
  return e1.cmp(e2) > 0;
}
/// operator >=
inline bool operator>=(const Expr& e1, const Expr& e2) {
  return e1.cmp(e2) >= 0;
}

/// return sign
inline int sign(const Expr& e) {
  return e.sign();
}
/// is zero?
inline bool isZero(const Expr& e) {
  return e.isZero();
}
/// compare
/** compare two Expr \a e1 and \a e2, return
 * \retval -1 if e1 < e2,
 * \retval 0 if e1 = e2,
 * \retval 1 if e1 > e2. */
inline int cmp(const Expr& e1, const Expr& e2) {
  return e1.cmp(e2);
}
/// absolute value
inline Expr abs(const Expr& x) {
  return x.abs();
}
/// absolute value (same as abs)
inline Expr fabs(const Expr& x) {
  return abs(x);
}
/// floor
inline BigInt floor(const Expr& e) {
  Expr tmp;
  return floor(e, tmp);
}
/// ceiling
inline BigInt ceil(const Expr& e) {
  return -floor(-e);
}
/// floorLg
inline long floorLg(const Expr& e) {
  Expr tmp;
  return floorLg(floor(e));
}
/// ceilLg
inline long ceilLg(const Expr& e) {
  Expr tmp;
  return ceilLg(ceil(e));
}
/// power
inline Expr power(const Expr& e, unsigned long p) {
  return pow(e, p);
}

/// divisibility predicate
/** We do not check if e2 is 0.
 * */
// NOTE:  The name "isDivisible" is not consistent
// 		with the analogous "divisible" predicate in BigInt!
inline bool isDivisible(const Expr& e1, const Expr& e2) {
  Expr result;
  floor(e1/e2, result);
  return (result.sign() == 0);
}

/// square root
inline Expr sqrt(const Expr& e) {
  if (e.sign() < 0) {
    core_error(" ERROR : sqrt of negative value ! ", __FILE__, __LINE__, false);
    if (get_static_AbortFlag())
      abort();
    get_static_InvalidFlag() = -5;
  }
  return Expr(new SqrtRep(e.Rep()));
}

//Following two have been added to make NT=Expr work for Polynomial<NT>
/// gcd
inline Expr gcd(const Expr& /*a*/, const Expr& /*b*/) {
  return Expr(1);
}
inline Expr div_exact(const Expr& x, const  Expr& y) {
  return x/y - x%y;
}

/// helper function for constructing Polynomial node (n-th node)
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, int n = 0) {
  return Expr(p, n);
}

/// helper function for constructing Polynomial node witb BFInterval
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, const BFInterval& I) {
  return Expr(p, I);
}
/// helper function for constructing Polynomial node with pair of BigFloats
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, const BigFloat& x,
		const BigFloat& y) {
  return Expr(p, BFInterval(x, y) );
}
/// helper function for constructing Polynomial node with pair of doubles
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, double x, double y) {
  return Expr(p, BFInterval(BigFloat(x), BigFloat(y)) );
}
/// helper function for constructing Polynomial node with pair of ints
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, int x, int y) {
  return Expr(p, BFInterval(BigFloat(x), BigFloat(y)) );
}

/// constructor for Polynomial node of the form x^m - n (i.e., radicals)
/** We assume that n >= 0 and m >= 1
 * */
template <class NT>
inline Expr radical(const NT& n, int m) {
  CGAL_assertion(n>=0 && m>=1);
  if (n == 0 || n == 1 || m == 1)
    return Expr(n);
  Polynomial<NT> Q(m);
  Q.setCoeff(0, -n);
  Q.setCoeff(m, 1);
  return Expr(Q, 0);
}

// We include this file here and not from inside Poly.h,
// because otherwise VC++.net2003 can't compile Expr.cpp
#include <CGAL/CORE/poly/Poly.tcc>

} //namespace CORE

#ifdef CGAL_HEADER_ONLY
#include <CGAL/CORE/Expr_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // _CORE_EXPR_H_
