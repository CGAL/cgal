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
 * File: Poly.h
 *  
 * Description: simple polynomial class
 * 
 *	REPRESENTATION:
 *	--Each polynomial has a nominal "degree" (this
 *		is an upper bound on the true degree, which
 *		is determined by the first non-zero coefficient).
 *	--coefficients are parametrized by some number type "NT".
 *	--coefficients are stored in the "coeff" array of
 *		length "degree + 1".  
 *		CONVENTION: coeff[i] is the coefficient of X^i.  So, a
 *			    coefficient list begins with the constant term.
 *	--IMPORTANT CONVENTION:
 *		the zero polynomial has degree -1
 *		while nonzero constant polynomials have degree 0.
 * 
 *	FUNCTIONALITY:
 *	--Polynomial Ring Operations (+,-,*)
 *	--Power
 *	--Evaluation
 *	--Differentiation
 *	--Remainder, Quotient 
 *      --GCD
 *	--Resultant, Discriminant (planned)
 *	--Polynomial Composition (planned)
 *	--file I/O (planned)
 *	
 * Author: Chee Yap 
 * Date:   May 28, 2002
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 * SPDX-License-Identifier: LGPL-3.0+
 ***************************************************************************/

#ifndef CORE_POLY_H
#define CORE_POLY_H

#include <CGAL/CORE/BigFloat.h>
#include <CGAL/CORE/Promote.h>
#include <vector>
#include <CGAL/assertions.h>
#include <CGAL/tss.h>

namespace CORE { 
using namespace std;
class Expr;
// ==================================================
// Typedefs
// ==================================================

//typedef std::vector<Expr>	VecExpr;
//typedef std::pair<Expr, Expr>	Interval;
//typedef std::vector<Interval>	VecInterval;
typedef std::pair<BigFloat, BigFloat>	BFInterval;
// NOTE: an error condition is indicated by
// the special interval (1, 0)
typedef std::vector<BFInterval>	BFVecInterval;


// ==================================================
// Polynomial Class
// ==================================================

template <class NT>
class Polynomial {
 private:
  //The following are used in the constructor from strings.
  //For more details see the related constructor.

public:
  typedef std::vector<NT> VecNT;
  typedef NT coeffType;

  int degree;	// This is the nominal degree (an upper bound
  // on the true degree)
  NT * coeff;	// coeff is an array of size degree+1;
  //	This remark holds even when degree = -1.
  // Notes:
  // (1) coeff[i] is the coefficient of x^i
  // (2) The Zero Polynomial has degree -1
  // (3) Nonzero Constant Polynomials has degree 0

  // STATIC MEMBERS
  // static NT ccc_; // THIS IS A TEMPORARY HACK
  static const int COEFF_PER_LINE = 4;		// pretty print parameters
  static const char INDENT_SPACE[];  		// pretty print parameters
  static const Polynomial<NT> & polyZero();
  static const Polynomial<NT> & polyUnity();

  // Constructors:
  Polynomial(void);	// the Zero Polynomial
  Polynomial(int n);	// construct the Unit Polynomial of nominal deg n>=0
  Polynomial(int n, const NT * coef);
  Polynomial(const Polynomial &);
  Polynomial(const VecNT &);
  Polynomial(int n, const char* s[]);
  Polynomial(const string & s, char myX='x');
  Polynomial(const char* s, char myX='x');
  ~Polynomial();

  private:
  void constructX(int n, Polynomial<NT>& P);
  void constructFromString(string & s, char myX='x');
  int getnumber(const char* c, int i, unsigned int len, Polynomial<NT> & P);
  bool isint(char c);
  int getint(const char* c, int i, unsigned int len, int & n);
  int matchparen(const char* cstr, int start);
  int getbasicterm(string & s, Polynomial<NT> & P);
  int getterm(string & s, Polynomial<NT> & P);


  public:
  //Returns a Polynomial corresponding to s, which is supposed to
  //contain as place-holders the chars 'x' and 'y'.
  Polynomial<NT> getpoly(string & s);

  // Assignment:
  Polynomial & operator=(const Polynomial&);

  // Expand and Contract
  //  -- they are semi-inverses: i.e., Contract(expand(p))=p
  int expand(int n);	// Change the current degree to n
  // Helper function for polynomial arithmetic
  int contract();	// get rid of leading zeros

  // Polynomial arithmetic (these are all self-modifying):
  Polynomial & operator+=(const Polynomial&);	// +=
  Polynomial & operator-=(const Polynomial&);	// -=
  Polynomial & operator*=(const Polynomial&);	// *=
  Polynomial & operator-();			// unary minus
  Polynomial & power (unsigned int n) ;		// power (*this is changed!)

  Polynomial & mulScalar (const NT & c);	// return (*this) * (c)
  Polynomial & mulXpower(int i); // If i >= 0, then this is equivalent
                                 // to multiplying by X^i.
                                 // If i < 0 to dividing by X^i
  Polynomial pseudoRemainder (const Polynomial& B, NT& C); // C = return value
  Polynomial pseudoRemainder (const Polynomial& B);        // no C version
  // The pseudo quotient of (*this) mod B
  //	is returned, but (*this) is transformed
  //	into the pseudo remainder.  If the argument C is provided,
  //	Then C*(*this) = B*pseudo-quotient + pseudo-remainder.
  Polynomial & negPseudoRemainder (const Polynomial& B);  // negative remainder
  Polynomial reduceStep (const Polynomial& B ); //helper for pseudoRemainder
  // One step of pseudo remainder
  // What is returned is a special polynomial C + X*M  (not "C+M")
  //    telling us the initial constant C and
  //    the quotient M of C*(THIS) divided by p.
  Polynomial testReduceStep(const Polynomial& A, const Polynomial& B); //helper

  // Get functions
  int getDegree() const;        // nominal degree
  int getTrueDegree() const;    // true degree
  NT getCoeffi(int i) const;
  const NT & getLeadCoeff() const;      // get TRUE leading coefficient
  const NT & getTailCoeff() const;      // get last non-zero coefficient
  NT** getCoeffs() ;		// get all coefficients
  const NT& getCoeff(int i) const;      // Get single coefficient of X^i
                                        // NULL pointer if invalid i
  // Set functions
  bool setCoeff(int i, const NT & cc);  // Make cc the coefficient of X^i
                                        // Return FALSE if invalid i
                                        // !! User's responsibility to
                                        // delete the old coefficient if
                                        // necessary !!
  // Helper Functions:
  /// Reverse reverses the coefficients
  void reverse();		
  /// Negation of a polynomial (multiplication by -1)
  /// Useful for Sturm
  Polynomial & negate();	
  /// Suppressing Zero Roots
  /// It amounts to dividing (*this) by X^k, so that the
  /// the tail coeff is non-zero. Returns the value of k.
  int makeTailCoeffNonzero();

  // Evaluation Functions:
  /// Polynomial evaluation where the coefficients are approximated first
  /// Returns a BigFloat with error that contains the value
  BigFloat evalApprox(const BigFloat& f, 
    const extLong& r=get_static_defRelPrec(), 
    const extLong& a=get_static_defAbsPrec()) const;
  /// Polynomial evaluation at a BigFloat value.
  /// The returned BigFloat (with error) has the exact sign.  
  /// In particular, if the value is 0, we return 0.
  /// @param oldMSB is any estimate of the negative log of the evaluation
  BigFloat evalExactSign(const BigFloat& val, const extLong& oldMSB=54) const;
  /// Polynomial evaluation that return the same type as its argument
  /// Caution: The type T must be greater or equal to the type NT
  /// 	NOTE: Eventually, we will remove this restriction by
  /// 	introduce MaxType(NT,T) for the return type.
  template <class T>
  MAX_TYPE(NT, T) eval(const T&) const;	

  // Bounds
  BigFloat CauchyUpperBound() const;  // Cauchy Root Upper Bound
  BigFloat CauchyLowerBound() const;  // Cauchy Root Lower Bound
  BigInt CauchyBound() const;  // Cauchy Root Bound from Erich Kaltofen
  BigInt UpperBound() const;  // Another Cauchy Root Bound; an improvement over
                               //Erich Kaltofen
  BigFloat sepBound() const;	// separation bound (multiple roots allowed)
  BigFloat height() const;	// height return type BigFloat
  BigFloat length() const;	// length return type BigFloat

  // Differentiation
  Polynomial & differentiate() ;		// self-differentiation
  Polynomial & differentiate(int n) ;		// multi self-differentiation

  // Reductions of polynomials (NT must have gcd function)
  Polynomial sqFreePart(); // Square free part of P is P/gcd(P,P'). Return gcd
  Polynomial & primPart();   // Primitive Part of *this (which is changed)

  //////////////////////////////////////////////////////////////////
  // Resultant and discriminant
  // NT & resultant() ;				// resultant
  // NT & discriminant() ;			// discriminant

  // Composition of Polynomials:
  // NOT yet implemented

  //////////////////////////////////////////////////////////////////
  // Polynomial Dump
  void dump(std::ofstream & ofs, std::string msg="",
         std::string com="# ", std::string com2="# ") const;  // dump to file
  void dump(std::string msg="", std::string com="# ",
	 std::string com2="# ") const; // dump to cout
  void filedump(std::ostream & os, std::string msg="", std::string com="# ",
	 std::string com2="# ") const; // dump workhorse (called by dump())
  void mapleDump() const;              // dump of maple code for Polynomial

}; //Polynomial Class

// template < class NT >
//     NT Polynomial<NT>::ccc_;

// ==================================================
// Static Constants
//	Does this belong here?
// ==================================================

template < class NT >
CORE_INLINE
const Polynomial<NT> & Polynomial<NT>::polyZero() {
  CGAL_STATIC_THREAD_LOCAL_VARIABLE_0(Polynomial<NT>, zeroP);
  return zeroP;
}

template < class NT >
CORE_INLINE
const Polynomial<NT> & Polynomial<NT>::polyUnity() {
  static const NT c[] = {1};
  CGAL_STATIC_THREAD_LOCAL_VARIABLE_2(Polynomial<NT>, unityP, 0, c);
  return unityP;
}

// ==================================================
// Useful functions for Polynomial class
// ==================================================

// polynomial arithmetic:
template < class NT >
Polynomial<NT> operator+(const Polynomial<NT>&, const Polynomial<NT>&);// +
template < class NT >
Polynomial<NT> operator-(const Polynomial<NT>&, const Polynomial<NT>&);// -
template < class NT >
Polynomial<NT> operator*(const Polynomial<NT>&, const Polynomial<NT>&);// *
template < class NT >
Polynomial<NT> power(const Polynomial<NT>&, int n);		// power
template < class NT >
Polynomial<NT> differentiate(const Polynomial<NT>&);		// differentiate
template < class NT >
Polynomial<NT> differentiate(const Polynomial<NT>&, int n);	// multi-differ.
//Content of a Polynomial
template < class NT >
NT content(const Polynomial<NT>& p);

template <class NT>
bool isDivisible(Polynomial<NT> p, Polynomial<NT> q);

// GCD of two polynomials
template < class NT >
Polynomial<NT> gcd(const Polynomial<NT>& p, const Polynomial<NT>& q);

//Resultant of two polynomials
template < class NT >
NT res( Polynomial<NT> p,  Polynomial<NT> q);

//Principal Subresultant Coefficient (psc) of two polynomials
template < class NT >
NT psc(int i, Polynomial<NT> p,  Polynomial<NT> q);

//Returns the polynomial which contains only the real roots
//of P which have multiplicity d
template < class NT >
Polynomial<NT> factorI(Polynomial<NT> p, int d);

// comparisons
template < class NT >
bool operator==(const Polynomial<NT>&, const Polynomial<NT>&); // ==
template < class NT >
bool operator!=(const Polynomial<NT>&, const Polynomial<NT>&); // !=
template < class NT >
bool zeroP(const Polynomial <NT>&);			// =Zero Poly?
template < class NT >
bool unitP(const Polynomial <NT>&);			// =Unit Poly?

// stream i/o
template < class NT >
std::ostream& operator<<(std::ostream&, const Polynomial<NT>&);
template < class NT >
std::istream& operator>>(std::istream&, Polynomial<NT>&);

// ==================================================
// Inline Functions
// ==================================================

// friend polynomial arithmetic:
template < class NT >
CORE_INLINE
Polynomial<NT> operator+(const Polynomial<NT>& p1,
                         const Polynomial<NT>& p2) {	// +
  return Polynomial<NT>(p1) += p2;
}
template < class NT >
CORE_INLINE
Polynomial<NT> operator-(const Polynomial<NT>& p1,
                         const Polynomial<NT>& p2) {	// -
  return Polynomial<NT>(p1) -= p2;
}
template < class NT >
CORE_INLINE
Polynomial<NT> operator*(const Polynomial<NT>& p1,
                         const Polynomial<NT>& p2) {	// *
  return Polynomial<NT> (p1) *= p2;
}
template < class NT >
CORE_INLINE
Polynomial<NT> power(const Polynomial<NT>& p, int n) {				// power
  return Polynomial<NT>(p).power(n);
}

// equal to zero poly?
template < class NT >
CORE_INLINE
bool zeroP(const Polynomial <NT>& p) {			// =Zero Poly?
  return (p.getTrueDegree()== -1);
}
template < class NT >
CORE_INLINE
bool unitP(const Polynomial <NT>& p) {			// =Unit Poly?
  int d = p.getTrueDegree();
  return ((d == 0) && p.coeff[0]==1 );
}

// get functions
template < class NT >
CORE_INLINE
int Polynomial<NT>::getDegree() const {
  return degree;
}
// get TRUE leading coefficient
template < class NT >
CORE_INLINE
const NT & Polynomial<NT>::getLeadCoeff() const {
  return getCoeff(getTrueDegree());
}

// get last non-zero coefficient
template < class NT >
CORE_INLINE
const NT & Polynomial<NT>::getTailCoeff() const {
  for (int i = 0; i<= getTrueDegree(); i++)
    if (coeff[i] != 0)
      return coeff[i];
  // This ought to be an error (user should check this) :
  NT * zero = new NT(0);
  return *zero;
}

template < class NT >
CORE_INLINE
NT** Polynomial<NT>::getCoeffs() {
  return &coeff;
}
template < class NT >
CORE_INLINE
const NT& Polynomial<NT>::getCoeff(int i) const {
  //if (i > degree) return NULL;
  CGAL_assertion(i <= degree);
  return coeff[i];
}
// set functions
template < class NT >
CORE_INLINE
bool Polynomial<NT>::setCoeff(int i, const NT & cc) {
  if ((i<0) || (i > degree))
    return false;
  coeff[i] = cc;
  return true;
}

// IMPLEMENTATIONS ARE FOUND IN
//#include <CGAL/CORE/poly/Poly.tcc>
//
// We include this file from CORE/Expr.h, AFTER the definition
// of class Expr, because otherwise VC++.net2003 can'y compile Expr.cpp

} //namespace CORE
#endif
