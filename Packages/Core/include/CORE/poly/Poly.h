/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project

 File: Poly.h

 Description: simple polynomial class

	REPRESENTATION:
	--Each polynomial has a nominal "degree" (this
		is an upper bound on the true degree, which
		is determined by the first non-zero coefficient).
	--coefficients are parametrized by some number type "NT".
	--coefficients are stored in the "coeff" array of
		length "degree + 1".  
		CONVENTION: coeff[i] is the coefficient of X^i.  So, a
			    coefficient list begins with the constant term.
	--IMPORTANT CONVENTION:
		the zero polynomial has degree -1
		while nonzero constant polynomials have degree 0.

	FUNCTIONALITY:
	--Polynomial Ring Operations (+,-,*)
	--Power
	--Evaluation
	--Differentiation
	--Remainder, Quotient 
	--Resultant, Discriminant (planned)
	--Polynomial Composition (planned)
	--file I/O (planned)
	
 Author: Chee Yap 
 Date:   May 28, 2002
 Core Library 1.4.1
 $Id$
 ************************************** */ 

#ifndef CORE_POLY_H
#define CORE_POLY_H

#include "../BigFloat.h"
#include <vector>

CORE_BEGIN_NAMESPACE

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
public:
  typedef std::vector<NT> VecNT;

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
  static  int COEFF_PER_LINE;		// pretty print parameters
  static const char * INDENT_SPACE;		// pretty print parameters

  static const Polynomial<NT> & polyZero();
  static const Polynomial<NT> & polyUnity();
  static Polynomial polyWilkinson;	      // a sample polynomial

  // Constructors:
  Polynomial(void);	// the Zero Polynomial
  Polynomial(int n);	// construct the Unit Polynomial of nominal deg n>=0
  Polynomial(int n, NT * coef);
  Polynomial(const Polynomial &);
  Polynomial(const VecNT &);
  Polynomial(int n, const char* s[]);
  ~Polynomial();

  // Assignment:
  Polynomial & operator=(const Polynomial&);

  // Expand and Contract  
  //  -- they are semi-inverses: i.e., Contract(expand(p))=p
  int expand(int n);	// Change the current degree to n
			// Helper function for polynomial arithmetic
  int contract();		// get rid of leading zeros

  // Polynomial arithmetic (these are all self-modifying):
  Polynomial & operator+=(const Polynomial&);	// +=
  Polynomial & operator-=(const Polynomial&);	// -=
  Polynomial & operator*=(const Polynomial&);	// *=
  Polynomial & operator-();			// unary minus
  Polynomial & power (unsigned int n) ;		// power (*this is changed!)

  Polynomial & mulScalar (const NT & c);	// return (*this) * (c)
  Polynomial & mulXpower(int i); 	// if i >= 0, then this is equivalent
					//	to multiplying by X^i
					// if i < 0 to dividing by X^i
  Polynomial pseudoRemainder (Polynomial& p);
			// the pseudo quotient of (*this) mod p
			//	is returned, but (*this) is transformed
			//	into the pseudo remainder
  Polynomial reduceStep (Polynomial& p );
			// One step of pseudo remainder 
			// What is returned is a special polynomial C+M
			//	telling us the initial constant C and					//	the quotient M of C*(THIS) divided by p.
			
  // Polynomial gcd(Polynomial p); // This may not be defined for some NT domains

  // Get functions
  int getDegree() const;	// nominal degree
  int getTrueDegree() const;	// true degree
  const NT & getLeadCoeff() const;	// get TRUE leading coefficient
  const NT & getTailCoeff() const;	// get last non-zero coefficient
  NT** getCoeffs() ;		// get all coefficients
  const NT& getCoeff(int i) const;	// Get single coefficient of X^i
				// NULL pointer if invalid i
  // Set functions
  bool setCoeff(int i, const NT & cc);	// Make cc the coefficient of X^i
					// Return FALSE if invalid i
					// !! User's responsibility to 
					// delete the old coefficient if
					// necessary !!
  // Helper Functions
  void reverse();		// reverse the coefficients;  useful when
				// input coefficients are in reversed

  // Evaluation
  Expr eval(const Expr&) const;		// evaluation
  BigFloat eval(const BigFloat&) const;	// evaluation

  // Bounds
  BigFloat CauchyUpperBound() const;  // Cauchy Root Upper Bound
  BigFloat CauchyLowerBound() const;  // Cauchy Root Lower Bound
  BigFloat sepBound() const;	// separation bound (multiple roots allowed)
  BigFloat height() const;	// height return type BigFloat
  BigFloat length() const;	// length return type BigFloat

  // Differentiation
  Polynomial & differentiate() ;		// self-differentiation
  Polynomial & differentiate(int n) ;		// multi self-differentiation

  // Reductions of polynomials
  Polynomial & squareFreePart(); // P/gcd(P,P')
  Polynomial & primPart(); 	//  Primitive Part of *this (which is changed)

  //////////////////////////////////////////////////////////////////
  // Resultant and discriminant
  // NT & resultant() ;				// resultant
  // NT & discriminant() ;				// discriminant

  // Composition of Polynomials:
  // NOT yet implemented

  //////////////////////////////////////////////////////////////////
  // Polynomial Dump
  void dump() const;			// plain dump
  void dump(std::string m1) const;	// dump with message string
  void mapleDump() const;		// dump of maple code for Polynomial

}; //Polynomial Class

// template < class NT >
//     NT Polynomial<NT>::ccc_;

// ==================================================
// Static Constants
//	Does this belong here?
// ==================================================

template < class NT >
  CORE_INLINE
  const Polynomial<NT> & Polynomial<NT>::polyZero(){
	static Polynomial<NT> zeroP;
	return zeroP;
  }

template < class NT >
  CORE_INLINE
  const Polynomial<NT> & Polynomial<NT>::polyUnity() {
	static NT c[] = {1};
	static Polynomial<NT> unityP(0, c);
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
				const Polynomial<NT>& p2){	// +
    return Polynomial<NT>(p1) += p2;
  }
template < class NT >
  CORE_INLINE
  Polynomial<NT> operator-(const Polynomial<NT>& p1,
				const Polynomial<NT>& p2){	// -
    return Polynomial<NT>(p1) -= p2;
  }
template < class NT >
  CORE_INLINE
  Polynomial<NT> operator*(const Polynomial<NT>& p1,
				const Polynomial<NT>& p2){	// *
    return Polynomial<NT> (p1) *= p2;
  }
template < class NT >
  CORE_INLINE
  Polynomial<NT> power(const Polynomial<NT>& p, int n){				// power
    return Polynomial<NT>(p).power(n);
  }

// equal to zero poly?
template < class NT >
  CORE_INLINE
  bool zeroP(const Polynomial <NT>& p){			// =Zero Poly?
	return (p.getTrueDegree()== -1);
  }
template < class NT >
  CORE_INLINE
  bool unitP(const Polynomial <NT>& p){			// =Unit Poly?
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
	if (coeff[i] != 0) return coeff[i];
     // This ought to be an error (user should check this :
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
	assert(i <= degree);
	return coeff[i];
  }
// set functions
template < class NT >
  CORE_INLINE
  bool Polynomial<NT>::setCoeff(int i, const NT & cc) {
	if ((i<0) || (i > degree)) return false;
	coeff[i] = cc; return true;
  }

// IMPLEMENTATIONS ARE FOUND IN
#include "Poly.tcc"

CORE_END_NAMESPACE
#endif
