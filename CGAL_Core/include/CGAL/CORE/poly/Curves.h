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
 * File: Curves.h
 *
 * Description: 
 * 	Two templated classes are defined here:
 *		Curve and BiPoly
 *	These classes are parametrized by the number type
 *		(called NT) which represents the
 *		domain of the coefficients of the underlying
 *		polynomials.  Standard default is NT=BigInt, but
 *		we will allow NT=int, NT=BigRat, NT=BigFloat, NT=Expr.
 *	BiPoly represents the class of bivariate polynomials,
 *		i.e.,  each BiPoly object is an element of NT[X,Y].
 *		We store each BiPoly as a list of polynomials in X.
 *	Curve represents the class of plane curves whose equation
 *		is A(X,Y)=0, for some BiPoly A(X,Y).
 *	Features:
 *		--Constructor from strings such as
 *			"3 x^2 + 7 xy^2 - 4 x + 13".
 *		--Basic plot functions
 *
 *	To Do:
 *	  --Dump should produce human readable strings like
 *	  	"3 x^2 + 7 xy^2 - 4 x + 13".
 *	  --String constructor generalizations:
 *	  	(1) allow one "=" sign (e.g., "3 x^2 = y^2 - xy")(DONE)
 *		(2) allow general parenthesis
 *		(3) allow X and Y (DONE)
 *	  --We should be able to read/write
 *	  	curve definitions from/to files
 *	  --Plot should be more efficient (use previous roots
 *	  	to help find the next roots, there should be
 *	  	a "plot structure" that is persistent)
 *	  --Plot should refine in both x- and y-increments.
 *	  --Plot should have some option to show the
 *	  	x- and y-axes, and to label some points.
 *	  --verticalIntersect(...) should be implemented using
 *	        Polynomial<BigFloat>, not Polynomial<Expr> for efficiency
 *	  --the plot parameters (eps,xmin,xmax,ymin,ymax) must be
 *	        made part of the Curve class (static members).
 *	        Incorporate the "setParams" method into class.
 *
 *  Author:  Vikram Sharma and Chee Yap
 *  Date:    April 12, 2004
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/


#ifndef CORE_CURVES_H
#define CORE_CURVES_H

#include <fstream>
#include <list>
#include "CGAL/CORE/poly/Poly.h"

namespace CORE { 

// ==================================================
// Curve Class
// ==================================================

//typedef BigInt NT;
//typedef Expr   NT;
//typedef Polynomial<NT>        PolyNT;
//typedef std::vector<Expr>	VecExpr;
//typedef std::vector<BigInt>	VecBigInt;
//typedef std::vector<NT>       VecNT;
//typedef std::vector<Polynomial<NT> >	VecPoly;

template <class NT>
class Monomial{
  //Helper class to store the coefficients for given x-deg and y-deg 
  //Used by string input routines
 public:
  NT coeff;
  int xdeg;
  int ydeg;

  Monomial(){
  }

  Monomial(NT& cf, int& dx, int& dy){
    coeff = cf;
    xdeg = dx;
    ydeg = dy;
  }

  void dump(){
    std::cout << coeff << "X^" << xdeg << " Y^" << ydeg;
  }
};


//Class of Bivariate polynomials
//	Viewed as a polynomial in Y with
//	coefficients which are polynomials in X
template <class NT>
class BiPoly{
 private:
  //The following are used in the constructor from strings.
  //For more details see the related constructor.
  void constructFromString(string& s, char myX='x', char myY='y');
  void constructX(int n, BiPoly<NT>& P);
  void constructY(int n, BiPoly<NT>& P);
  int getnumber(const char* c, int i, unsigned int len, BiPoly<NT> & P);
  bool isint(char c);
  int getint(const char* c, int i, unsigned int len, int & n);
  int matchparen(const char* cstr, int start);
  int getbasicterm(string s, BiPoly<NT> & P);
  int getterm(string s, BiPoly<NT> & P);

 public:
  int ydeg; //Y-degree of the polynomial
  std::vector<Polynomial<NT> > coeffX; //vector of (1+ydeg) polynomials in X
	// If ydeg = d, then the polynomial is F(X,Y) =
        //   (Y^d * coeffX[d]) + (Y^{d-1} * coeffX[d-1]) +...+ (coeffX[0]).

  ////////////////////////////////////////////////////////
  //Constructors
  ////////////////////////////////////////////////////////

  //BiPoly()
  BiPoly(); // zero polynomial

  //BiPoly(n)
  BiPoly(int n);// creates a BiPoly with nominal y-degree equal to n.

  //BiPoly(vp)
  BiPoly(std::vector<Polynomial<NT> > vp); // From vector of Polynomials

  //BiPoly(p, flag):
  //	if true, it converts polynomial p(X) into P(Y)
  // 	if false, it creates the polynomial Y-p(X)
  BiPoly(Polynomial<NT> p, bool flag=false);
  
  //BiPoly(deg, d[], C[]):
  //	Takes in a list of list of coefficients.
  //	Each cofficient list represents a polynomial in X
  //
  //  deg - ydeg of the bipoly
  //  d[] - array containing the degrees of each coefficient (i.e., X poly)
  //  C[] - list of coefficients, we use array d to select the
  //      coefficients
  BiPoly(int deg, int *d, NT *C);

  //BiPoly(String s, char myX, char myY)
  //  myX and myY are names of the two variables.
  //  Default values of myX and myY are 'x' and 'y'.
  //  The string s has the form "3 x^2 + 7 xy^2 - 4 x + 13"
  //
  //  For now, we assume no parentheses, * or =.
  
  BiPoly(const string& s, char myX='x', char myY='y');
  BiPoly(const char* s, char myX='x', char myY='y');

  // copy constructor
  BiPoly(const BiPoly<NT>&);

  //Destructor
  ~BiPoly();
  //Destructor helper
  void deleteCoeffX();


  ////////////////////////////////////////////////////////
  // METHODS
  ////////////////////////////////////////////////////////
  
  // filedump (msg, ofs, com, com2)
  // 	where msg, com, com2 are strings.
  // 	msg is an message and com, com2 are the strings
  // 	preceding each output line
  // 	(e.g., msg="BiVariate Polynomial"  and com=com2="# ")
  // This is called by the other dump functions
  void dump(std::ostream & os, std::string msg = "",
      std::string com="# ", std::string com2 = "# ") const;
  // dump(ofs, msg, com) -- dump to file
  //void dump(std::ofstream & ofs, std::string msg,
  //    std::string com="# ", std::string com2="# ") const;

  // dump(msg, com) -- dump to std output
  void dump(std::string msg="", std::string com="",
      std::string com2="") const;

  /*Cannot work with these two functions right now.
    BiPoly as per now can only handle BigInt and int since
    Expr cannot be handled by Polynomial class.*/
  
  // yPolynomial(x) 
  //   returns the polynomial (in Y) when we substitute X=x
  
  /* BiPoly<NT> yPolynomial(const Expr & x) {

    VecExpr vE;

    for (int i=0; i<= ydeg; i++) {
      vE.push_back(coeffX[i].eval(x));
    }
    
    return BiPoly<NT>(vE);
  }//yPolynomial
  */

  Polynomial<NT> yPolynomial(const NT & x);

  // Expr version of yPoly (temporary hack)
  Polynomial<Expr> yExprPolynomial(const Expr & x);

  // BF version of yPoly (temporary hack)
  Polynomial<BigFloat> yBFPolynomial(const BigFloat & x);

  // xPolynomial(y) 
  //   returns the polynomial (in X) when we substitute Y=y
  //   
  //   N.B. May need the
  //   		Polynomial<Expr> xExprPolynomial(Expr y)
  //   version too...
  //
  Polynomial<NT> xPolynomial(const NT & y) ;
  
  // getYdegree()
  int getYdegree() const;
  
  // getXdegree()
  int getXdegree();

  // getTrueYdegree
  int getTrueYdegree();

  //eval(x,y)
  Expr eval(Expr x, Expr y);//Evaluate the polynomial at (x,y)

  ////////////////////////////////////////////////////////
  // Polynomial arithmetic (these are all self-modifying)
  ////////////////////////////////////////////////////////
  
  // Expands the nominal y-degree to n;
  //	Returns n if nominal y-degree is changed to n
  //	Else returns -2

  int expand(int n);

  // contract() gets rid of leading zero polynomials
  //	and returns the new (true) y-degree;
  //	It returns -2 if this is a no-op

  int contract();

  // Self-assignment
  BiPoly<NT> & operator=( const BiPoly<NT>& P);

  // Self-addition
  BiPoly<NT> & operator+=( BiPoly<NT>& P);
   
  // Self-subtraction
  BiPoly<NT> & operator-=( BiPoly<NT>& P);

  // Self-multiplication
  BiPoly<NT> & operator*=( BiPoly<NT>& P);
  
  // Multiply by a polynomial in X
  BiPoly<NT> & mulXpoly( Polynomial<NT> & p);

  //Multiply by a constant
  BiPoly<NT> & mulScalar( NT & c);

  // mulYpower: Multiply by Y^i (COULD be a divide if i<0)
  BiPoly<NT> & mulYpower(int s);
  
  // Divide by a polynomial in X.
  // We replace the coeffX[i] by the pseudoQuotient(coeffX[i], P)
  BiPoly<NT> & divXpoly( Polynomial<NT> & p);
  
  //Using the standard definition of pseudRemainder operation.
  //	--No optimization!
  BiPoly<NT>  pseudoRemainderY (BiPoly<NT> & Q);

  //Partial Differentiation
  //Partial Differentiation wrt Y
  BiPoly<NT> & differentiateY();

  BiPoly<NT> & differentiateX();
  BiPoly<NT> & differentiateXY(int m, int n);//m times wrt X and n times wrt Y

  //Represents the bivariate polynomial in (R[X])[Y] as a member
  //of (R[Y])[X].
  //But since our polynomials in X can only have NT coeff's thus
  // to represent the above polynomial we switch X and Y once
  // the conversion has been done.
  //NOTE: This is different from replacing X by Y which was
  //      done in the case of the constructor from a polynomial in X
  //Need to calculate resultant wrt X.
  BiPoly<NT> & convertXpoly();

  //Set Coeffecient to the polynomial passed as a parameter
  bool setCoeff(int i, Polynomial<NT> p);

  void reverse();
  Polynomial<NT> replaceYwithX();

  //Binary-power operator
  BiPoly<NT>& pow(unsigned int n);

  //Returns a Bipoly corresponding to s, which is supposed to
  //contain as place-holders the chars 'x' and 'y'.
  BiPoly<NT> getbipoly(string s);
};//BiPoly Class

  ////////////////////////////////////////////////////////
  // Helper Functions
  ////////////////////////////////////////////////////////
//Experimental version of constructor from strings containing general 
//parentheses


// zeroPinY(P)
//	checks whether a Bi-polynomial is a zero Polynomial
template <class NT>
bool zeroPinY(BiPoly<NT> & P);

// gcd(P,Q)
//   This gcd is based upon the subresultant PRS to avoid
//   exponential coeffecient growth and gcd computations, both of which 
//   are expensive since the coefficients are polynomials

template <class NT>
BiPoly<NT> gcd( BiPoly<NT>& P ,BiPoly<NT>& Q);

// resY(P,Q):
//      Resultant of Bi-Polys P and Q w.r.t. Y.
//      So the resultant is a polynomial in X
template <class NT>
Polynomial<NT>  resY( BiPoly<NT>& P ,BiPoly<NT>& Q);

// resX(P,Q):
//      Resultant of Bi-Polys P and Q w.r.t. X.
//      So the resultant is a polynomial in Y
//	We first convert P, Q to polynomials in X. Then 
// 	call resY and then turn it back into a polynomial in Y
//	QUESTION: is this last switch really necessary???
template <class NT>
BiPoly<NT>  resX( BiPoly<NT>& P ,BiPoly<NT>& Q);

//Equality operator for BiPoly
template <class NT>
bool operator==(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Addition operator for BiPoly
template <class NT>
 BiPoly<NT> operator+(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Subtraction operator for BiPoly
template <class NT>
 BiPoly<NT> operator-(const BiPoly<NT>& P, const BiPoly<NT>& Q);

//Multiplication operator for BiPoly
template <class NT>
 BiPoly<NT> operator*(const BiPoly<NT>& P, const BiPoly<NT>& Q);


  ////////////////////////////////////////////////////////
  //Curve Class
  //  	extends BiPoly Class
  ////////////////////////////////////////////////////////

template < class NT >
class Curve : public BiPoly<NT> {
public:
  // Colors for plotting curves

  static const int NumColors=7;
  static double red_comp(int i){
  	static const double RED_COMP[] = {0.9, 0.8, 0.7, 0.6, 0.8, 0.8, 0.7};
	return RED_COMP[i % NumColors];
  }
  static double green_comp(int i){
  	static const double GREEN_COMP[] = {0.5, 0.9, 0.3, 0.9, 0.7, 0.55, 0.95};
	return GREEN_COMP[i % NumColors];
  }
  static double blue_comp(int i){
  	static const double BLUE_COMP[] = {0.8, 0.3, 0.8, 0.5, 0.4, 0.85, 0.35};
	return BLUE_COMP[i % NumColors];
  }

  Curve(); // zero polynomial
  
  //Curve(vp):
  //    construct from a vector of polynomials
  Curve(std::vector<Polynomial<NT> > vp);
  //	  : BiPoly<NT>(vp){
  //}
  
  //Curve(p):
  //	Converts a polynomial p(X) to a BiPoly in one of two ways:
  // 	    (1) if flag is false, the result is Y-p(X) 
  // 	    (2) if flag is true, the result is p(Y) 
  //    The default is (1) because we usually want to plot the
  //        graph of the polynomial p(X)
  Curve(Polynomial<NT> p, bool flag=false);
  //	  : BiPoly<NT>(p, flag){
  //}

  //Curve(deg, d[], C[]):
  //	Takes in a list of list of coefficients.
  //	Each cofficient list represents a polynomial in X
  //
  //  deg - ydeg of the bipoly
  //  d[] - array containing the degrees of each coefficient (i.e., X poly)
  //  C[] - list of coefficients, we use array d to select the
  //      coefficients
  Curve(int deg, int *d, NT *C);
  //	  : BiPoly<NT>(deg, d, C){
  //}

  Curve(const BiPoly<NT> &P);
  //	  : BiPoly<NT>(P){
  //}

  //Curve(n) -- the nominal y-degree is n
  Curve(int n);

  //Creates a curve from a string (no parentheses, no *, no =)
  Curve(const string & s, char myX='x', char myY='y');
  Curve(const char* s, char myX='x', char myY='y');

  /////////////////////////////////////////////////////////////////////////
  // verticalIntersections(x, vecI, aprec=0):
  //    The list vecI is passed an isolating intervals for y's such that (x,y)
  //    lies on the curve.
  //    If aprec is non-zero (!), the intervals have with < 2^{-aprec}.
  //    Return is -2 if curve equation does not depend on Y
  //    	-1 if infinitely roots at x,
  //    	0 if no roots at x
  //    	1 otherwise

  int verticalIntersections(const BigFloat & x, BFVecInterval & vI,
			    int aprec=0);
  
  // TO DO: 
  // 		horizontalIntersections(...)
  
  /////////////////////////////////////////////////////////////////////////
  // plot(eps, x1, y1, x2, y2)
  //
  // 	All parameters have defaults
  //
  //    Gives the points on the curve at resolution "eps".  Currently,
  //    eps is viewed as delta-x step size (but it could change).
  //    The display is done in the rectangale 
  //    defined by [(x1, y1), (x2, y2)].
  //    The output is written into a file in the format specified
  //    by our drawcurve function (see COREPATH/ext/graphics).
  //
  //    Heuristic: the open polygonal lines end when number of roots
  //    changes...
  //
  int  plot( BigFloat eps=0.1, BigFloat x1=-1.0,
	     BigFloat y1=-1.0, BigFloat x2=1.0, BigFloat y2=1.0, int fileNo=1);

// selfIntersections():
//   this should be another member function that lists
//   all the self-intersections of a curve
//
//  template <class NT>
//  void selfIntersections(BFVecInterval &vI){
//  ...
//  }

};// Curve class


  ////////////////////////////////////////////////////////
  // Curve helper functions
  ////////////////////////////////////////////////////////


//Xintersections(C, D, vI):
//  returns the list vI of x-ccordinates of possible intersection points.
//  Assumes that C & D are quasi-monic.(or generally aligned)
template <class NT>
void  Xintersections( Curve<NT>& P ,Curve<NT>& Q, BFVecInterval &vI);

//Yintersections(C, D, vI):
//	similar to Xintersections
template <class NT>
void  Yintersections( Curve<NT>& P ,Curve<NT>& Q, BFVecInterval &vI);

// Display Intervals
template <class NT>
void showIntervals(char* s, BFVecInterval &vI);

// Set Display Parameters
// ...

////////////////////////////////////////////////////////
// IMPLEMENTATIONS ARE FOUND IN Curves.tcc
////////////////////////////////////////////////////////
#include <CGAL/CORE/poly/Curves.tcc>


} //namespace CORE
#endif
/*************************************************************************** */
// END
/*************************************************************************** */
