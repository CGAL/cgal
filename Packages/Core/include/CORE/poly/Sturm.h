/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/); you may
 * redistribute it under the terms of the Q Public License version 1.0.
 * See the file LICENSE.QPL distributed with CORE.
 *
 * Licensees holding a valid commercial license may use this file in
 * accordance with the commercial license agreement provided with the
 * software.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 *  File: Sturm.h
 * 
 *  Description: 
 *  The templated class Sturm implements Sturm sequences.
 *  Basic capabilities include:
 *     counting number of roots in an interval, 
 *     isolating all roots in an interval
 *     isolating the i-th largest (or smallest) root in interval
 *  It is based on the Polynomial class.
 * 
 *   BigFloat intervals are used for this (new) version.
 *   It is very important that the BigFloats used in these intervals
 *   have no error at the beginning, and this is maintained
 *   by refinement.  Note that if x, y are error-free BigFloats,
 *   then (x+y)/2 may not be error-free (in current implementaion.
 *   We have to call a special "exact divide by 2" method,
 *   (x+y).div2() for this purpose.
 *
 *   CONVENTION: an interval defined by a pair of BigFloats x, y
 *   has this interpretation:
 *       (1) if x>y,  it represents an invalid interval.
 *       (2) if x=y,  it represents a unique point x.
 *       (3) if x<y,  it represents the open interval (x,y).
 *           In this case, we always may sure that x, y are not zeros.
 * 
 *   TODO LIST and Potential Bugs:
 *   (1) Split an isolating interval to give definite sign (done)
 *   (2) Should have a test for square-free polynomials (done)
 * 
 *  Author:  Chee Yap and Sylvain Pion, Vikram Sharma
 *  Date:    July 20, 2002
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Source$
 * $Revision$ $Date$
 ***************************************************************************/


#ifndef CORE_STURM_H
#define CORE_STURM_H

#include "CORE/BigFloat.h"
#include "CORE/poly/Poly.h"

CORE_BEGIN_NAMESPACE

// ==================================================
// Sturm Class
// ==================================================

template < class NT >
class Sturm {
public:
  int len;      // len is 1 less than number of non-zero entries in array seq.
  		//     I.e., len + 1 = length of the Sturm Sequence
                // N.B. When len = -1 or len = 0 are special,
                //     the array seq is not used!
                //     Hence, one must test these special cases
  Polynomial<NT> * seq;      // array of polynomials of length "len+1"
  Polynomial<NT> g;//GCD of input polynomial P and it's derivative P'
  NT cont;//Content of the square-free part of input polynomial P
  //Thus P = g * cont * seq[0]
  static int N_STOP_ITER;    // Stop IterE after this many iterations. This
  // is initialized below, outside the Newton class
  bool NEWTON_DIV_BY_ZERO;   // This is set to true when there is divide by
  // zero in Newton iteration (at critical value)
  // User is responsible to check this and to reset.
  typedef Polynomial<NT> PolyNT;

  // ===============================================================
  // CONSTRUCTORS
  // ===============================================================
  // Null Constructor
  Sturm() : len(0), NEWTON_DIV_BY_ZERO(false) {}

  // Constructor from a polynomial
  Sturm(Polynomial<NT> pp) : NEWTON_DIV_BY_ZERO(false) {
    len = pp.getTrueDegree();
    if (len <= 0) return; // hence, seq is not defined in these cases
    seq = new Polynomial<NT> [len+1];
    seq[0] = pp;
    g = seq[0].sqFreePart();
    cont = content(seq[0]);
    seq[0].primPart();
    seq[1] = differentiate(seq[0]);
    int i;
    for (i=2; i <= len; i++) {
      seq[i] = seq[i-2];
      seq[i].negPseudoRemainder(seq[i-1]);
      if (zeroP(seq[i])){
	len = i;
        break;
      }
      seq[i].primPart(); // Primitive part is important to speed
      // up large polynomials! However, for first 2 polymials,
      // we MUST NOT take primitive part, because we
      // want to use them in Newton Iteration
    }
  }

  // Chee: 7/31/04
  // 	We need BigFloat version of Sturm(Polynomial<NT>pp) because
  // 	of curve verticalIntersection() ... .  We also introduce
  // 	various support methods in BigFloat.h (exact_div, gcd, etc).
  // Constructor from a BigFloat polynomial 
  //	Need the fake argument to avoid compiler overloading errors
  Sturm(Polynomial<BigFloat> pp, bool fake) : NEWTON_DIV_BY_ZERO(false) {
    len = pp.getTrueDegree();
    if (len <= 0) return; // hence, seq is not defined in these cases
    seq = new Polynomial<NT> [len+1];
    seq[0] = pp;
    g = seq[0].sqFreePart();
    cont = content(seq[0]);
    seq[0].primPart();
    seq[1] = differentiate(seq[0]);
    int i;
    for (i=2; i <= len; i++) {
      seq[i] = seq[i-2];
      seq[i].negPseudoRemainder(seq[i-1]);
      if (zeroP(seq[i])){
	len = i;
        break;
      }
      seq[i].primPart(); // Primitive part is important to speed
      // up large polynomials! However, for first 2 polymials,
      // we DO NOT take primitive part, because we
      // want to use them in Newton Iteration
    }
  }

  // Constructor from an array of NT's
  //   -- this code is identical to constructing from a polynomial...
  Sturm(int n, NT * c) : NEWTON_DIV_BY_ZERO(false) {
    Polynomial<NT> pp(n, c); // create the polynomial pp first and call the
    (*this) = Sturm<NT>(pp);//constructor from a polynomial
  }

  // copy constructor
  Sturm(const Sturm&s) : len(s.len), NEWTON_DIV_BY_ZERO(s.NEWTON_DIV_BY_ZERO) {
    if (len <= 0) return;
    seq = new Polynomial<NT> [len];
    for (int i=0; i<len; i++)
      seq[i] = s.seq[i];
  }

  // assignment operator
  const Sturm& operator=(const Sturm& o) {
    if (this == &o)
      return *this;
    if (len > 0)
      delete[] seq;
    NEWTON_DIV_BY_ZERO = o.NEWTON_DIV_BY_ZERO;
    len = o.len;
    if (len > 0) {
      seq = new Polynomial<NT>[len];
      for (int i=0; i<len; i++)
        seq[i] = o.seq[i];
    }
    return *this;
  }

  // destructor
  ~Sturm() {
    if (len != 0)
      delete[] seq;
  }

  // METHODS

  // dump functions
  void dump(std::string msg) const {
    std::cerr << msg << std::endl;
    if (len <= 0) std::cerr << " len = " << len << std::endl;
    else
       for (int i=0; i<len; i++)
         std::cerr << " seq[" << i << "] = " << seq[i] << std::endl;
  }
  void dump() const {
    dump("");
  }

  // signVariations(x, sx)
  //   where sx = sign of evaluating seq[0] at x
  //   PRE-CONDITION: sx != 0  and len > 0
  int signVariations(const BigFloat & x, int sx) const {
    assert((sx != 0) && (len >0));
    int cnt = 0;
    int last_sign = sx;
    for (int i=1; i<=len; i++) {// Chee (4/29/04): Bug fix,
        // should start iteration at i=1, not i=0.  Potential error
        // if seq[0].eval(x)=0 (though not in our usage).
      int sgn = sign(seq[i].eval(x));
      if (sgn*last_sign < 0) {
        cnt++;
        last_sign *= -1;
      }
    }
    return cnt;
  }

  // signVariations(x)
  //   --the first polynomial eval is not yet done
  //   --special return value of -1, indicating x is root!
  int signVariations(const BigFloat & x) const {
    if (len <= 0) return len; 
    int signx = sign(seq[0].eval(x));
    if (signx == 0)
      return (-1);    // THIS indicates that x is a root...
    		      // REMARK: in our usage, this case does not arise
    return signVariations(x, signx);
  }//signVariations(x)

  // signVariation at +Infinity
  int signVariationsAtPosInfty() const {
    if (len <= 0) return len;
    int cnt = 0;
    int last_sign = sign(seq[0].coeff[seq[0].getTrueDegree()]);
    assert(last_sign != 0);
    for (int i=1; i<len; i++) {
      int sgn = sign(seq[i].coeff[seq[i].getTrueDegree()]);
      if (sgn*last_sign < 0)
        cnt++;
      if (sgn != 0)
        last_sign = sgn;
    }
    return cnt;
  }

  // signVariation at -Infinity
  int signVariationsAtNegInfty() const {
    if (len <= 0) return len;
    int cnt = 0;
    int last_sign = sign(seq[0].coeff[seq[0].getTrueDegree()]);
    if (seq[0].getTrueDegree() % 2 != 0)
      last_sign *= -1;
    assert(last_sign != 0);
    for (int i=1; i<len; i++) {
      int parity = (seq[i].getTrueDegree() % 2 == 0) ? 1 : -1;
      int sgn = parity * sign(seq[i].coeff[seq[i].getTrueDegree()]);
      if (sgn*last_sign < 0)
        cnt++;
      if (sgn != 0)
        last_sign = sgn;
    }
    return cnt;
  }

  // numberOfRoots(x,y):
  //   COUNT NUMBER OF ROOTS in the close interval [x,y]
  //   IMPORTANT: Must get it right even if x, y are roots
  //   Assert("x and y are exact")
  //       [If the user is unsure of this assertion, do
  //        "x.makeExact(); y.makeExact()" before calling].
  ///////////////////////////////////////////
  int numberOfRoots(const BigFloat &x, const BigFloat &y) const {
    assert(x <= y);   // we allow x=y
    if (len <= 0) return len;  // return of -1 means infinity of roots!
    int signx = sign(seq[0].eval(x));
    if (x == y) return ((signx == 0) ? 1 : 0);
    int signy = sign(seq[0].eval(y));
    // easy case: THIS SHOULD BE THE OVERWHELMING MAJORITY

    if (signx != 0 && signy != 0)
      return (signVariations(x, signx) - signVariations(y, signy));
    // harder case: THIS SHOULD BE VERY INFREQUENT
    BigFloat sep = (seq[0].sepBound()).div2();
    BigFloat newx, newy;
    if (signx == 0)
      newx = x - sep;
    else
      newx = x;
    if (signy == 0)
      newy = y + sep;
    else
      newy = y;
    return (signVariations(newx, sign(seq[0].eval(newx)))
            - signVariations(newy, sign(seq[0].eval(newy))) );
  }//numberOfRoots

  // numberOfRoots():
  //   Counts the number of real roots of a polynomial
  ///////////////////////////////////////////
  int numberOfRoots() const {
    if (len <= 0) return len;  // return of -1 means infinity of roots!
    BigFloat bd = seq[0].CauchyUpperBound();
    return numberOfRoots(-bd, bd);
  }

  // numberOfRoots above x:
  //   Default value x=0 (i.e., number of positive roots)
  //   assert(len >= 0)
  ///////////////////////////////////////////
  int numberOfRootsAbove(const BigFloat &x = 0) const {
    if (len <= 0) return len;  // return of -1 means infinity of roots!
    int signx = sign(seq[0].eval(x));
    if (signx != 0)
      return signVariations(x, signx) - signVariationsAtPosInfty();
    BigFloat newx = x - (seq[0].sepBound()).div2();
    return signVariations(newx, sign(seq[0].eval(newx)))
           - signVariationsAtPosInfty();
  }

  // numberOfRoots below x:
  //   Default value x=0 (i.e., number of negative roots)
  //   assert(len >= 0)
  ///////////////////////////////////////////
  int numberOfRootsBelow(const BigFloat &x = 0) const {
    if (len <= 0) return len;  // return of -1 means infinity of roots!
    int signx = sign(seq[0].eval(x));
    if (signx != 0)
      return signVariationsAtNegInfty() - signVariations(x, signx);
    BigFloat newx = x + (seq[0].sepBound()).div2();
    return signVariationsAtNegInfty()
           - signVariations(newx, sign(seq[0].eval(newx)));
  }


  /// isolateRoots(x, y, v)
  ///             Assertion(x, y are exact BigFloats)
  ///   isolates all the roots in [x,y] and returns them in v.
  /**   v is a list of intervals
   *    [x,y] is the initial interval to be isolated
   *
   *    Properties we guarantee in the return values:
   *
   *    (0) All the intervals have exact BigFloats as endpoints
   *    (1) If 0 is a root, the corresponding isolating interval will be
   *        exact, i.e., we return [0,0].
   *    (2) If an interval is [0,x], it contains a positive root
   *    (3) If an interval is [y,0], it contains a negative root
   */
  void isolateRoots(const BigFloat &x, const BigFloat &y,
                    BFVecInterval &v) const {
    assert(x<=y);

    int n = numberOfRoots(x,y);
    if (n == 0) return;
    if (n == 1) {
      if ((x > 0) || (y < 0)) // usual case: 0 is not in interval
        v.push_back(std::make_pair(x, y));
      else { // if 0 is inside our interval (this extra
	     // service is not strictly necessary!)
        if (seq[0].coeff[0] == 0)
          v.push_back(std::make_pair(BigFloat(0), BigFloat(0)));
        else if (numberOfRoots(0,y) == 0)
          v.push_back(std::make_pair(x, BigFloat(0)));
        else
          v.push_back(std::make_pair(BigFloat(0), y));
      }
    } else { // n > 1
      BigFloat mid = (x+y).div2(); // So mid is exact.
      if (sign(seq[0].eval(mid)) != 0)  { // usual case: mid is non-root
      	isolateRoots(x, mid, v);
      	isolateRoots(mid, y, v); 
      } else { // special case: mid is a root
	BigFloat tmpEps = (seq[0].sepBound()).div2();  // this is exact!
	if(mid-tmpEps > x )//Since otherwise there are no roots in (x,mid)
	  isolateRoots(x, (mid-tmpEps).makeCeilExact(), v);
	v.push_back(std::make_pair(mid, mid));
	if(mid+tmpEps < y)//Since otherwise there are no roots in (mid,y)
	  isolateRoots((mid+tmpEps).makeFloorExact(), y, v);
      }
    }
  }//isolateRoots(x,y,v)

  // isolateRoots(v)
  ///   isolates all roots and returns them in v
  /**   v is a vector of isolated intervals
   */
  void isolateRoots(BFVecInterval &v) const {
    if (len <= 0) {
       v.clear(); return;
    }
    BigFloat bd = seq[0].CauchyUpperBound();
    // Note: bd is an exact BigFloat (this is important)
    isolateRoots(-bd, bd, v);
  }

  // isolateRoot(i)
  ///   Isolates the i-th smallest root 
  ///         If i<0, isolate the (-i)-th largest root
  ///   Defaults to i=0 (i.e., the smallest positive root a.k.a. main root)
  BFInterval isolateRoot(int i = 0) const {
    if (len <= 0) 
       return BFInterval(1,0);   // ERROR CONDITION
    if (i == 0)
      return mainRoot();
    BigFloat bd = seq[0].CauchyUpperBound();
    return isolateRoot(i, -bd, bd);
  }

  // isolateRoot(i, x, y)
  ///   isolates the i-th smallest root in [x,y]
  /**   If i is negative, then we want the i-th largest root in [x,y]
   *    We assume i is not zero.
   */
  BFInterval isolateRoot(int i, BigFloat x, BigFloat y) const {
    int n = numberOfRoots(x,y);
    if (i < 0) {//then we want the n-i+1 root
      i += n+1;
      if (i <= 0)
        return BFInterval(1,0); // ERROR CONDITION
    }
    if (n < i)
      return BFInterval(1,0);  // ERROR CONDITION INDICATED
    //Now 0< i <= n
    if (n == 1) {
      if ((x>0) || (y<0)) return BFInterval(x, y);
      if (seq[0].coeff[0] == NT(0)) return BFInterval(0,0);
      if (numberOfRoots(0, y)==0) return BFInterval(x,0);
      return BFInterval(0,y);
    }
    BigFloat m = (x+y).div2();
    n = numberOfRoots(x, m);
    if (n >= i)
	    return isolateRoot(i, x, m);
    // Now (n < i) but we have to be careful if m is a root
    if (sign(seq[0].eval(m)) != 0)   // usual case
      return isolateRoot(i-n, m, y);
    else
      return isolateRoot(i-n+1, m, y);
  }

  // same as isolateRoot(i).
  BFInterval diamond(int i) const {
    return isolateRoot(i);
  }

  // First root above
  BFInterval firstRootAbove(const BigFloat &e) const {
    if (len <= 0)
       return BFInterval(1,0);   // ERROR CONDITION
    return isolateRoot(1, e, seq[0].CauchyUpperBound());
  }

  // Main root (i.e., first root above 0)
  BFInterval mainRoot() const {
    if (len <= 0)
       return BFInterval(1,0);   // ERROR CONDITION
    return isolateRoot(1, 0, seq[0].CauchyUpperBound());
  }

  // First root below
  BFInterval firstRootBelow(const BigFloat &e) const {
    if (len <= 0)
       return BFInterval(1,0);   // ERROR CONDITION
    BigFloat bd = seq[0].CauchyUpperBound(); // bd is exact
    int n = numberOfRoots(-bd, e);
    if (n <= 0)
      return BFInterval(1,0);
    BigFloat bdBF = BigFloat(ceil(bd));
    if (n == 1)
      return BFInterval(-bdBF, e);
    return isolateRoot(n, -bdBF, e);
  }

  // Refine an interval I to absolute precision 2^{-aprec}
  //   THIS USES bisection only!  Use only for debugging (it is too slow)
  //
  BFInterval refine(const BFInterval& I, int aprec) const {
    // assert( There is a unique root in I )
    // We repeat binary search till the following holds
    //      width/2^n <= eps             (eps = 2^(-aprec))
    //   => log(width/eps) <= n
    //   => n = ceil(log(width/eps)) this many steps of binary search
    //   will work.
    // At each step we verify
    //           seq[0].eval(J.first) * seq[0].eval(J.second) < 0

    BigFloat width = I.second - I.first;
    if (width <= 0) return I;  // Nothing to do if the
                               //   interval I is exact or inconsistent
    BigFloat eps = BigFloat::exp2(-aprec);   //  eps = 2^{-aprec}
    extLong n =  width.uMSB() + (extLong)aprec;


    BFInterval J = I;           // Return value is the Interval J
    BigFloat midpoint;
    while(n >= 0) {
      midpoint = (J.second + J.first).div2();
      BigFloat m = seq[0].eval(midpoint);
      if (m == 0) {
        J.first = J.second = midpoint;
        return J;
      }
      if (seq[0].eval(J.first) * m < 0) {
        J.second = midpoint;
      } else {
        J.first = midpoint;
      }

      n--;
    }

    return J;
  }//End Refine

  // Refine First root above
  BFInterval refinefirstRootAbove(const BigFloat &e, int aprec) const {
    BFInterval I = firstRootAbove(e);
    return refine(I,aprec);
  }

  // Refine First root below
  BFInterval refinefirstRootBelow(const BigFloat &e, int aprec) const {
    BFInterval I = firstRootBelow(e);
    return refine(I,aprec);
  }

  // refineAllRoots(v, aprec)
  //     will modify v so that v is a list of isolating intervals for
  //     the roots of the polynomial in *this.  The size of these intervals
  //     are at most 2^{-aprec}.
  // If v is non-null, we assume it is a list of initial isolating intervals.
  // If v is null, we will first call isolateRoots(v) to set this up.
  void refineAllRoots( BFVecInterval &v, int aprec) {
    BFVecInterval v1;
    BFInterval  J;
    if (v.empty())
      isolateRoots(v);

    for (BFVecInterval::const_iterator it = v.begin();
         it != v.end(); ++it) {        // Iterate through all the intervals
      //refine them to the given precision aprec
      J = refine(BFInterval(it->first, it->second), aprec);
      v1.push_back(std::make_pair(J.first, J.second));
    }
    v.swap(v1);
  }//End of refineAllRoots

  // This is the new version of "refineAllRoots"
  //    	based on Newton iteration
  // It should be used instead of refineAllRoots!
  void newtonRefineAllRoots( BFVecInterval &v, int aprec) {

    BFVecInterval v1;
    BFInterval  J;

    if (v.empty())
      isolateRoots(v);
    for (BFVecInterval::iterator it = v.begin();
         it != v.end(); ++it) {        // Iterate through all the intervals
      //refine them to the given precision aprec
      J = newtonRefine(*it, aprec);

      if (NEWTON_DIV_BY_ZERO) {
        J.first = 1;
        J.second = 0;   // indicating divide by zero
      }
      v1.push_back(std::make_pair(J.first, J.second));
    }
    v.swap(v1);
  }//End of newtonRefineAllRoots


  /////////////////////////////////////////////////////////////////
  //  DIAGNOSTIC TOOLS
  /////////////////////////////////////////////////////////////////
  // Polynomial tester:   P is polynomial to be tested
  //          prec is the bit precision for root isolation
  //          n is the number of roots predicted

  static void testSturm(PolyNT &P, int prec, int n = -1) {
    Sturm<NT> SS (P);
    BFVecInterval v;
    SS.refineAllRoots(v, prec);
    std::cout << "   Number of roots is " << v.size() <<std::endl;
    if ((n >= 0) & (v.size() == (unsigned)n))
      std::cout << " (CORRECT!)" << std::endl;
    else
      std::cout << " (ERROR!) " << std::endl;
    int i = 0;
    for (BFVecInterval::const_iterator it = v.begin();
         it != v.end(); ++it) {
      std::cout << ++i << "th Root is in ["
      << it->first << " ; " << it->second << "]" << std::endl;
    }
  }// testSturm

  // testNewtonSturm( Poly, aprec, n)
  //   will run the Newton-Sturm refinement to isolate the roots of Poly
  //         until absolute precision aprec.
  //   n is the predicated number of roots
  //      (will print an error message if n is wrong)
  static void testNewtonSturm(PolyNT &P, int prec, int n = -1) {
    Sturm<NT> SS (P);
    BFVecInterval v;
    SS.newtonRefineAllRoots(v, prec);
    std::cout << "   Number of roots is " << v.size();
    if ((n >= 0) & (v.size() == (unsigned)n))
      std::cout << " (CORRECT!)" << std::endl;
    else
      std::cout << " (ERROR!) " << std::endl;

    int i = 0;
    for (BFVecInterval::iterator it = v.begin();
         it != v.end(); ++it) {
      std::cout << ++i << "th Root is in ["
      << it->first << " ; " << it->second << "]" << std::endl;
      if(it->second - it->first <= (1/power(BigFloat(2), prec)))
	std::cout << " (CORRECT!) Precision attained" << std::endl;
      else
	std::cout << " (ERROR!) Precision not attained" << std::endl;
    }
  }// testNewtonSturm

  // val = newtonIterN(n, bf, del, err)
  // 
  //    val is the root after n iterations of Newton
  //       starting from initial value of bf and is exact.
  //    Return by reference, "del" (difference between returned val and value
  //       in the previous Newton iteration)
  //
  //    Also, "err" is returned by reference and bounds the error in "del".
  //
  //    IMPORTANT: we assume that when x is an exact BigFloat,
  //    then Polynomial<NT>::eval(x) will be exact!
  //    
  BigFloat newtonIterN(long n, const BigFloat& bf, BigFloat& del,
		  unsigned long & err) {
    if (len <= 0) return bf;   // Nothing to do!  User must
                               // check this possibility!
    BigFloat val = bf;  
    // val.makeExact();    // val is exact

    // newton iteration
    for (int i=0; i<n; i++) {
      BigFloat ff = seq[1].eval(val);
      			// ff= f'(val), exact result as val is exact!
      if (ff == 0) {
        NEWTON_DIV_BY_ZERO = true;
        del = 0;
        core_error("Zero divisor in Newton Iteration",
                __FILE__, __LINE__, false);
        return 0;
      }
      BigFloat f = seq[0].eval(val);
      			// f = f(val), exact result as val is exact!
      if (f == 0) {
        NEWTON_DIV_BY_ZERO = false;
        del = 0;    // Indicates that we have reached the exact root
		    //    This is because eval(val) is exact!!!
        return val; // val is the exact root, before the last iteration
      }
      del = (f/ff); // But the accuracy of "f/ff" must be controllable
		    // by the caller...
      err = del.err();
      del.makeExact(); // makeExact() is necessary
      val -= del;
      // val.makeExact();  // -- unnecessary...
    }
    return val;
  }//newtonIterN

  // v = newtonIterE(prec, bf, del)
  //
  //    return the value v which is obtained by Newton iteration
  //    until del.uMSB < -prec, starting from initial value of bf.
  //    Returned value is an exact BigFloat.
  //    We guarantee at least one Newton step (so del is defined).
  //
  //    Return by reference "del" (difference between returned val and value
  //       in the previous Newton iteration).  This "del" is an upper bound
  //       on the last (f/f')-value in Newton iteration.
  //
  //    IN particular, if v is in the Newton zone of a root z^*, then z^* is
  //       guaranteed to lie inside [v-del, v+del].
  //
  //    Note that this is dangerous unless you know that bf is already
  //       in the Newton zone.  So we use the global N_STOP_ITER to
  //       prevent infinite loop.

  BigFloat newtonIterE(int prec, const BigFloat& bf, BigFloat& del) {
    // usually, prec is positive
    int count = N_STOP_ITER; // upper bound on number of iterations
    int stepsize = 1;
    BigFloat val = bf;
    unsigned long err = 0;

    do {
      val = newtonIterN(stepsize, val, del, err);
      count -= stepsize;
      stepsize++; // heuristic
    } while ((del != 0) && ((del.uMSB() >= -prec) && (count >0))) ;

    if (count == 0) core_error("newtonIterE: reached count=0",
		    	__FILE__, __LINE__, true);
    del = BigFloat(core_abs(del.m()), err, del.exp() );
    del.makeCeilExact();
    return val;
  }

  // A Smale bound which is an \'a posteriori condition. Applying 
  // Newton iteration to any point z satisfying this condition we are 
  // sure to converge to the nearest root in a certain interval of z.
  // The condition is for all k >= 2,
  //    | \frac{p^(k)(z)}{k!p'(z)} |^{1\(k-1)} < 1/8 * |\frac{p'(z)}{p(z)}|
  // Note: code below has been streamlined (Chee)
  /*
    bool smaleBound(const Polynomial<NT> * p, BigFloat z){
    int deg = p[0].getTrueDegree();
    BigFloat max, temp, temp1, temp2;
    temp2 = p[1].eval(z);
    temp = core_abs(temp2/p[0].eval(z))/8;
    BigInt fact_k = 2;
    for(int k = 2; k <= deg; k++){
      temp1 = core_abs(p[k].eval(z)/(fact_k*temp2)); 
      if(k-1 == 2)
	temp1 = sqrt(temp1);
      else
	temp1 = root(temp1, k-1);
      if(temp1 >= temp) return false; 
    }
    return true;
    }
   */

  //An easily computable Smale's point estimate for Newton as compared to the
  //one above. The criterion is
  //
  // ||f||_{\infty} * \frac{|f(z)|}{|f'(z)|^2} 
  //                * \frac{\phi'(|z|)^2}{\phi(|z|)}  < 0.03
  // where
  //           \phi(r) = \sum_{i=0}{m}r^i,
  //           m = deg(f)
  //
  //It is given as Theorem B in [Smale86].
  //Reference:- Chapter 8 in Complexity and Real Computation, by
  //            Blum, Cucker, Shub and Smale
  //
  //For our implementation we calculate an upper bound on
  //the second fraction in the inequality above.  For r>0,
  //
  //    \phi'(r)^2     m^2 (r^m + 1)^2
  //     ---------  <  -------------------          (1)
  //    \phi(r)        (r-1) (r^{m+1} - 1)
  //
  // Alternatively, we have
  // 
  //    \phi'(r)^2     (mr^{m+1} + 1)^2
  //     ---------  <  -------------------          (2)
  //    \phi(r)        (r-1)^3 (r^{m+1} - 1)
  //
  // The first bound is better when r > 1.
  // The second bound is better when r << 1.
  // Both bounds (1) and (2) assumes r is not equal to 1.
  // When r=1, the exact value is
  //
  //    \phi'(r)^2     m^2 (m + 1)
  //     ---------  =  -----------                  (3)
  //    \phi(r)            4
  // 
  bool smaleBoundTest(const BigFloat& z){
    assert(z.isExact());   // the bound only makes sense for exact z

#ifdef CORE_DEBUG
    std::cout <<"Computing Smale's bound = " <<  std::endl;
#endif

    if(seq[0].eval(z) == 0)// Reached the exact root.
      return true;

    BigFloat fprime = seq[1].eval(z);
    if (fprime == 0) return false;  // z is a critical value!
    BigFloat temp =        // assume eval(z) return exact values!
	    (seq[0].eval(z)/power(fprime, 2)).makeCeilExact();
    temp = core_abs(temp)*seq[0].height();  // remains exact
    //Thus, temp >=  ||f||_{\infty} |\frac{f(z)}{f'(z)^2}|

    int m = seq[0].getTrueDegree();    
    BigFloat x = core_abs(z);
    if (x==1)   // special case, using (3)
	    return (temp * BigFloat(m*m*(m+1)).div2().div2() < 0.03);

    BigFloat temp1;
    if (x>1) { // use formula (1)
      temp1 = power(m* (power(x, m)+1), 2);          // m^2*(x^m + 1)^2
      temp1 /= ((x - 1)*(power(x, m+1) - 1));        // formula (1)
    } else {  // use formula (2)
      temp1 = power(m*(power(x, m+1) +1), 2);        // (m*x^{m+1} + 1)^2
      temp1 /= (power(x - 1,3)*(power(x, m+1) -1));  // formula (2)
    }

#ifdef CORE_DEBUG
    std::cout <<"Value returned by Smale bound = " << temp * temp1.makeCeilExact() << std::endl;
#endif

    if(temp * temp1.makeCeilExact() < 0.03)          // make temp1 exact!
      return true;
    else
      return false;
  }//smaleBoundTest


  // yapsBound(p)
  // 	returns a bound on size of isolating interval of polynomial p
  // 	which is guaranteed to be in the Newton Zone.
  //    N.B. p MUST be square-free
  //
  //   Reference: Theorem 6.37, p.184 of Yap's book
  //   	   [Fundamental Problems of Algorithmic Algebra]

  BigFloat yapsBound(const Polynomial<NT> & p) const {
    int deg = p.getTrueDegree();
    return  1/(1 + pow(BigFloat(deg), 3*deg+9)
               *pow(BigFloat(2+p.height()),6*deg));
  }

  //newtonRefine(J, a) 
  //
  //    ASSERT(J is an isolating interval for some root x^*)
  //
  //    ASSERT(J.first and J.second are exact BigFloats)
  //
  //    Otherwise, the boundaries of the interval are not well defined.
  //    We will return a refined interval with exact endpoints,
  //    still called J, containing x^* and
  //
  // 			|J| < 2^{-a}.
  //
  // 	TO DO: write a version of newtonRefine(J, a, sign) where
  // 	sign=J.first.sign(), since you may already know the sign
  // 	of J.first.  This will skip the preliminary stuff in the
  // 	current version.
  //
  BFInterval newtonRefine(BFInterval &J, int aprec) {

#ifdef CORE_DEBUG_NEWTON
std::cout << "In newtonRefine, input J=" << J.first
	<< ", " << J.second << std::endl;
#endif

    if (len <= 0) return J;   // Nothing to do!  User must
                               // check this possibility!
      

    if((J.second - J.first).uMSB() < -aprec){
      return (J);
    }
    int xSign, leftSign, rightSign;

    leftSign = sign(seq[0].eval(J.first));
    if (leftSign == 0) {
      J.second = J.first;
      return J;
    }

    rightSign = sign(seq[0].eval(J.second));
    if (rightSign == 0) {
      J.first = J.second;
      return J;
    }

    assert( leftSign * rightSign < 0 );

    //N is number of times Newton is called without checking
    // whether the result is still in the interval or not
    #define NO_STEPS 1
    int N = NO_STEPS;

    BigFloat x, del, olddel, temp;
    unsigned long err;
    BigFloat yap = yapsBound(seq[0]);

    x = (J.second + J.first).div2();

    //MAIN WHILE LOOP. We always make sure that J contains the root
    while ( !smaleBoundTest(x) && 
	    (J.second - J.first) > yap &&
	   (J.second - J.first).uMSB() >= -aprec) {
      x = newtonIterN(N, x, del, err);
      if ((del == 0)&&(NEWTON_DIV_BY_ZERO == false)) {  // reached exact root!
        J.first = J.second = x;
        return J;
      }

      BigFloat left(x), right(x);
      if (del>0) {
      	left -= del; right += del;
      } else {
      	left += del; right -= del;
      }
      //left and right are exact, since x is exact.
      if ((  (J.first <= x && x <= J.second) &&
	       (J.first < left || right < J.second ))
          && (NEWTON_DIV_BY_ZERO == false)
	   ) {
                                  // we can get a better root:
        if (left > J.first) {
	  int lSign = sign(seq[0].eval(left));
          if (lSign == leftSign)  // leftSign=sign of J.first
            J.first = left;
	  else if (lSign == 0) {
            J.first = J.second = left;
            return J;
          }
        }
        if (right < J.second) {
	  int rSign = sign(seq[0].eval(right));
          if (rSign == rightSign)
            J.second = right;
	  else if (rSign == 0) {
            J.first = J.second = right;
            return J;
          }
        }

	// No, it is not necessary to update x to
	// the midpoint of the new interval J.
	// REASON: basically, it is hard to be smarter than Newton's method!
	// Newton might bring x very close to one endpoint, but it can be
	// because the root is near there!  In any case,
	// by setting x to the center of J, you only gain at most
	// one bit of accuracy, but you stand to loose an
	// arbitrary amount of bits of accuracy if you are unlucky!
	// So I will comment out the next line.  --Chee (Aug 9, 2004).
	// 
	// x = (J.second + J.first).div2();

        N ++;      // be more confident or aggressive
	           //  (perhaps we should double N)
		   //
      } else {// Either NEWTON_DIV_BY_ZERO=true
	      // Or Newton took us out of interval J: so we need to backup
	x = (J.second + J.first).div2();//Reset x to midpoint since it was the
	                                //value from a failed Newton step
	xSign = sign(seq[0].eval(x));
	if (xSign == rightSign) {
	  J.second = x;
	} else if (xSign == leftSign) {
	  J.first = x;
	} else { // xSign must be 0
	  J.first = J.second = x; return J;
	}
	x = (J.second + J.first).div2();
	// reduce value of N:
	N = core_max(N-1, NO_STEPS);   // N must be at least NO_STEPS
      }
    }//MAIN WHILE LOOP

    if((J.second - J.first).uMSB() >= -aprec){ // The interval J
                    //still hasn't reached the required precision.
                    //But we know the current value of x (call it x_0)
		    //is in the strong Newton basin of the
		    //root x^* (because it passes Smale's bound)
      //////////////////////////////////////////////////////////////////
      //Both x_0 and the root x^* are in the interval J.
      //Let NB(x^*) be the strong Newton basin of x^*.  By definition,
      //x_0 is in NB(x^*) means that:
      //
      //    x_0 is in NB(x^*) iff |x_n-x^*| \le 2^{-2^{n}+1} |x_0-x^*|
      //    
      // where x_n is the n-th iterate of Newton.  
      //    
      //  LEMMA 1: if x_0  \in NB(x^*) then 
      //               |x_0 - x^*| <= 2|del|      (*)
      //  and
      //               |x_1 - x^*| <= |del|       (**)
      //
      //  where del = -f(x_0)/f'(x_0) and x_1 = x_0 + del
      //Proof:
      //Since x_0 is in the strong Newton basin, we have
      //         |x_1-x^*| <= |x_0-x^*|/2.        (***)
      //The bound (*) is equivalent to
      //         |x_0-x^*|/2 <= |del|.
      //This is equivalent to
      //         |x_0-x^*| - |del| <= |x_0-x^*|/2,
      //which follows from
      //         |x_0-x^* + del| <= |x_0-x^*|/2,
      //which is equivalent to (***).  
      //The bound (**) follows from (*) and (***).
      //QED
      //
      //  COMMENT: the above derivation refers to the exact values,
      //  but what we compute is X_1 where X_1 is an approximation to
      //  x_1.  However, let us write X_1 = x_0 - DEL, where DEL is
      //  an approximation to del.  
      //
      //  LEMMA 2:  If |DEL| >= |del|,
      //  then (**) holds with X_1 and DEL in place of x_1 and del.
      //    

#ifdef CORE_DEBUG
      std::cout << "Inside Newton Refine: Refining Part " << std::endl;

      if((J.second - J.first) > yap)
	std::cout << "Smales Bound satisfied " << std::endl;
      else
	std::cout << "Chees Bound satisfied " << std::endl;
#endif
      xSign = sign(seq[0].eval(x));
      if(xSign == 0){
	J.first = J.second = x; 
	return J; // missing before!
      }

      //int k = clLg((-(J.second - J.first).lMSB() + aprec).asLong());
      x = newtonIterE(aprec, x, del);
      xSign = sign(seq[0].eval(x));

      if(xSign == leftSign){//Root is greater than x
	J.first = x;
	J.second = x + del;  // justified by Lemma 2 above
      }else if(xSign == rightSign){//Root is less than x
	J.first = x - del;   // justified by Lemma 2 above
	J.second = x ;
      }else{//x is the root
	J.first = J.second = x;
      }
    }



#ifdef CORE_DEBUG
    std::cout << " Returning from Newton Refine: J.first = " << J.first << " J.second = " << J.second << " aprec = " << aprec << " Sign at the interval endpoints = " << sign(seq[0].eval(J.first)) << " : " << sign(seq[0].eval(J.second)) << " Err at starting = " << J.first.err() << " Err at end = " << J.second.err() << std::endl;
#endif

    if(seq[0].eval(J.first) * seq[0].eval(J.second) > 0)
      std::cerr <<" ERROR! Root is not in the Interval " << std::endl;

#ifdef CORE_DEBUG_NEWTON
    if(J.second - J.first >  BigFloat(1).exp2(-aprec))
      std::cerr << "ERROR! Newton Refine failed to achieve the desired precision" << std::endl;
#endif

      return(J);
 }//End of newton refine

};// Sturm class

// ==================================================
// Static initialization
// ==================================================
template <class NT>
int Sturm<NT>:: N_STOP_ITER = 10000;   // stop IterE after this many loops
// Reset this as needed

// ==================================================
// Helper Functions 
// ==================================================

// isZeroIn(I):
//          returns true iff 0 is in the closed interval I
inline bool isZeroIn(BFInterval I) {
	return ((I.first <= 0.0) && (I.second >= 0.0));
}

CORE_END_NAMESPACE

#endif
