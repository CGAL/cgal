/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2003 Exact Computation Project

   File: Sturm.h

   Description: 
	The templated class Sturm implements Sturm sequences.
	Basic capabilities include:
		counting number of roots in an interval, 
		isolating all roots in an interval
		isolating the i-th largest (or smallest) root in interval
	It is based on the Polynomial class.

> 	BigFloat intervals are used for this (new) version.
> 	It is very important that the BigFloats used in these intervals
> 	have no error at the beginning, and this is maintained
> 	by refinement.  Note that if x, y are error-free BigFloats,
> 	then (x+y)/2 may not be error-free (in current implementaion.
> 	We have to call a special "exact divide by 2" method,
> 	(x+y).div2() for this purpose.
> 
>    TODO LIST and Potential Bugs:
> 	(1) Split an isolating interval to give definite sign
> 	(2) Should have a test for square-free polynomials
> 	(3) The copy constructor seems to fail sometimes?

   Author:  Chee Yap and Sylvain Pion, Vikram Sharma
   Date:    July 20, 2002

   Core Library  Version 1.5
   $Id$
 ************************************** */

#include <CORE/BigFloat.h>
#include <CORE/poly/Poly.h>

CORE_BEGIN_NAMESPACE

// ==================================================
// Sturm Class
// ==================================================

template < class NT >
class Sturm {
public:
  int len;			// length of array
  Polynomial<NT> * seq;		// array of polynomials of length "len"
  static int N_STOP_ITER;	// Stop IterE after this many iterations.  This is
  				//    initialized below, outside the Newton class 
  bool NEWTON_DIV_BY_ZERO;
  				// This is set to true when there is divide by
				//    zero in Newton iteration (at critical value)
				// User is responsible to check this and to reset.
  typedef Polynomial<NT> PolyNT;
  // CONSTRUCTORS
  Sturm() : len(0), NEWTON_DIV_BY_ZERO(false){		// null constructor
  }

  Sturm(Polynomial<NT> pp) : NEWTON_DIV_BY_ZERO(false) { // constructor from polynomial
    int ell = pp.getTrueDegree();
    seq = new Polynomial<NT> [ell+1];
    seq[0] = pp;	// Primitive part is important
    			// to speed up large polynomials!
			// However, for first 2 polymials, we
			// DO NOT take primitive part, because we
			// want to use them in Newton Iteration
    seq[1] = differentiate(pp);
    int i;
    for (i=2; i <= ell; i++) {
        seq[i] = seq[i-2];
	seq[i].pseudoRemainder(seq[i-1]);
	if (zeroP(seq[i])) break;
	seq[i].primPart();
    }
    len = i;
  }

  // copy constructor
  Sturm(const Sturm&s) : len(s.len), NEWTON_DIV_BY_ZERO(s.NEWTON_DIV_BY_ZERO) {
    seq = new Polynomial<NT> [len];
    for (int i=0; i<len; i++)
      seq[i] = s.seq[i];
  }

  const Sturm &operator=(const Sturm &o){
    if (this == &o) return *this;
    if (len != 0) delete[] seq;
    seq=NULL;
    NEWTON_DIV_BY_ZERO=o.NEWTON_DIV_BY_ZERO;
    len= o.len;
    if (len !=0){
      seq = new Polynomial<NT> [len];
      for (int i=0; i<len; i++)
        seq[i] = o.seq[i];
    }
    return *this;
  }


  ~Sturm() {
    if (len != 0) delete[] seq;
  }

  // METHODS

  // dump functions
  void dump(std::string msg) const {
    std::cerr << msg << std::endl;
    for (int i=0; i<len; i++)
      std::cerr << " seq[" << i << "] = " << seq[i] << std::endl;
  }
  void dump() const {
    dump("");
  }

  // signVariations(x, sx)
  //	where sx = sign of evaluating the first polynomial at x
  //	PRE-CONDITION: sx != 0
  int signVariations(const BigFloat & x, int sx) const
  {
    assert(sx != 0);
    int cnt = 0;
    int last_sign = sx;
    for (int i=0; i<len; i++)
    {
      int sgn = seq[i].eval(x).sign();
      if (sgn*last_sign < 0)
        cnt++;
      if (sgn != 0)
        last_sign = sgn;
    }
    return cnt;
  }

  // signVariations(x)
  //	--the first polynomial eval is not yet done
  int signVariations(const BigFloat & x) const
  {
	int signx = seq[0].eval(x).sign();
	if (signx == 0)
		return (-1);	// error!
	return signVariations(x, signx);
  }//signVariations(x)

  // signVariation at +Infinity
  int signVariationsAtPosInfty() const
  {
    int cnt = 0;
    int last_sign = seq[0].coeff[seq[0].getTrueDegree()].sign();
    assert(last_sign != 0);
    for (int i=1; i<len; i++)
    {
      int sgn = seq[i].coeff[seq[i].getTrueDegree()].sign();
      if (sgn*last_sign < 0)
        cnt++;
      if (sgn != 0)
        last_sign = sgn;
    }
    return cnt;
  }

  // signVariation at -Infinity
  int signVariationsAtNegInfty() const
  {
    int cnt = 0;
    int last_sign = seq[0].coeff[seq[0].getTrueDegree()].sign();
    if (seq[0].getTrueDegree() % 2 != 0) last_sign *= -1;
    assert(last_sign != 0);
    for (int i=1; i<len; i++)
    {
      int parity = (seq[i].getTrueDegree() % 2 == 0) ? 1 : -1;
      int sgn = parity * seq[i].coeff[seq[i].getTrueDegree()].sign();
      if (sgn*last_sign < 0)
        cnt++;
      if (sgn != 0)
        last_sign = sgn;
    }
    return cnt;
  }

  // numberOfRoots(x,y):
  // 	COUNT NUMBER OF ROOTS in the close interval [x,y]
  //	IMPORTANT: Must get it right even if x, y are roots
  int numberOfRoots(const BigFloat &x, const BigFloat &y) const
  {
    // assert(x <= y);   // we allow x=y (probabably never used)
    int signx = seq[0].eval(x).sign();
    int signy = seq[0].eval(y).sign();
    // easy case: THIS SHOULD BE THE OVERWHELMING MAJORITY
    if (signx != 0 && signy != 0) 
	return (signVariations(x, signx) - signVariations(y, signy));
    // harder case: THIS SHOULD BE VERY INFREQUENT
    BigFloat sep = (seq[0].sepBound())/2;
    BigFloat newx, newy;
    if (signx == 0) newx = x - sep; else newx = x;
    if (signy == 0) newy = y + sep;  else newy = y;
    return (signVariations(newx, seq[0].eval(newx).sign())
    		- signVariations(newy, seq[0].eval(newy).sign()) );
  }//numberOfRoots

  // numberOfRoots above x:
  //   NOTE: should do a special case for x=0
  ///////////////////////////////////////////
  int numberOfRootsAbove(const BigFloat &x) const
  {
    int signx = seq[0].eval(x).sign();
    if (signx != 0)
       return signVariations(x, signx) - signVariationsAtPosInfty();
    BigFloat newx = x - (seq[0].sepBound())/2;
    return signVariations(newx, seq[0].eval(newx).sign())
		- signVariationsAtPosInfty();
  }

  // numberOfRoots below x:
  //   NOTE: should do a special case for x=0
  ///////////////////////////////////////////
  int numberOfRootsBelow(const BigFloat &x) const
  {
    int signx = seq[0].eval(x).sign();
    if (signx != 0) 
       return signVariationsAtNegInfty() - signVariations(x, signx);
    BigFloat newx = x + (seq[0].sepBound())/2;
    return signVariationsAtNegInfty() 
		- signVariations(newx, seq[0].eval(newx).sign());
  }


  // isolateRoots(x, y, v) 
  ///   isolates all the roots in [x,y] and returns them in v.
  /** 	v is a list of intervals
   *    [x,y] is the initial interval to be isolated
   */
  void isolateRoots(const BigFloat &x, const BigFloat &y,
                    BFVecInterval &v) const
  {
    assert(x<=y);

    int n = numberOfRoots(x,y);
    if (n == 0)
      return;
    if (n == 1) {
       if ((x <= 0) && (y >= 0)) { 	// if 0 is inside our interval
	   if (seq[0].coeff[0] == 0) 
      		v.push_back(std::make_pair(BigFloat(0), BigFloat(0)));
	   else if (numberOfRoots(0,y) == 0)
      		v.push_back(std::make_pair(x, BigFloat(0)));
	   else
      		v.push_back(std::make_pair(BigFloat(0), y));
	}	 else
      		v.push_back(std::make_pair(x, y));
    } else {
      BigFloat mid = (x+y).div2();
      isolateRoots(x, mid, v);
      isolateRoots(mid, y, v);
    }
  }

  // isolateRoots(v)
  ///   isolates all roots and returns them in v
  /**  	v is a vector of isolated intervals
   */
  void isolateRoots(BFVecInterval &v) const
  {
    BigFloat bd = seq[0].CauchyUpperBound();
    	// Note: bd is an exact BigFloat (this is important)
    isolateRoots(-bd, bd, v);
  }

  // isolateRoot(i)
  ///   isolates the i-th smallest root
  BFInterval isolateRoot(int i) const
  {
    if (i == 0) return mainRoot();
    BigFloat bd = seq[0].CauchyUpperBound();
    return isolateRoot(i, -bd, bd);
  }

   // isolateRoot(i, x, y)
  ///	isolates the i-th smallest root in [x,y]
  /**   If i is negative, then we want the i-th largest root in [x,y]
   *    We assume i is not zero.
   */
  BFInterval isolateRoot(int i, BigFloat x, BigFloat y) const
  {
    int n = numberOfRoots(x,y);
    if (i < 0) {
      i += n+1;
      if (i <= 0) return BFInterval(1,0); // ERROR CONDITION
    }
    if (n < i) return BFInterval(1,0);  // ERROR CONDITION INDICATED
    if (n == 1) return BFInterval(x, y);
    BigFloat m = (x+y).div2();
    n = numberOfRoots(x, m);
    if (n < i)
      return isolateRoot(i-n, m, y);
    else
      return isolateRoot(i, x, m);
  }

  // same as isolateRoot(i).
  BFInterval diamond(int i) const
  {
    return isolateRoot(i);
  }

  // First root above
  BFInterval firstRootAbove(const BigFloat &e) const
  {
    return isolateRoot(1, e, seq[0].CauchyUpperBound());
  }

  // Main root (i.e., first root above 0)
  BFInterval mainRoot() const
  {
    return isolateRoot(1, 0, seq[0].CauchyUpperBound());
  }
  
  // First root below
  BFInterval firstRootBelow(const BigFloat &e) const
  {
    BigFloat bd = seq[0].CauchyUpperBound(); // bd is exact
    int n = numberOfRoots(-bd, e);
   if (n <= 0) return BFInterval(1,0);
    BigFloat bdBF = BigFloat(ceil(bd));
    if (n == 1) return BFInterval(-bdBF, e);
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

     BigFloat eps = BigFloat::exp2(-aprec);	//  eps = 2^{-aprec}
     extLong n =  width.uMSB() + (extLong)aprec; 
     

     BFInterval J = I;           // Return value is the Interval J
     BigFloat midpoint;
     while(n >= 0){
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
  BFInterval refinefirstRootAbove(const BigFloat &e, int aprec) const
  {
    BFInterval I = firstRootAbove(e);
    return refine(I,aprec);
  }

  // Refine First root below
  BFInterval refinefirstRootBelow(const BigFloat &e, int aprec) const
  {
    BFInterval I = firstRootBelow(e);
    return refine(I,aprec);
  }

  // refineAllRoots(v, aprec)
  //     will modify v so that v is a list of isolating intervals for
  //     the roots of the polynomial in *this.  The size of these intervals
  //     are at most 2^{-aprec}. 
  // If v is non-null, we assume it is a list of initial isolating intervals.
  // If v is null, we will first call isolateRoots(v) to set this up.
  void refineAllRoots( BFVecInterval &v, int aprec){

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

  // This is the experimental version of "refineAllRoots"
  // 	based on Newton iteration
  //
  void newtonRefineAllRoots( BFVecInterval &v, int aprec){

    BFVecInterval v1;
    BFInterval  J;

    if (v.empty()) 
	    isolateRoots(v);
    for (BFVecInterval::const_iterator it = v.begin();
         it != v.end(); ++it) {        // Iterate through all the intervals
        //refine them to the given precision aprec

      J = newtonRefine(*it, aprec);
      if (NEWTON_DIV_BY_ZERO) {
	      J.first = 1; J.second = 0;   // indicating divide by zero
      }
      v1.push_back(std::make_pair(J.first, J.second));
    }
    v.swap(v1);
  }//End of newtonRefineAllRoots


  /////////////////////////////////////////////////////////////////
  //  DIAGNOSTIC TOOLS 
  /////////////////////////////////////////////////////////////////
  // Polynomial tester:	P is polynomial to be tested
  // 			prec is the bit precision for root isolation
  // 			n is the number of roots predicted
  
  static void testSturm(PolyNT &P, int prec, int n = -1)
  {
    Sturm<NT> SS (P);
    BFVecInterval v;
    SS.refineAllRoots(v, prec);
    std::cout << "   Number of roots is " << v.size();
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
  // 	will run the Newton-Sturm refinement to isolate the roots of Poly
  // 		until absolute precision aprec.
  //	n is the predicated number of roots
  //		(will print an error message if n is wrong)
  static void testNewtonSturm(PolyNT &P, int prec, int n = -1)
  {
    Sturm<NT> SS (P);
    BFVecInterval v;
    SS.newtonRefineAllRoots(v, prec);
    std::cout << "   Number of roots is " << v.size();
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
  }// testNewtonSturm

  // Newton(prec, bf):
  // 	iterate until epsilon step size (i.e., epsilon < 2^{-prec})
  //	starting from initial bf
  //
  BigFloat newton(long prec, const BigFloat& bf) const {
	  // usually, prec is positive
    int count = 0;
    BigFloat val = bf;
    extLong uMSB;
    BigFloat del;
    // newton iteration
    do { 
      del = seq[0].eval(val)/seq[1].eval(val);
      del.makeExact();
      uMSB = del.uMSB();
      val -= del;       
      val.makeExact();
      count ++;
    } while ((uMSB >= -prec) && (count < N_STOP_ITER)) ;
			// N_STOP_ITER to prevent infinite loops
    if (count == N_STOP_ITER) {
	    core_error("Newton Iteration reach limit",
			    __FILE__, __LINE__, true);
    }
    return val;
  }

  // v = newtonIterN(n, bf, del)
  //    v is the root after n iterations of Newton
  // 		starting from initial value of bf
  // 	Return by reference, "del" (difference between returned val and value
  // 		in the previous Newton iteration)
  BigFloat newtonIterN(long n, const BigFloat& bf, BigFloat& del) {
    int count = 0;
    BigFloat val = bf;

    extLong uMSB = val.uMSB();
    // newton iteration
    for (int i=0; i<n; i++) { 
      BigFloat ff = seq[1].eval(val).makeExact();
      if (ff == 0) {
	 NEWTON_DIV_BY_ZERO = true;
	 del = 0;
	 core_error("Zero divisor in Newton Iteration", __FILE__, __LINE__, false);
	 return 0;
      }
      BigFloat f = seq[0].eval(val).makeExact();
      if (f == 0) {
		 NEWTON_DIV_BY_ZERO = false;
		 del = 0;    // indicates that we have reached the exact root
		 return val; // val is the exact root, before the last iteration
      }
      del = (f/ff).makeExact();
      uMSB = del.uMSB();
      val -= del;
      val.makeExact();
      count ++;
    }
    return val;
  }


  // v = newtonIterE(prec, bf, del)
  //
  // 	return the value v which is obtained by Newton iteration until del.uMSB < -prec
  // 		starting from initial value of bf
  // 	Return by reference "del" (difference between returned val and value
  // 		in the previous Newton iteration)

  BigFloat newtonIterE(int prec, const BigFloat& bf, BigFloat& del) const{
	// usually, prec is positive
    int count = 0;
    BigFloat val = bf;

    do { 
      val = newtonIterN(1, val, del);
      count++;
    } while ((del != 0) && ((del.uMSB() >= -prec) && (count < N_STOP_ITER))) ;
				// N_STOP_ITER to prevent infinite loops
    return val;
  }

  
  // Smale's bound
  //   smalesBound(p, q) where q is the derivative of p
  //   returns Smale's original bound for convergence in range space
  
  BigFloat smaleBound(const Polynomial<NT> & p, const Polynomial<NT> & q) const{
    int deg = p.getTrueDegree();
    BigFloat rho = 1/(1 + pow(q.length(), deg)
		    	  *pow(p.length() + 1, deg-1) );
    return rho/13;
  }

  // Yap's bound
  //   yapsBound(p) returns the bound on size of isolating interval
  //   which will guarantee being in Newton zone
  
  BigFloat yapsBound(const Polynomial<NT> & p) const{
    int deg = p.getTrueDegree();
    return  1/(1 + pow(BigFloat(deg), 3*deg+9)
	    	       *pow(BigFloat(2+p.height()),6*deg));
  }

  // newtonRefine(I, a) assumes I is an isolating interval for a root x^*,
  //     and will return an approximate root x which is guaranteed
  //     to be in the Newton basin of the root, and such that
  //	|x-x^*| < 2^{-a}.

  BFInterval newtonRefine(const BFInterval I, int aprec) {
    BFInterval J = I;
    int leftSign = seq[0].eval(I.first).sign();
    if (leftSign == 0) { J.first = J.second = I.first; return J;}

    int rightSign = seq[0].eval(I.second).sign();
    if (rightSign == 0) { J.first = J.second = I.second; return J;}

    assert( leftSign * rightSign < 0 );

    int steps  = 1; //The number of times Newton is called without checking 
    	// whether the iterations are still in the interval
    BigFloat x;
    BigFloat smale = smaleBound(seq[0], seq[1]);
    BigFloat yap = yapsBound(seq[0]);
    int N = steps;

    BigFloat oldval = 0.0, temp;
    BigFloat del;

    x = (J.second + J.first).div2();
    // x.makeExact();
    while ( //seq[0].eval(x).abs() > smale && (J.second - J.first) > yap && 
	   del.uMSB() >= - aprec ){
	oldval = x;

	x = newtonIterN(N, x, del);
	if (del == 0) {  // either reached exact root or Divide-by-zero.  It is
			// the user's responsibility to check which is the case.
		J.first = x; J.second = x; 
		return J;	
	}
	del = core_abs(del);
	if  ( 
	     (J.first <= x && x <= J.second)
		        &&
	     ((J.first < x - del) || (J.second > x+ del)))
	{		// we can get a better root:
	  if (x-del > J.first) {
		  int  lSign = seq[0].eval(x - del).sign();
		  if (lSign == 0) {
			  J.first = J.second = x-del; return J;
		  }
		  if (lSign == leftSign)
			  J.first = x - del;
		  else {
			  J.second = x - del;
			  x = (J.second + J.first).div2();
		  }
	  }
	  if (x+del < J.second) {
		  int  rSign = seq[0].eval(x + del).sign();
		  if (rSign == 0) {
			  J.first = J.second = x+del; return J;
		  }
		  if (rSign == rightSign)
			  J.second = x+del;
		  else {
		 	J.first = x + del;
			x = (J.second + J.first).div2();
		  }
	  }

	  N ++;		// be more confident or aggressive

	} else {
          N --;

	  if (N == 0) {   // do bisection in this case
	    x = (J.second + J.first).div2();  // x.makeExact();
       		if(seq[0].eval(J.first) * seq[0].eval(x) < 0){
	 		J.second = x;  
       		}else if(seq[0].eval(J.second) * seq[0].eval(x) < 0){
	 		J.first = x;
       		}else{
	 		J.first = x;
	 		J.second = x;
			return J;
       		}
       		x = (J.second + J.first).div2();  // x.makeExact();
				//Reset x to be the midpoint of the interval
       		N = steps;  	// restart

	  } else {
	    // Reset x to the oldval before doing Newton in this iteration
	    	x = oldval;
          }
        }
    }//end of outer while

    /* For future work. We should exploit the fact that we have reached
       the bounds. Could possible speed up things for high precision 
    if(del.uMSB() >= -aprec){// If any of the bounds are satisfied 
                             //before reaching the required precision
      x = newtonIterE(aprec, x, del);
      J.first = x; 
      J.second = J.first;
    }
    */

    return(J);

  }//End of newton refine

};// Sturm class

// ==================================================
// Static initialization
// ==================================================
template <class NT>
int Sturm<NT>:: N_STOP_ITER = 10000;	// stop IterE after this many loops
  					// Reset this as needed
CORE_END_NAMESPACE

