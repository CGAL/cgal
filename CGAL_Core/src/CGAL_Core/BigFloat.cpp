/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CORE (http://cs.nyu.edu/exact/core/).
 * You can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation,
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
 * File: BigFloat.cpp
 * Synopsis:
 *       BigFloat numbers with error bounds 
 *
 *       EXACTNESS PROPERTY:
 *       ==================
 *       For BigFloats that are exact (i.e., error=0),
 *       addition/subtraction and multiplication return the
 *       exact result (i.e., error=0).  We also introduce the operation
 *       div2(), which simply divides a BigFloat by 2,
 *       but this again preserves exactness.  Such exactness
 *       properties are used in our Newton iteration/Sturm Sequences.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#include <ctype.h>
#include <CGAL/CORE/BigFloat.h>
#include <CGAL/CORE/Expr.h>

namespace CORE { 


////////////////////////////////////////////////////////////
// Misc Helper Functions
////////////////////////////////////////////////////////////

BigInt FiveTo(unsigned long exp) {
  if (exp == 0)
    return BigInt(1);
  else if (exp == 1)
    return BigInt(5);
  else {
    BigInt x = FiveTo(exp / 2);

    x = x * x;

    if (exp & 1)
      x *= 5;

    return x;
  }
}

////////////////////////////////////////////////////////////
//  class BigFloat
////////////////////////////////////////////////////////////

// STATIC BIGFLOAT CONSTANTS
// ZERO
const BigFloat& BigFloat::getZero() {
  static BigFloat Zero(0);
  return Zero;
}
// ONE
const BigFloat& BigFloat::getOne() {
  static BigFloat One(1);
  return One;
}

// A special constructor for BigFloat from Expr
// -- this method is somewhat of an anomaly (we normally do not expect
//    BigFloats to know about Expr).
BigFloat::BigFloat(const Expr& E, const extLong& r, const extLong& a)
  : RCBigFloat(new BigFloatRep()) {
  *this = E.approx(r, a).BigFloatValue(); // lazy implementaion, any other way?
}

////////////////////////////////////////////////////////////
//  class BigFloatRep
////////////////////////////////////////////////////////////

BigFloatRep::BigFloatRep(double d) : m(0), err(0), exp(0) {
  if (d != 0.0) {
    int isNegative = 0;

    if (d < 0.0) {
      isNegative = 1;
      d          = - d;
    }

    int    binExp;
    double f = frexp(d, &binExp);

    exp = chunkFloor(binExp);

    long s = binExp - bits(exp);

    long   stop = 0;
    double intPart;

    // convert f into a BigInt
    while (f != 0.0 && stop < DBL_MAX_CHUNK) {
      f =   ldexp(f, (int)CHUNK_BIT);
      f =   modf(f, &intPart);
      m <<= CHUNK_BIT;
      m +=  (long)intPart;
      exp--;
      stop++;
    }
#ifdef CORE_DEBUG
    assert (s >= 0);
#endif

    if (s)
      m <<= s;
    if (isNegative)
      negate(m);
  }
}//BigFloatRep constructor

//  approximation
void BigFloatRep::trunc(const BigInt& I, const extLong& r, const extLong& a) {
  if (sign(I)) {
    long tr = chunkFloor((- r + bitLength(I)).asLong());
    long ta = chunkFloor(- a.asLong());
    long t;

    if (r.isInfty() || a.isTiny())
      t = ta;
    else if (a.isInfty())
      t = tr;
    else
      t = ta < tr ? tr : ta;

    if (t > 0) {  // BigInt remainder;
      m   = chunkShift(I, - t);
      err = 1;
      exp = t;
    } else { //  t <= 0
      m   = I;
      err = 0;
      exp = 0;
    }
  } else {//  I == 0
    m   = 0;
    err = 0;
    exp = 0;
  }
}

void BigFloatRep :: truncM(const BigFloatRep& B, const extLong& r, const extLong& a) {
  if (sign(B.m)) {
    long tr = chunkFloor((- 1 - r + bitLength(B.m)).asLong());
    long ta = chunkFloor(- 1 - a.asLong()) - B.exp;
    long t;

    if (r.isInfty() || a.isTiny())
      t = ta;
    else if (a.isInfty())
      t = tr;
    else
      t = ta < tr ? tr : ta;

    if (t >= chunkCeil(clLg(B.err))) {
      m   = chunkShift(B.m, - t);
      err = 2;
      exp = B.exp + t;
    } else //  t < chunkCeil(clLg(B.err))
      core_error(std::string("BigFloat error: truncM called with stricter")
	  + "precision than current error.", __FILE__, __LINE__, true);
  } else {//  B.m == 0
    long t = chunkFloor(- a.asLong()) - B.exp;

    if (t >= chunkCeil(clLg(B.err))) {
      m   = 0;
      err = 1;
      exp = B.exp + t;
    } else //  t < chunkCeil(clLg(B.err))
      core_error(std::string("BigFloat error: truncM called with stricter")
	  + "precision than current error.", __FILE__, __LINE__, true);
  }
}

// This is the main approximation function
// REMARK: would be useful to have a self-modifying version
// 		of this function (e.g., for Newton).
void BigFloatRep::approx(const BigFloatRep& B,
              const extLong& r, const extLong& a) {
  if (B.err) {
    if (1 + clLg(B.err) <= bitLength(B.m))
      truncM(B, r + 1, a);
    else //  1 + clLg(B.err) > lg(B.m)
      truncM(B, CORE_posInfty, a);
  } else {//  B.err == 0
    trunc(B.m, r, a - bits(B.exp));
    exp += B.exp;
  }
  // Call normalization globally     -- IP 10/9/98
  normal();
}

void BigFloatRep::div(const BigInt& N, const BigInt& D,
              const extLong& r, const extLong& a) {
  if (sign(D)) {
    if (sign(N)) {
      long tr = chunkFloor((- r + bitLength(N) - bitLength(D) - 1).asLong());
      long ta = chunkFloor(- a.asLong());

      if (r.isInfty() || a.isTiny())
        exp = ta;
      else if (a.isInfty())
        exp = tr;
      else
        exp = ta < tr ? tr : ta;

      BigInt remainder;

      // divide(chunkShift(N, - exp), D, m, remainder);
      div_rem(m, remainder, chunkShift(N, - exp), D);

      if (exp <= 0 && sign(remainder) == 0)
        err = 0;
      else
        err = 1;
    } else {//  N == 0
      m   = 0;
      err = 0;
      exp = 0;
    }
  } else //  D == 0
    core_error( "BigFloat error: zero divisor.", __FILE__, __LINE__, true);

  // Call normalization globally     -- IP 10/9/98
  normal();
}//div

//  error-normalization
void BigFloatRep::normal() {
  long le = flrLg(err);

  if (le >= CHUNK_BIT + 2) { // so we do not carry more than 16 = CHUNK_BIT + 2
	                     // bits of error
    long f = chunkFloor(--le); // f is roughly equal to floor(le/CHUNK_BIT)
    long bits_f = bits(f);   // f chunks will have bits_f many bits
#ifdef CORE_DEBUG
    assert (bits_f >= 0);
#endif

    m   >>= bits_f;  // reduce mantissa by bits_f many bits
    err >>= bits_f;  // same for err
    err +=  2;       // why 2?
    exp +=  f;       
  }
  if (err == 0)      // unlikely, if err += 2 above
    eliminateTrailingZeroes();
}

// bigNormal(err) 
//     convert a bigInt error value (=err) into an error that fits into
//     a long number.  This is done by
//     by increasing the exponent, and corresponding decrease
//     in the bit lengths of the mantissa and error.
//
void BigFloatRep::bigNormal(BigInt& bigErr) {
  long le = bitLength(bigErr);

  if (le < CHUNK_BIT + 2) {
    err = ulongValue(bigErr);
  } else {
    long f = chunkFloor(--le);
    long bits_f = bits(f);
#ifdef CORE_DEBUG
    assert(bits_f >= 0);
#endif

    m      >>= bits_f;
    bigErr >>= bits_f;
    err    = ulongValue(bigErr) + 2; // you need to add "2" because "1" comes
    		// from truncation error in the mantissa, and another
		// "1" comes from the truncation error in the bigErr.
		// (But there is danger of overflow...)
    exp    += f;
  }

  if (err == 0)
    eliminateTrailingZeroes();
}

// ARITHMETIC:
//  Addition
void BigFloatRep::add(const BigFloatRep& x, const BigFloatRep& y) {
  long expDiff = x.exp - y.exp;

  if (expDiff > 0) {//  x.exp > y.exp
    if (!x.err) {
      m   = chunkShift(x.m, expDiff) + y.m;
      err = y.err;
      exp = y.exp;
    } else {//  x.err > 0
      m   = x.m + chunkShift(y.m, - expDiff); // negative shift!
      err = x.err + 5; // To account for y.err (but why 5?)
      exp = x.exp;     // 
      // normal();
    }
  } else if (!expDiff) {//  x.exp == y.exp
    m   = x.m + y.m;
    err = x.err + y.err;
    exp = x.exp;
    // normal();
  } else {//  x.exp < y.exp
    if (!y.err) {
      m   = x.m + chunkShift(y.m, - expDiff);
      err = x.err;
      exp = x.exp;
    } else {//  y.err > 0
      m   = chunkShift(x.m, expDiff) + y.m;
      err = y.err + 5;
      exp = y.exp;
      // normal();
    }
  }
  // Call normalization globally     -- IP 10/9/98
  normal();
}

//  Subtraction
void BigFloatRep::sub(const BigFloatRep& x, const BigFloatRep& y) {
  long expDiff = x.exp - y.exp;

  if (expDiff > 0) {//  x.exp > y.exp
    if (!x.err) {
      m   = chunkShift(x.m, expDiff) - y.m;
      err = y.err;
      exp = y.exp;
    } else {//  x.err > 0
      m   = x.m - chunkShift(y.m, - expDiff);
      err = x.err + 5;
      exp = x.exp;
      // normal();
    }
  } else if (!expDiff) {
    m   = x.m - y.m;
    err = x.err + y.err;
    exp = x.exp;
    // normal();
  } else { //  x.exp < y.exp
    if (!y.err) {
      m   = x.m - chunkShift(y.m, - expDiff);
      err = x.err;
      exp = x.exp;
    } else {//  y.err > 0
      m   = chunkShift(x.m, expDiff) - y.m;
      err = y.err + 5;
      exp = y.exp;
      // normal();
    }
  }
  // Call normalization globally     -- IP 10/9/98
  normal();
}

void BigFloatRep::mul(const BigFloatRep& x, const BigFloatRep& y) {
  m = x.m * y.m;
  exp = x.exp + y.exp;
  // compute error (new code, much faster. Zilin Du, Nov 2003)
  if (x.err == 0 && y.err == 0) {
    err = 0;
    eliminateTrailingZeroes();
  } else {
    BigInt bigErr(0);
    if (y.err != 0)
      bigErr += abs(x.m)*y.err;
    if (x.err != 0)
      bigErr += abs(y.m)*x.err;
    if (x.err !=0 && y.err != 0)
      bigErr += x.err*y.err;
    bigNormal(bigErr);
  }
}
// BigFloat div2 will half the value of x, exactly with NO error
// 	REMARK: should generalize this to dividing by any power of 2
// 	We need this in our use of BigFloats to maintain isolation
// 	intervals (e.g., in Sturm sequences)	--Chee/Vikram 4/2003
//
void BigFloatRep :: div2(const BigFloatRep& x) {
  if (isEven(x.m)) {
    m = (x.m >> 1);
    exp = x.exp ;
  } else {
    m = (x.m << static_cast<unsigned long>(CHUNK_BIT-1));
    exp = x.exp -1;
  }
}

// Converts a BigFloat interval into one BigFloat with almost same error bound
// This routine ignores the errors in inputs a and b.
// But you cannot really ignore them since, they are taken into account
// when you compute "r.sub(a,b)"...
void BigFloatRep::centerize(const BigFloatRep& a, const BigFloatRep& b) {
  if ((a.m == b.m) && (a.err == b.err) && (a.exp == b.exp)) {
    m = a.m;
    err = a.err;
    exp = a.exp;
    return;
  }

  BigFloatRep r;
  r.sub(a, b);
  r.div2(r);

  //setup mantissa and exponent, but not error bits
  // But this already sets the error bits? Chee
  add(a,b);
  div2(*this);
  // error bits = ceil ( B^{-exp}*|a-b|/2 )

  // bug fixed: possible overflow on converting
  // Zilin & Vikram, 08/24/04
  // err = 1 + longValue(chunkShift(r.m, r.exp - exp));
  BigInt E = chunkShift(r.m, r.exp - exp);
  bigNormal(E);
}

// BigFloat Division, computing x/y:
//      Unlike +,-,*, this one takes a relative precision bound R
//	Note that R is only used when x and y are error-free!
//	(This remark means that we may be less efficient than we could be)
//
//  	Assert( R>0  && R< CORE_Infty )
//
void BigFloatRep :: div(const BigFloatRep& x, const BigFloatRep& y,
                        const extLong& R) {
  if (!y.isZeroIn()) { //  y.m > y.err, so we are not dividing by 0
    if (!x.err && !y.err) {
      if (R < 0 || R.isInfty()) //Oct 9, 2002: fixed major bug! [Zilin/Chee]
        div(x.m, y.m, defBFdivRelPrec, CORE_posInfty);
      else
        div(x.m, y.m, R, CORE_posInfty);
      exp += x.exp - y.exp; // chen: adjust exp.
    } else {//  x.err > 0 or y.err > 0
      BigInt bigErr, errRemainder;

      if (x.isZeroIn()) { //  x.m <= x.err
        m   = 0;
        exp = x.exp - y.exp;

        div_rem(bigErr, errRemainder, abs(x.m) + static_cast<long>(x.err),
                abs(y.m) - static_cast<long>(y.err));
      } else { //  x.m > x.err
        long lx = bitLength(x.m);
        long ly = bitLength(y.m);
        long r;

        if (!x.err) //  x.err == 0 and y.err > 0
          r = ly + 2;
        else if(!y.err)  //  x.err > 0 and y.err == 0
          r = lx + 2;
        else  //  x.err > 0 and y.err > 0
          r = lx < ly ? lx + 2: ly + 2;

        long   t = chunkFloor(- r + lx - ly - 1);
        BigInt remainder;

        div_rem(m, remainder, chunkShift(x.m, - t), y.m);
        exp = t + x.exp - y.exp;

        long delta = ((t > 0) ? 2 : 0);

        // Chen Li: 9/9/99
        // here again, it use ">>" operator with a negative
        // right operand. So the result is not well defined.
        // Erroneous code:
        //   divide(abs(remainder) + (static_cast<long>(x.err) >> bits(t))
        //                     + delta + static_cast<long>(y.err) * abs(m),
        //                     abs(y.m) - static_cast<long>(y.err),
        //                     bigErr,
        //                     errRemainder);
        // New code:
        BigInt errx_over_Bexp = x.err;
        long bits_Bexp = bits(t);
        if (bits_Bexp >= 0) {
          errx_over_Bexp >>= bits_Bexp;
        } else {
          errx_over_Bexp <<= (-bits_Bexp);
        }

        // divide(abs(remainder) + errx_over_Bexp
        //        + delta + static_cast<long>(y.err) * abs(m),
        //        abs(y.m) - static_cast<long>(y.err),
        //        bigErr,
        //        errRemainder);
        div_rem(bigErr, errRemainder,
                abs(remainder) + errx_over_Bexp + delta + static_cast<long>(y.err) * abs(m),
                abs(y.m) - static_cast<long>(y.err));
      }

      if (sign(errRemainder))
        ++bigErr;

      bigNormal(bigErr);
    }
  } else {//  y.m <= y.err
    core_error("BigFloat error: possible zero divisor.",
		    __FILE__, __LINE__, true);
  }

  // Call normalization globally     -- IP 10/9/98
  // normal(); -- Chen: after calling bigNormal, this call is redundant.
}// BigFloatRep::div

//  squareroot for BigInt argument, without initial approximation
//  sqrt(x,a) computes sqrt of x to absolute precision a.
//      -- this is where Newton is applied
//      -- this is called by BigFloatRep::sqrt(BigFloat, extLong)
void BigFloatRep::sqrt(const BigInt& x, const extLong& a) {
  sqrt(x, a, BigFloat(x, 0, 0));
} // sqrt(BigInt x, extLong a) , without initial approx

//  sqrt(x,a,A) where
//      x = bigInt whose sqrt is to be computed
//      a = absolute precision bound
//      A = initial approximation in BigFloat
//  -- this is where Newton is applied
//  -- it is called by BigFloatRep::sqrt(BigFloatRep, extLong, BigFloat)
void BigFloatRep::sqrt(const BigInt& x, const extLong& a, const BigFloat& A) {
  if (sign(x) == 0) {
    m = 0;
    err = 0;
    exp = 0;
  } else if (x == 1) {
    m = 1;
    err = 0;
    exp = 0;
  } else  {// main case
    // here is where we use the initial approximation
    m = A.m();
    err = 0;
    exp = A.exp();

    BigFloatRep q, z;
    extLong     aa;
    // need this to make sure that in case the
    // initial approximation A is less than sqrt(x)
    // then Newton iteration will still proceed at
    // least one step.
    bool firstTime = true;
    for (;;) {
      aa    = a - bits(exp);
      q.div(x, m, CORE_posInfty, aa);
      q.err = 0;
      q.exp -= exp;

      z.sub(*this, q);  // this=current approximation, so z = this - q
      /*if (sign(z.m) <= 0 || z.MSB() < - a)  // justification: see Koji's
          break;                              // thesis (p. 28) which states
                                              // that we can exit when
                                              // " (*this) <= q + 2**(-a)"
      */
      // The preceding code is replaced by what follows:
      if (z.MSB() < -a)
        break;
      if (sign(z.m) <= 0) {
        if (firstTime)
          firstTime = false;
        else
          break;
      }

      z.add(*this, q);
      // Chen Li: a bug fixed here.
      //      m   = z.m >> 1;
      //      err = 0;
      //      exp = z.exp;
      if ((z.m > 1) && isEven(z.m)) {
        m = z.m >> 1;               // exact division by 2
        err = 0;
        exp = z.exp;
      } else {                      // need to shift left before division by 2
        m = chunkShift(z.m, 1) >> 1;
        err = 0;
        exp = z.exp - 1;
      }//else
    }//for
  }//else
} // sqrt of BigInt, with initial approx

// MAIN ENTRY INTO SQRT FUNCTION (BIGFLOAT ARGUMENT, WITHOUT INITIAL APPROX)
void BigFloatRep::sqrt(const BigFloatRep& x, const extLong& a) {
  sqrt(x, a, BigFloat(x.m, 0, x.exp));
} //sqrt(BigFloat, extLong a)

// MAIN ENTRY INTO SQRT FUNCTION (BIGFLOAT ARGUMENT WITH INITIAL APPROXIMATION)
void BigFloatRep::sqrt(const BigFloatRep& x, const extLong& a, const BigFloat& A) {
  // This computes the sqrt of x to absolute precision a, starting with
  // the initial approximation A
  if (sign(x.m) >= 0) {          //  x.m >= 0
    int delta = x.exp & 1;    // delta=0 if x.exp is even, otherwise delta=1

    if (x.isZeroIn()) {         //  x.m <= x.err
      m = 0;
      if (!x.err)
        err = 0;
      else {                  //  x.err > 0
        err = (long)(std::sqrt((double)x.err));
        err++;
        err <<= 1;
        if (delta)
          err <<= HALF_CHUNK_BIT;
      }
      exp = x.exp >> 1;
      normal();
    } else {
      long aExp = A.exp() - (x.exp >> 1);
      BigFloat AA( chunkShift(A.m(), delta), 0, aExp);

      if (!x.err) {             //  x.m > x.err = 0 (ERROR FREE CASE)
        BigFloatRep z;
        extLong ppp;
        if (a.isInfty())        //Oct 9, 2002: fixed major bug! [Zilin/Chee]
          ppp = defBFsqrtAbsPrec;
        else
          ppp = a + EXTLONG_ONE;
        extLong absp  = ppp + bits(x.exp >> 1);

        z.sqrt(chunkShift(x.m, delta), absp, AA); // call sqrt(BigInt, a, AA)

        long p = (absp + bits(z.exp)).asLong();

        // Next, normalize the error:
        if (p <= 0) {
          m = z.m;
          // Chen Li: a bug fixed
          //    BigInt bigErr = 1 << (-p);
          BigInt bigErr(1);
          bigErr = bigErr << static_cast<unsigned long>(-p);
          exp = z.exp + (x.exp >> 1);
          bigNormal(bigErr);
        } else {                  //  p > 0
          m = chunkShift(z.m, chunkCeil(p));
          long r = CHUNK_BIT - 1 - (p + CHUNK_BIT - 1) % CHUNK_BIT;
#ifdef CORE_DEBUG
          assert(r >= 0);
#endif

          err = 1 >> r;
          exp = - chunkCeil(ppp.asLong());
          normal();
        }
      } else {                      //  x.m > x.err > 0 (mantissa has error)
        BigFloatRep z;
        extLong absp=-flrLg(x.err)+bitLength(x.m)-(bits(delta) >> 1)+EXTLONG_FOUR;

        z.sqrt(chunkShift(x.m, delta), absp, AA);

        long qqq = - 1 + (bitLength(x.m) >> 1) - delta * HALF_CHUNK_BIT;
        long qq  = qqq - clLg(x.err);
        long q   = qq + bits(z.exp);

        if (q <= 0) {
          m = z.m;
          long   qqqq   = - qqq - bits(z.exp);
          // Chen Li (09/08/99), a bug fixed here:
          //        BigInt bigErr = x.err << - qqqq;
          // when (-qqqq) is negative, the result is not correct.
          // how "<<" and ">>" process negative second operand is
          // not well defined. Seems it just take it as a unsigned
          // integer and extract the last few bits.
          // x.err is a long number which easily overflows.
          // From page 22 of Koji's paper, I think the exponent is
          // wrong here. So I rewrote it as:
          BigInt bigErr = x.err;
          if (qqqq >= 0) {
            bigErr <<= qqqq;
          } else {
            bigErr >>= (-qqqq);
            ++bigErr; // we need to keep its ceiling.
          }

          exp = z.exp + (x.exp >> 1);
          bigNormal(bigErr);
        } else {         //  q > 0
          m = chunkShift(z.m, chunkCeil(q));
          long r = CHUNK_BIT - 1 - (q + CHUNK_BIT - 1) % CHUNK_BIT;
#ifdef CORE_DEBUG
          assert(r >= 0);
#endif

          err = 1 >> r;
          exp = (x.exp >> 1) - chunkCeil(qq);
          normal();
        }
      }  // end of case with error in mantissa
    }//else
  } else
    core_error("BigFloat error: squareroot called with negative operand.",
		    __FILE__, __LINE__, true);
} //sqrt with initial approximation

//  compareMExp(x)
//    returns  1 if *this > x
//             0 if *this = x,
//            -1 if *this < x,
//
//      Main comparison method for BigFloat
//      This is called by BigFloat::compare()
//      BE CAREFUL:  The error bits are ignored!
//      Need another version if we want to take care of error bits

int BigFloatRep :: compareMExp(const BigFloatRep& x) const {
  int st = sign(m);
  int sx = sign(x.m);

  if (st > sx)
    return 1;
  else if (st == 0 && sx == 0)
    return 0;
  else if (st < sx)
    return - 1;
  else { //  need to compare m && exp
    long expDiff = exp - x.exp;

    if (expDiff > 0) //  exp > x.exp
      return cmp(chunkShift(m, expDiff), x.m);
    else if (!expDiff)
      return cmp(m, x.m);
    else  //  exp < x.exp
      return cmp(m, chunkShift(x.m, - expDiff));
  }
}

// 3/6/2000:
// This is a private function used by BigFloatRep::operator<<
// to get the exact value
// of floor(log10(M * 2^ e)) where E is an initial guess.
// We will return the correct E which satisfies
//              10^E <= M * 2^e < 10^{E+1}
// But we convert this into
//              mm <= M < 10.mm

long BigFloatRep :: adjustE( long E, BigInt M, long ee) const {
  if (M<0)
    M=-M;
  BigInt mm(1);
  if (ee > 0)
    M = (M<<static_cast<unsigned long>(ee));
  else
    mm = (mm << static_cast<unsigned long>(-ee));
  if (E > 0)
    mm *= (FiveTo(E)<< static_cast<unsigned long>(E));
  else
    M *= (FiveTo(-E) << static_cast<unsigned long>(-E));

  if (M < mm) {
    do {
      E--;
      M *= 10;
    } while (M < mm);
  } else if (M >= 10*mm) {
    mm *= 10;
    do {
      E++;
      mm *= 10;
    } while (M >= mm);
  }
  return E;
}

BigFloatRep::DecimalOutput
BigFloatRep::toDecimal(unsigned int width, bool Scientific) const {
  BigFloatRep::DecimalOutput decOut;                // to be returned
  if (err > 0) {
    decOut.isExact = false;
  } else { // err == 0
    decOut.isExact = true;
  }

  if (err > 0 && err >= abs(m)) {
    // if err is larger than mantissa, sign and significant values
    // can not be determined.
    core_error("BigFloat error: Error is too big!",
		    __FILE__, __LINE__, false);
    decOut.rep = "0.0e0";          // error is too big
    decOut.isScientific = false;
    decOut.noSignificant = 0;
    decOut.errorCode = 1;          // sign of this number is unknown
    return decOut;
  }

  decOut.sign = sign(m);
  decOut.errorCode = 0;

  BigInt M(m);                  // temporary mantissa
  long lm = bitLength(M);       // binary length of mantissa
  long e2 = bits(exp);          // binary shift length represented by exponent
  long le = clLg(err);          // binary length of err
  if (le == -1)
    le = 0;

  long L10 = 0;
  if (M != 0) {
    L10 = (long)std::floor((lm + e2) / lgTenM);
    L10 = adjustE(L10, m, e2);     // L10: floor[log10(M 2^(e2))], M != 0
  } else {
    L10 = 0;
  }
  // Convention: in the positional format, when the output is
  // the following string of 8 characters:
  //             (d0, d1, d2, d3, ".", d4, d5, d6, d7)
  // then the decimal point is said to be in the 4th position.
  // E.g., (d0, ".", d1, d2) has the decimal point in the 1st position.
  // The value of L10 says that the decimal point of output should be at
  // the (L10 + 1)st position. This is
  // true regardingless of whether M = 0 or not. For zero, we output
  // {0.0*} so L10=0.  In general, the |value| is less than 10
  // if and only if L10 is 0 and the
  // decimal point is in the 1st place.  Note that L10 is defined even if
  // the output is an integer (in which case it does not physically appear
  // but conceptually terminates the sequence of digits).

  // First, get the decimal representaion of (m * B^(exp)).
  if (e2 < 0) {
    M *= FiveTo(-e2); // M = x * 10^(-e2)
  } else if (e2 > 0) {
    M <<= e2;         // M = x * 2^(e2)
  }

  std::string decRep = M.get_str();
  // Determine the "significant part" of this string, i.e. the part which
  // is guaranteed to be correct in the presence of error,
  // except that the last digit which might be subject to +/- 1.

  if (err != 0) {     // valid = number of significant digits
    unsigned long valid = floorlg10(m) - (long)std::floor(std::log10(float(err)));
    if (decRep.length() > valid) {
      decRep.erase(valid);
    }
  }

  // All the digits in decM are correct, except the last one might
  // subject to an error +/- 1.

  if ((decRep[0] == '+') || (decRep[0] == '-')) {
    decRep.erase(0, 1);
  }

  // Second, make choice between positional representation
  // and scientific notation.  Use scientific notation when:
  // 0) if scientific notation flag is on
  // 1) err * B^exp >= 1, the error contribute to the integral part.
  // 2) (1 + L10) >= width, there is not have enough places to hold the
  //    positional representation, not including decimal point.
  // 3) The distance between the first significant digit and decimal
  //    point is too large for the width limit. This is equivalent to
  //            Either ((L10 >= 0 and (L10 + 1) > width))
  //            Or  ((L10 < 0) and (-L10 + 1) > width).

  if (Scientific ||
      ((err > 0) && (le + e2) >= 0) ||          // if err*B^exp >= 1
      ((L10 >= 0) && (L10 + 1 >= (long)width )) ||
      ((L10 < 0) && (-L10 + 1 > (long)width ))) {
    // use scientific notation
    decRep = round(decRep, L10, width);
    decOut.noSignificant = width;
    decRep.insert(1, ".");
    if (L10 != 0) {
      decRep += 'e';
      if (L10 > 0) {
        decRep += '+';
      } else { // L10 < 0
        decRep += '-';
      }
      char eBuf[48]; // enought to hold long number L10
      int ne = 0;
      if ((ne = sprintf(eBuf, "%ld", labs(L10))) >= 0) {
        eBuf[ne] = '\0';
      } else {
        //perror("BigFloat.cpp: Problem in outputing the exponent!");
        core_error("BigFloat error: Problem in outputing the exponent",
			__FILE__, __LINE__, true);
      }
      decRep += eBuf;
      decOut.isScientific = true;
    }
  } else {
    // use conventional positional notation.
    if (L10 >= 0) { // x >= 1 or x == 0 and L10 + 1 <= width
      // round when necessary
      if (decRep.length() > width ) {
        decRep = round(decRep, L10, width );
        if (decRep.length() > width ) {
          // overflow happens! use scientific notation
          return toDecimal(width, true);
        }
      }
      decOut.noSignificant = decRep.length();
      if (L10 + 1 < (long)width ) {
        decRep.insert(L10 + 1, ".");
      } else { // L10 + 1 == width
        // do nothing
      }
    } else { // L10 < 0, 0 < x < 1
      // (-L10) leading zeroes, including one to the left of decimal dot
      // need to be added in beginning.
      decRep = std::string(-L10, '0') + decRep;
      // then round when necessary
      if (decRep.length() > width ) {
        decRep = round(decRep, L10, width );
        // cannot overflow since there are L10 leading zeroes.
      }
      decOut.noSignificant = decRep.length() - (-L10);
      decRep.insert(1, ".");
    }
    decOut.isScientific = false;
  }
#ifdef CORE_DEBUG
  assert(decOut.noSignificant >= 0);
#endif

  decOut.rep = decRep;
  return decOut;
}//toDecimal

std::string BigFloatRep::round(std::string inRep, long& L10, unsigned int width) const {
  // round inRep so that the length would not exceed width.
  if (inRep.length() <= width)
    return inRep;

  int i = width; // < length
  bool carry = false;

  if ((inRep[i] >= '5') && (inRep[i] <= '9')) {
    carry = true;
    i--;
    while ((i >= 0) && carry) {
      if (carry) {
        inRep[i] ++;
        if (inRep[i] > '9') {
          inRep[i] = '0';
          carry = true;
        } else {
          carry = false;
        }
      }
      i-- ;
    }

    if ((i < 0) && carry) { // overflow
      inRep.insert(inRep.begin(), '1');
      L10 ++;
      width ++;
    }
  }

  return inRep.substr(0, width);
}//round(string,width)


// This function fromString(str, prec) is similar to the
//      constructor Real(char * str, extLong prec)
// See the file Real.cc for the differences

void BigFloatRep :: fromString(const char *str, const extLong & prec ) {
  // NOTE: prec defaults to defBigFloatInputDigits (see BigFloat.h)
  // check that prec is not INFTY
  if (prec.isInfty())
    core_error("BigFloat error: infinite precision not allowed",
			__FILE__, __LINE__, true);

  const char *e = strchr(str, 'e');
  int dot = 0;
  long e10 = 0;
  if (e != NULL)
    e10 = atol(e+1);    // e10 is decimal precision of the input string
  // i.e., input is A/10^{e10}.
  else {
    e = str + strlen(str);
#ifdef CORE_DEBUG
    assert(*e == '\0');
#endif

  }

  const char *p = str;
  if (*p == '-' || *p == '+')
    p++;
  m = 0;
  exp = 0;

  for (; p < e; p++) {
    if (*p == '.') {
      dot = 1;
      continue;
    }
    m = m * 10 + (*p - '0');
    if (dot)
      e10--;
  }

  BigInt one = 1;
  long t = (e10 < 0) ? -e10 : e10;
  BigInt ten = FiveTo(t) * (one << static_cast<unsigned long>(t));

  // HERE IS WHERE WE USE THE SYSTEM CONSTANT
  //           defBigFloatInputDigits
  // Note: this constant is rather similar to defInputDigits which
  //     is used by Real and Expr for controlling
  //     input accuracy.  The difference is that defInputDigits can
  //     be CORE_INFTY, but defBigFloatInputDigits must be finite.

  if (e10 < 0)
    div(m, ten, CORE_posInfty, 4 * prec);
  else
    m *= ten;
  if (*str == '-')
    m = -m;
}//BigFloatRep::fromString

std::istream& BigFloatRep :: operator >>(std::istream& i) {
  int size = 20;
  char *str = new char[size];
  char *p = str;
  char c;
  int d = 0, e = 0, s = 0;
  // d=1 means dot is found
  // e=1 means 'e' or 'E' is found
  //  int done = 0;

  // Chen Li: fixed a bug, the original statement is
  //  for (i.get(c); c == ' '; i.get(c));
  // use isspace instead of testing c == ' ', since it must also
  // skip tab, catridge/return, etc.
  // Change to:
  //  int status;
  do {
    c = i.get();
  } while (isspace(c)); /* loop if met end-of-file, or
                               char read in is white-space. */
  // Chen Li, "if (c == EOF)" is unsafe since c is of char type and
  // EOF is of int tyep with a negative value -1
  if (i.eof()) {
    i.clear(std::ios::eofbit | std::ios::failbit);
    return i;
  }

  // the current content in "c" should be the first non-whitespace char
  if (c == '-' || c == '+') {
    *p++ = c;
    i.get(c);
  }

  for (; isdigit(c) || (!d && c=='.') ||
       (!e && ((c=='e') || (c=='E'))) || (!s && (c=='-' || c=='+')); i.get(c)) {
    if (!e && (c == '-' || c == '+'))
      break;
    // Chen Li: put one more rule to prohibite input like
    //  xxxx.xxxe+xxx.xxx:
    if (e && (c == '.'))
      break;
    if (p - str == size) {
      char *t = str;
      str = new char[size*2];
      memcpy(str, t, size);
      delete [] t;
      p = str + size;
      size *= 2;
    }
#ifdef CORE_DEBUG
    assert((p-str) < size);
#endif

    *p++ = c;
    if (c == '.')
      d = 1;
    // Chen Li: fix a bug -- the sign of exponent can not happen before
    // the character "e" appears! It must follow the "e' actually.
    //    if (e || c == '-' || c == '+') s = 1;
    if (e)
      s = 1;
    if ((c == 'e') || (c=='E'))
      e = 1;
  }

  // chenli: make sure that the p is still in the range
  if (p - str >= size) {
    int len = p - str;
    char *t = str;
    str = new char[len + 1];
    memcpy(str, t, len);
    delete [] t;
    p = str + len;
  }

#ifdef CORE_DEBUG
  assert(p - str < size);
#endif

  *p = '\0';
  i.putback(c);
  fromString(str);
  delete [] str;
  return i;
}//operator >>


// BigFloatRep::toDouble()
//      converts the BigFloat to a machine double
//      This is a dangerous function as the method
//      is silent when it does not fit into a machine double!
// ToDo: fix this by return a machine NaN, +/- Infinity, +/- 0,
//      when appropriate.
//      Return NaN when error is larger than mantissa
//      Return +/- Infinity if BigFloat is too big
//      Return +/- 0 if BigFloat is too small
#ifdef _MSC_VER
#pragma warning(disable: 4723)
#endif
double BigFloatRep :: toDouble() const {
  if (m == 0)
    return (sign(m) * 0.0);

  long e2 = bits(exp);
  long le = clLg(err);  // if err=0, le will be -1
  if (le == -1)
    le = 0;

  BigInt M = m >> static_cast<unsigned long>(le);// remove error bits in mantissa

  // Below, we want to return NaN by computing 0.0/0.0.
  // To avoid compiler warnings about divide by zero, we do this:

  double foolCompilerZero;
  foolCompilerZero = 0.0;

  // COMMENT: we should directly store the
  //    special IEEE values NaN, +/-Infinity, +/-0 in the code!!

  if (M == 0)
    return ( 0.0/foolCompilerZero ) ; // return NaN

  e2 += le;             // adjust exponent corresponding to error bits

  int len = bitLength(M) - 53;  // this is positive if M is too large

  if (len > 0) {
    M >>= len;
    e2 += len;
  }

  double tt = doubleValue(M);

  int ee = e2 + bitLength(M) - 1; // real exponent.

  if (ee >= 1024)       // overflow!
    return (  sign(m)/foolCompilerZero  );      // return a signed infinity

  if (ee <= -1075)      // underflow!
    // NOTE: if (-52 < ee <= 0) get denormalized number
    return ( sign(m) * 0.0 );  // return signed zero.

  // Execute this loop if e2 < 0;
  for (int i = 0; i > e2; i--)
    tt /= 2;

  // Execute this loop if e2 > 0;
  for (int j = 0; j < e2; j++)
    tt *= 2;

  return tt;
}//toDouble
#ifdef _MSC_VER
#pragma warning(default: 4723)
#endif
BigInt BigFloatRep::toBigInt() const {
  long e2 = bits(exp);
  long le = clLg(err);
  if (le == -1)
    le = 0;
#ifdef CORE_DEBUG
  assert (le >= 0);
#endif

  BigInt M = m >> static_cast<unsigned long>(le); // discard the contaminated bits.
  e2 += le;           // adjust the exponent

  if (e2 < 0)
    return M >> static_cast<unsigned long>(-e2);
  else if (e2 > 0)
    return M << static_cast<unsigned long>(e2);
  else
    return M;
}

long BigFloatRep :: toLong() const {
  // convert a BigFloat to a long integer, rounded toward -\infty.
  long e2 = bits(exp);
  long le = clLg(err);
#ifdef CORE_DEBUG
  assert (le >= 0);
#endif

  BigInt M = m >> static_cast<unsigned long>(le); // discard the contaminated bits.
  e2 += le;           // adjust the exponent
  long t;
  if (e2 < 0)
    t = ulongValue(M >> static_cast<unsigned long>(-e2));
  else if (e2 > 0)
    t = ulongValue(M << static_cast<unsigned long>(e2));
  else
    t = ulongValue(M);
  // t = M.as_long();
  // Note: as_long() will return LONG_MAX in case of overflow.

  return t;
}

// pow(r,n) function for BigFloat
// Note: power(r,n) calls pow(r,n)
BigFloat pow(const BigFloat& r, unsigned long n) {
  if (n == 0)
    return BigFloat(1);
  else if (n == 1)
    return r;
  else {
    BigFloat x = r;
    while ((n % 2) == 0) { // n is even
      x *= x;
      n >>= 1;
    }
    BigFloat u = x;
    while (true) {
      n >>= 1;
      if (n == 0)
        return u;
      x *= x;
      if ((n % 2) == 1) // n is odd
        u *= x;
    }
    //return u; // unreachable
  }
}//pow

// experimental
BigFloat root(const BigFloat& x, unsigned long k,
         const extLong& a, const BigFloat& A) {
  if (x.sign() == 0) {
    return BigFloat(0);
  } else if (x == 1) {
    return BigFloat(1);
  } else  {
    BigFloat q, del, zz;
    BigFloat z = A;
    BigFloat bk = long(k);
    for (; ;) {
      zz = pow(z, k-1);
      q = x.div(zz, a);
      q.makeExact();
      del = z - q;
      del.makeExact();
      if (del.MSB() < -a)
        break;
      z = ((bk-1)*z + q).div(bk, a);
          // newton's iteration: z_{n+1}=((k-1)z_n+x/z_n^{k-1})/k
      z.makeExact();
    }
    return z;
  }
}//root

} //namespace CORE
