/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: Real.cpp
 * 
 * Synopsis: The Real class is a superclass for all the number 
 *           systems in the Core Library (int, long, float, double,
 *           BigInt, BigRat, BigFloat, etc)
 *           
 * Written by 
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu> 
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#include "Real.h"

CORE_BEGIN_NAMESPACE

#ifndef CORE_ENABLE_INLINES
#include "Real.inl"
#endif

const Real& Real::getZero() {
  static Real Zero(0);
  return Zero;
}

BigInt floor(const Real& r, Real &sub) {
  BigInt f = r.approx(CORE_INFTY, 2).getBigFloat().toBigInt();
  sub = r-f;
  // Adjustment
  if (sub<0)
    ++sub, --f;
  if (sub>=1)
    --sub, ++f;
  assert(sub >=0 && sub<1);
  return f;
}

Real pow(const Real& r, unsigned long n) {
  if (n == 0) 
    return Real(1);
  else if (n == 1)
    return r;
  else {
    Real x = r;
    while ((n % 2) == 0) { // n is even
      x *= x;
      n >>= 1;
    }
    Real u = x;
    while (true) {
      n >>= 1;
      if (n == 0) return u;
      x *= x;
      if ((n % 2) == 1) // n is odd
        u *= x;
    }
    return u;
  }
}//pow

extern BigInt FiveTo(unsigned long exp);

// Constructor for Real from String
//   --Default value of the argument "prec" is defInputDigits
//   --If prec is CORE_posInfty, then the input is
//	read in exactly.  Otherwise, we convert to a RealBigFloat
//	with absolute error at most 10^{-prec}

// Constructor Real( char *str, extLong& prec)
//	is very similar to
//		BigFloatRep::fromString( char *str, extLong& prec);
// Differences:
//	In BigFloat(str, prec), the value of prec cannot be infinity, and
//		it defaults to defBigFloatInputDigits;
//	In Real(str, prec), the value of prec is allowed to be infinity, and
//		it defaults to defInputDigits.
//
// Why do we have the two versions?  It allows us to use the BigFloat class
//	directly, without relying on Real class.

Real::Real(const char *str, const extLong& prec )
 	// NOTE: prec defaults to defInputDigits (see Real.h)
{ 
//	8/8/01, Chee and Zilin: add a new rational string format:
//		this format is indicated by the presence of a slash "/"
//		Moreover, the value of prec is ignored (basically
//		assumed to be infinity).

  if (strchr(str, '/') != NULL) {	// this is a rational number
    rep = new RealBigRat(BigRat(str));
    return;
  }

  const char *e = strchr(str, 'e');
  int dot = 0;
  long e10 = 0;
  if (e != NULL)
    e10 = atol(e+1);	// e10 is decimal precision of the input string
			// i.e., input is A/10^{e10}.
  else {
    e = str + strlen(str);
#ifdef DEBUG
    assert(*e == '\0');
#endif 
  }

  const char *p = str;
  if (*p == '-' || *p == '+') p++;
  BigInt m = 0;
  
  for (; p < e; p++) {
    if (*p == '.') {
      dot = 1;
      continue;
    }
    m = m * 10 + (*p - '0');
    if (dot) e10--;
  }
  
  long t = (e10 < 0) ? -e10 : e10;
  BigInt one = 1;
  BigInt ten = FiveTo(t) * (one << t);
  if (*str == '-') m = -m;
  if (e10 >= 0) {
    // convert exactly from integer numbers 
    m *= ten;
    rep = new RealBigInt(m);
  } else { // e10 < 0,  fractional numbers
    // HERE IS WHERE WE USE THE SYSTEM CONSTANT
    //	       defInputDigits
    // Note: defInputDigits should be at least log_2(10).
    //       We default defInputDigits to 4.
    BigRat r(m, ten);
    if (prec.isInfty()) { // convert exactly! to a big rational
      rep = new RealBigRat(r);
    } else { 
      // convert approximately, to a BigFloat within the 
      // specified precision:     
      // BigFloat bf(r, CORE_posInfty, prec * lgTenM) ;
      BigFloat bf(r, CORE_posInfty, prec * 4) ;
      rep = new RealBigFloat(bf);
    }
  }
}// Real(str, prec)

// The operator >>(i,x) calls the constructor Real(char*)
std::istream& operator >>(std::istream& i, Real& x)
{
  int size = 20;
  char *str = new char[size];
  char *p = str;
  char c;
  int d = 0, e = 0, s = 0;
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
  // Chen Li,
  // original "if (c == EOF) ..." is unsafe since c is of char type and
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
	 (!e && c=='e') || (!s && (c=='-' || c=='+')); i.get(c)) {
    if (!e && (c == '-' || c == '+')) break;
    // Chen Li: put one more rule to prohibite input like 
    //  xxxx.xxxe+xxx.xxx:
    if (e && (c == '.')) break;
    if (p - str == size) {
      char *t = str;
      str = new char[size*2];
      memcpy(str, t, size);
      delete [] t;
      p = str + size;
      size *= 2;
    }
#ifdef DEBUG
    assert((p-str) < size);
#endif
    *p++ = c;
    if (c == '.')                  d = 1;
    // Chen Li: fix a bug -- the sign of exponent can not happen before
    // the character "e" appears! It must follow the "e' actually. 
    //    if (e || c == '-' || c == '+') s = 1;
    if (e) s = 1;
    if (c == 'e')                  e = 1;
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

#ifdef DEBUG
  assert(p - str < size);
#endif
  *p = '\0';
  i.putback(c);
  // old: x = Real(str, i.precision()); // use precision of input stream.
  x = Real(str);  // default precision = defInputDigits
  delete [] str;
  return i;
}//operator >> (std::istream&, Real&)

//  stream : This function is common to all Realbase_for<> types.
template <class T>
std::ostream& Realbase_for<T>::operator <<(std::ostream& o) const
{
  o << ker;
  return o;
}

//  constructor for RealLong
template<>
Realbase_for<long>::Realbase_for(const long &l)
  : ker(l)
{
  mostSignificantBit = ((ker != 0) ? extLong(flrLg(ker)) : CORE_negInfty); 
  //  This computes the bit length of "ker" minus 1,
  //  i.e., floor(log_2(|ker|)) .
}

//  constructor for RealDouble
template<>
Realbase_for<double>::Realbase_for(const double& d)
  : ker(d)
{
  mostSignificantBit = BigFloat(ker).MSB();
}

//  constructor for RealBigInt
template<>
Realbase_for<BigInt>::Realbase_for<BigInt>(const BigInt& I)
  : ker(I)
{
  mostSignificantBit = (sign(ker)? extLong(floorLg(ker)) : CORE_negInfty);
}

//  constructor for RealBigFloat
template<>
Realbase_for<BigFloat>::Realbase_for(const BigFloat& B)
  : ker(B)
{
  mostSignificantBit = ker.MSB();
}

//  constructor for RealBigRat
template<>
Realbase_for<BigRat>::Realbase_for(const BigRat& R) : ker(R)
{
  // MSB of a rational x/y is given by floorLg(|x/y|)
  BigInt x = ker.numerator();
  BigInt y = ker.denominator();
  if (ker.sign()) {
    mostSignificantBit = extLong(floorLg(x) - floorLg(y));
    x.abs();
    if ((y << mostSignificantBit.asLong()) > x) 
      mostSignificantBit = mostSignificantBit - 1;
  } else
    mostSignificantBit = CORE_negInfty;
  /*
  mostSignificantBit = ker.sign() ? \
       extLong(floorLg(x) - floorLg(y)) : CORE_negInfty;

  // This gives us an approximation to msb that could off by 1 in
  // one direction.  So we next adjust for this possibility:
  // The exact value of msb(x/y) is given by
  //   y.2^msb <= x < y.2^{msb+1}.

  // 5/16/02: fixed a bug in logic (Pion/Zilin/Chee)
  x.abs();
  if ((y << mostSignificantBit.asLong()) > x) 
       mostSignificantBit = mostSignificantBit - 1;
  */
}

CORE_END_NAMESPACE
