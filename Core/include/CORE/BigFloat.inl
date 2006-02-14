/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: BigFloat.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

// class BigFloatRep

// static member
CORE_INLINE void BigFloatRep::error(const char* msg)
{
  //  output error message to stderr
  std::cerr << "BigFloat Error : " << msg << std::endl;
  abort();
  exit(1);
}

CORE_INLINE long BigFloatRep::chunkCeil(long bits) {
  if (bits > 0)
    return (bits - 1) / CHUNK_BIT + 1;
  else
    return - (- bits) / CHUNK_BIT;
}//chunkCeil

CORE_INLINE long BigFloatRep::chunkFloor(long bits) {
  if (bits >= 0)
    return bits / CHUNK_BIT;
  else
    return - (- bits - 1) / CHUNK_BIT - 1;
}//chunkFloor

// bits(c) returns the number of bits in c chunks:
CORE_INLINE long BigFloatRep::bits(long chunks) {
  return CHUNK_BIT * chunks;
}//bits

CORE_INLINE BigInt BigFloatRep::chunkShift(const BigInt& x, long s) {
  if (!s || sign(x) == 0)
    return x;
  else if (s > 0)
    //  shift left
    if (sign(x) > 0)
      return x << bits(s);
    else //  x < 0
      return - ((-x) << bits(s));
  else //  shift right
    if (sign(x) > 0)
      return x >> bits(-s);
    else //  x < 0
      return - ((-x) >> bits(-s));
}//chunkShift

CORE_INLINE BigFloatRep* BigFloatRep::exp2(int e) {
  long ee;  // this is going to be the exponent
  if (e >= 0) 
    ee = e/CHUNK_BIT;
  else 
    ee = - ((-e + CHUNK_BIT -1)/CHUNK_BIT);
   
  int rem = e - (ee * CHUNK_BIT);     // Assert: 0 <= rem < CHUNK_BIT
 
  return new BigFloatRep((1<<rem), 0, ee); // assume CHUNK_BIT is less than int wid
}

//  constructor
CORE_INLINE BigFloatRep::BigFloatRep(const BigInt& I, unsigned long u, long l)
  : m(I), err(u), exp(l), refCount(0)
{}

CORE_INLINE BigFloatRep::BigFloatRep(const char *str)
  : m(0), err(0), exp(0), refCount(0)
{
  fromString(str);
}

CORE_INLINE BigRat BigFloatRep::BigRatize() const {
  if (exp >= 0)
    return BigRat(chunkShift(m, exp), 1);
  else
    return BigRat(m, chunkShift(1, - exp));  
}

//  the destructor
CORE_INLINE BigFloatRep::~BigFloatRep()
{}

CORE_INLINE void BigFloatRep::approx(const BigRat& R, const extLong& r, const extLong& a) {
  div(R.numerator(), R.denominator(), r, a);
}

//  eliminate trailing zeroes
CORE_INLINE void BigFloatRep::eliminateTrailingZeroes() {
  // eliminate trailing 0's    -- IP 10/9/98
  if (err == 0 && m != 0) {
    while ((m & ((1 << CHUNK_BIT) - 1)) == 0) {
      m >>= CHUNK_BIT;
      exp++;
    }
  }
}

//  bultin functions
CORE_INLINE extLong BigFloatRep::lMSB() const {
  if (!isZeroIn())
    return extLong(floorLg(abs(m) - err)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

CORE_INLINE extLong BigFloatRep::uMSB() const {
  return extLong(floorLg(abs(m) + err)) + bits(exp);
}

CORE_INLINE extLong BigFloatRep::MSB() const {
  // Note : MSB is undefined if it's not exact.
  if (sign(m))          // sign(m) is non-zero
    return extLong(floorLg(m)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

CORE_INLINE extLong BigFloatRep::flrLgErr() const {
  if (err)
    return extLong(flrLg(err)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

CORE_INLINE extLong BigFloatRep::clLgErr() const {
  if (err)
    return extLong(clLg(err)) + bits(exp);
  else
    return extLong(CORE_negInfty);
}

// isZero() = true iff zero is inside the interval of BigFloat:
CORE_INLINE bool BigFloatRep::isZeroIn() const {
  if (err == 0) return (m == 0);	//Nov 6, 2002: bug fix!
  long lm = m.bitLength();
  if (lm > CHUNK_BIT+2)
    return false;   // since err < 4 * 2^{CHUNK_BIT}
  else
    return (abs(m) <= BigInt(err));
}

CORE_INLINE int BigFloatRep::signM() const {
  return sign(m);
}

CORE_INLINE double BigFloatRep::lg10(BigInt x) {
  if (x == 0)
    return 0;
  
  BigInt t(abs(x));  
  long l = -1;
  double d = 0;
  
  while (t > 0) {
    l++;
    d /= 10;
    d += BigInt(t%10).ulongValue();
    t /= 10;
  }
  return log10(d) + l;
}

// this is a simpler form of lg10()
CORE_INLINE long BigFloatRep::floorlg10(BigInt x)
{
  if (x == 0)
    return 0;
  BigInt t(abs(x));
  long l = -1;
  
  while (t > 0) {
    l++;
    t /= 10;
  }
  return l;
}

CORE_INLINE std::ostream& BigFloatRep::operator<<(std::ostream& o) const {
  bool sci = (o.flags() & std::ios ::scientific);
  BigFloatRep::DecimalOutput r = toDecimal(o.precision(), sci);
  if (r.sign == -1)
    o << "-";
  o << r.rep;
  return o;
}

/* Returns a std::string with precision and format specified
   Works as cout << with the exception that if the output
   contains any error it returns a NULL
   Joaquin Grech 31/5/03
   */
CORE_INLINE std::string BigFloatRep::toString(long prec, bool sci) const {
  BigFloatRep::DecimalOutput r = toDecimal(prec, sci);
  
  if (r.errorCode == 0) {
	if (r.sign == -1)
		return std::string("-")+r.rep;
	else
		return r.rep;
  }
  return NULL;
}

CORE_INLINE void BigFloatRep::dump() const {
  std::cout << "---- BFRep: " << this << " ----" << std::endl;
  std::cout << "  BF value: ";
  this->operator<<(std::cout);
  std::cout <<  std::endl;
  std::cout << "  m = " << m << std::endl;
  std::cout << "  err = " << err << std::endl;
  std::cout << "  exp = " << exp << std::endl;
  std::cout << " -- End of BFRep " << this << " -- " << std::endl;
}


#endif

