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
  long lm = m.bit_length();
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
    d += bigIntToLong(t%10);
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


//  class BigFloat
//  constructors
CORE_INLINE BigFloat::BigFloat(BigFloatRep * r) { //private
  this->rep = r;
  rep->refCount ++;
}

CORE_INLINE BigFloat::BigFloat() {
  rep = new BigFloatRep(CORE_BIGINT_ZERO, 0, 0);
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(int i) {
  rep = new BigFloatRep(i);
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(long l) {
  rep = new BigFloatRep(l);
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(double d) {
  rep = new BigFloatRep(d);
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(const BigInt& I, unsigned long u, long l) {
                // default values: u=0, l=0 
  rep = new BigFloatRep(I, u, l);
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(const char *str) {
  rep = new BigFloatRep(str);
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(const BigFloat& x) {
  rep = x.rep;
  rep->refCount++;
}

CORE_INLINE BigFloat::BigFloat(const BigRat& R, const extLong& r, const extLong& a) {
        // default values: r= defRelPrec, a= defAbsPrec
  rep = new BigFloatRep();
  rep->refCount++;
  rep->approx(R, r, a);
}

CORE_INLINE BigRat BigFloat::BigRatize() const {
  return rep->BigRatize();
}

//  the destructor
CORE_INLINE BigFloat::~BigFloat() {
  if (--rep->refCount == 0)
    delete rep;
}

//  assignment operator
CORE_INLINE BigFloat& BigFloat::operator= (const BigFloat& x) {
  if (this == &x)               //Bug fixed: avoid self-assignment
    return *this;               //Zilin Du, 07/11/01
  if (--rep->refCount == 0)
  delete rep;
  
  rep = x.rep;
  rep->refCount++;
  
  return *this;
}

//  approximation
CORE_INLINE void BigFloat::approx(const BigInt& I, const extLong& r, const extLong& a) {
  if ((rep->refCount) > 1) {//  *rep is shared
    --rep->refCount;
    
    rep = new BigFloatRep();
    rep->refCount++;
  }
  rep->trunc(I, r, a);
}

CORE_INLINE void BigFloat::approx(const BigFloat& B, const extLong& r, const extLong& a) {
  if ((rep->refCount) > 1) {//  *rep is shared
    --rep->refCount;
    
    rep = new BigFloatRep();
    rep->refCount++;
  }
  
  rep->approx(*B.rep, r, a);
}

CORE_INLINE void BigFloat::approx(const BigRat& R, const extLong& r, const extLong& a) {
  if (rep->refCount > 1) {//  *rep is shared
      --rep->refCount;
    
      rep = new BigFloatRep();
      rep->refCount++;
  }
  
  rep->approx(R, r, a);
}

//  unary minus
CORE_INLINE BigFloat BigFloat::operator- () const {
  return BigFloat(- rep->m, rep->err, rep->exp);
}

//  arithmetics
CORE_INLINE BigFloat operator+ (const BigFloat& x, const BigFloat& y) {
  BigFloat z;    
  z.rep->add(*x.rep, *y.rep);
  return z;
}

CORE_INLINE BigFloat operator- (const BigFloat& x, const BigFloat& y) {
  BigFloat z;  
  z.rep->sub(*x.rep, *y.rep);
  return z;
}

CORE_INLINE BigFloat operator* (const BigFloat& x, const BigFloat& y) {
  BigFloat z;    
  z.rep->mul(*x.rep, *y.rep);
  return z;
}

CORE_INLINE BigFloat operator/ (const BigFloat& x, const BigFloat& y) {
  BigFloat z;  
  z.rep->div(*x.rep, *y.rep, defBFdivRelPrec);
  return z;
}

CORE_INLINE BigFloat BigFloat::div(const BigFloat& x, const extLong& r) const {
  BigFloat y;
  y.rep->div(*rep, *x.rep, r);
  return y;
}

//  squareroot
CORE_INLINE BigFloat BigFloat::sqrt(const extLong& a) const {
  BigFloat x;  
  x.rep->sqrt(*rep, a);
#ifdef DEBUG  
  assert((x * x - (*this)).abs() < BigFloat(ldexp(2.0, -a.asLong())));
#endif
  return x;
}

//  squareroot with initial approximation A
CORE_INLINE BigFloat BigFloat :: sqrt(const extLong& a, const BigFloat& A) const
{
  BigFloat x;  
  x.rep->sqrt(*rep, a, A);
#ifdef DEBUG  
  assert((x * x - (*this)).abs() < BigFloat(ldexp(2.0, -a.asLong())));
#endif
  return x;
}

// sqrt to defAbsPrec:
CORE_INLINE BigFloat sqrt(const BigFloat& x) {
  return x.sqrt(defBFsqrtAbsPrec);
}

//  comparisons
CORE_INLINE int BigFloat::compare(const BigFloat& x) const {
  return rep->compareMExp(*x.rep);
}

CORE_INLINE bool operator== (const BigFloat& x, const BigFloat& y) {
  return x.compare(y) == 0;
}

CORE_INLINE bool operator!= (const BigFloat& x, const BigFloat& y) {
  return x.compare(y) != 0;
}

CORE_INLINE bool operator< (const BigFloat& x, const BigFloat& y) {
  return x.compare(y) < 0;
}

CORE_INLINE bool operator<= (const BigFloat& x, const BigFloat& y) {
  return x.compare(y) <= 0;
}

CORE_INLINE bool operator> (const BigFloat& x, const BigFloat& y) {
  return x.compare(y) > 0;
}

CORE_INLINE bool operator>= (const BigFloat& x, const BigFloat& y) {
  return x.compare(y) >= 0;
}

//  arithmetic and assignment opeartors
CORE_INLINE BigFloat& BigFloat::operator+= (const BigFloat& x) {
  BigFloat z;
  z.rep->add(*rep, *x.rep);
  
  if (--rep->refCount == 0)
  	delete rep;
  
  rep = z.rep;
  rep->refCount++;
  
  return *this;
}

CORE_INLINE BigFloat& BigFloat::operator-= (const BigFloat& x) {
  BigFloat z;
  z.rep->sub(*rep, *x.rep);
  
  if (--rep->refCount == 0)
  	delete rep;
  
  rep = z.rep;
  rep->refCount++;
  
  return *this;
}

CORE_INLINE BigFloat& BigFloat::operator*= (const BigFloat& x) {
  BigFloat z;
  z.rep->mul(*rep, *x.rep);
  
  if (--rep->refCount == 0)
  	delete rep;
  
  rep = z.rep;
  rep->refCount++;
  
  return *this;
}

CORE_INLINE BigFloat& BigFloat::operator/= (const BigFloat& x) {
  BigFloat z;
  z.rep->div(*rep, *x.rep, defBFdivRelPrec);
  
  if (--rep->refCount == 0)
  	delete rep;
  
  rep = z.rep;
  rep->refCount++;
  
  return *this;
}

//  cast operators
CORE_INLINE double BigFloat::toDouble() const {
  return rep->toDouble();
}

CORE_INLINE float BigFloat::toFloat() const {
  return (float)rep->toDouble();
}

CORE_INLINE BigInt BigFloat::toBigInt() const {
  return rep->toBigInt();
}

CORE_INLINE long BigFloat::toLong() const {
  long l = rep->toLong();
  if ((l == LONG_MAX) || (l == LONG_MIN))
    return l; // return the overflown value.
  if ((sign() < 0) && (*this != BigFloat(l))) {
    // a negative value not exactly rounded.
    l--; // rounded to floor.
  }
  return l;
}

CORE_INLINE int BigFloat::toInt() const {
  return (int)rep->toLong();
}

//  builtin function
CORE_INLINE bool BigFloat::isExact() const {
  return rep->err == 0;
}

CORE_INLINE extLong BigFloat::lMSB() const {
  return rep->lMSB();
}

CORE_INLINE extLong BigFloat::uMSB() const {
  return rep->uMSB();
}

CORE_INLINE extLong BigFloat::MSB() const {
  return rep->MSB();
}

CORE_INLINE extLong BigFloat::flrLgErr() const {
  return rep->flrLgErr();
}

CORE_INLINE extLong BigFloat::clLgErr() const {
  return rep->clLgErr();
}

CORE_INLINE bool BigFloat::isZeroIn() const {
  return rep->isZeroIn();
}

CORE_INLINE int BigFloat::sign() const {
  return rep->signM();
}

CORE_INLINE int sign(const BigFloat& x) {
  return x.sign();
}

CORE_INLINE BigFloat BigFloat::abs() const {
  if (rep->signM() >= 0)
    return BigFloat(*this);
  else
    return BigFloat(- rep->m, rep->err, rep->exp);
}

CORE_INLINE BigFloat abs(const BigFloat& x) {
  return x.abs();
}

//  stream operators
CORE_INLINE std::ostream& operator<< (std::ostream& o, const BigFloat& x) {
  x.rep->operator<<(o);
  return o;
}

CORE_INLINE std::istream& operator>> (std::istream& i, BigFloat& x) {
  x.rep->operator>>(i);  
  return i;
}

CORE_INLINE void BigFloat::dump() const { 
  rep->dump(); 
}

#endif

