/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: Expr.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

// Constructors and Destructor
CORE_INLINE Expr::Expr()
{ rep = new ConstDoubleRep(); }

CORE_INLINE Expr::Expr(int i)
{ rep = new ConstDoubleRep(i); }

CORE_INLINE Expr::Expr(unsigned int i)
{ rep = new ConstDoubleRep(i); }

CORE_INLINE Expr::Expr(long i)
{ rep = new ConstRealRep(Real(i)); }

CORE_INLINE Expr::Expr(unsigned long int i)
{ rep = new ConstRealRep(Real(i)); }

/** \note {the results of this constructor can be somewhat unpredictable.
 *  One might assume that new Expr(.1) is exactly equal to .1, but it is
 *  actually equal to
 *      .1000000000000000055511151231257827021181583404541015625.
 *  This is so because .1 cannot be represented exactly as a double
 *  (or, for that matter, as a binary fraction of any finite length). 
 *  Thus, the long value that is being passed in to the constructor is not
 *  exactly equal to .1, appearances nonwithstanding. */
CORE_INLINE Expr::Expr(float d)
{ rep = new ConstDoubleRep(d); }

CORE_INLINE Expr::Expr(double d)
{ rep = new ConstDoubleRep(d); }

CORE_INLINE Expr::Expr(const BigInt& I)
{ rep = new ConstRealRep(Real(I)); }

CORE_INLINE Expr::Expr(const BigRat& r)
{ rep = new ConstRealRep(Real(r)); }

CORE_INLINE Expr::Expr(const BigFloat& B)
{ rep = new ConstRealRep(Real(B)); }

/** The (String) constructor, on the other hand, is perfectly predictable:
 * new Expr(".1") is exactly equal to .1, as one would expect. Therefore,
 * it is generally recommended that the (String) constructor be used in
 * preference to the (double) constructor. */
CORE_INLINE Expr::Expr(const char *s, const extLong & p)
{ rep = new ConstRealRep(Real(s, p)); }

CORE_INLINE Expr::Expr(const Real& r)
{ rep = new ConstRealRep(r); }

CORE_INLINE Expr::Expr(const Expr& e)
{ rep = e.rep;  rep->incRefCount(); }

CORE_INLINE Expr::~Expr()
{ rep->decRefCount(); }

// Approximation Functions
CORE_INLINE const Real & Expr::approx(const extLong& relPrec, 
                                      const extLong& absPrec) const
{ return rep->getAppValue(relPrec, absPrec); }

// Help Functions
CORE_INLINE long Expr::getExponent() const
{ return getBigFloat().exp(); }

CORE_INLINE BigInt Expr::getMantissa() const
{ return getBigFloat().m(); }

CORE_INLINE BigFloat Expr::getBigFloat() const
{ return rep->getBigFloat(); }

CORE_INLINE int Expr::sign() const
{ return getSign(); }

// Conversion Functions
CORE_INLINE int Expr::toInt() const
{ return (this->approx(64, 1024)).toInt(); }

CORE_INLINE long Expr::toLong() const
{ return (this->approx(64, 1024)).toLong(); }

CORE_INLINE float Expr::toFloat() const
{ return (this->approx(53, 1024)).toFloat(); }

/** chen - use equivalent precision (rel:53, abs: 1024)
    as in IEEE double. enforce an evaluation in case
    before this has been done before casting. */
CORE_INLINE double Expr::toDouble() const 
{ return (this->approx(53, 1024)).toDouble(); }

CORE_INLINE const char* Expr::toString() const
{ return rep->toString(); }

// Assignment Operators
CORE_INLINE Expr& Expr::operator=(const Expr& e) {
  if (this == &e) return *this;
  e.rep->incRefCount();
  rep->decRefCount();
  rep = e.rep;
  return *this;
}

CORE_INLINE Expr& Expr::operator+=(const Expr& e) {
  ExprRep *old = rep;
  rep = new AddRep(rep, e.rep);
  old->decRefCount();
  return *this;
}

CORE_INLINE Expr& Expr::operator-=(const Expr& e) {
  ExprRep *old = rep;
  rep = new SubRep(rep, e.rep);
  old->decRefCount();
  return *this;
}

CORE_INLINE Expr& Expr::operator*=(const Expr& e) {
  ExprRep *old = rep;
  rep = new MultRep(rep, e.rep);
  old->decRefCount();
  return *this;
}

CORE_INLINE Expr& Expr::operator/=(const Expr& e) {
  if (e == 0) {
    std::cerr << " ERROR : division by zero ! " << std::endl;
    abort();
  }
  ExprRep *old = rep;
  rep = new DivRep(rep, e.rep);
  old->decRefCount();
  return *this;
}

// Increment, Decrement, and Unary Minus Operators
CORE_INLINE Expr& Expr::operator++()
{ *this += 1; return *this; }

CORE_INLINE Expr Expr::operator++(int)
{ Expr t = *this; *this += 1; return t; }

CORE_INLINE Expr& Expr::operator--()
{ *this -= 1; return *this; }

CORE_INLINE Expr Expr::operator--(int)
{ Expr t = *this; *this -= 1; return t; }

CORE_INLINE Expr Expr::operator-() const 
{ return Expr(new NegRep(rep)); }

// Arithematic Operators
CORE_INLINE Expr operator+(const Expr& e1, const Expr& e2)
{ return Expr(new AddRep(e1.rep, e2.rep)); }

CORE_INLINE Expr operator-(const Expr& e1, const Expr& e2)
{ return Expr(new SubRep(e1.rep, e2.rep)); }

CORE_INLINE Expr operator*(const Expr& e1, const Expr& e2)
{ return Expr(new MultRep(e1.rep, e2.rep)); }

CORE_INLINE Expr operator/(const Expr& e1, const Expr& e2)
{
  if (e2 == 0) {
    std::cerr << " ERROR : division by zero ! " << std::endl;
    abort();
  }
  return Expr(new DivRep(e1.rep, e2.rep));
}

CORE_INLINE Expr sqrt(const Expr& e)
{ return Expr(new SqrtRep(e.rep)); }

/// Comparison Operators
/** this is inefficient if you compare to zero: 
 *  e.g., if (e != 0) {...} use sign(e) == 0 instead */
CORE_INLINE bool operator==(const Expr& e1, const Expr& e2)
{ return !compare(e1, e2); }

CORE_INLINE bool operator!=(const Expr& e1, const Expr& e2)
{ return compare(e1, e2); }

CORE_INLINE bool operator< (const Expr& e1, const Expr& e2)
{ return compare(e1, e2) < 0; }

CORE_INLINE bool operator<=(const Expr& e1, const Expr& e2)
{ return compare(e1, e2) <= 0; }

CORE_INLINE bool operator> (const Expr& e1, const Expr& e2)
{ return compare(e1, e2) > 0; }

CORE_INLINE bool operator>=(const Expr& e1, const Expr& e2)
{ return compare(e1, e2) >= 0; }

// Builtin Functions
CORE_INLINE int sign(const Expr& e)
{ return e.sign(); }

CORE_INLINE bool isZero(const Expr& e)
{ return sign(e) == 0; }

CORE_INLINE int compare(const Expr& e1, const Expr& e2)
{ return e1.rep == e2.rep ? 0 : sign(e1 - e2); }

CORE_INLINE BigInt floor(const Expr& e)
{ Expr tmp; return floor(e, tmp); }

CORE_INLINE BigInt ceil(const Expr& e)
{ return -floor(-e); }

CORE_INLINE Expr power(const Expr& e, unsigned long p)
{ return pow(e, p); }

CORE_INLINE Expr abs(const Expr& x)
{ return (sign(x)>=0) ? (x) : (-x); }

CORE_INLINE Expr fabs(const Expr& x)
{ return abs(x); }

// I/O Stream
CORE_INLINE std::ostream& operator<<(std::ostream& o, const Expr& e)
{ o << *(e.rep); return o; }

CORE_INLINE std::istream& operator>>(std::istream& i, Expr& e) {
  Real rVal;
  i >> rVal; 		// precision is = defInputDigits
  if (i)  e = rVal;	// only assign when reading is successful.
  return i;
}

// Deprecated Functions
CORE_INLINE int Expr::getSign() const
{ return rep->getSign(); }

CORE_INLINE int Expr::intValue() const
{ return toInt(); }

CORE_INLINE long Expr::longValue() const
{ return toLong(); }

CORE_INLINE float Expr::floatValue() const
{ return toFloat(); }

CORE_INLINE double Expr::doubleValue() const 
{ return toDouble(); }

CORE_INLINE int ExprRep::getExactSign() { 
  if (!nodeInfo) initNodeInfo();

  if (!flagsComputed()) {
    degreeBound();	
#ifdef DEBUG
    dagSize();
    fullClearFlag();
#endif
    computeExactFlags();
  } 
  return sign();
}

// Chee, 7/17/02: degreeBound() function is now
// taken out of "computeExactFlags()
CORE_INLINE int ExprRep::getSign() { 
  if (ffVal.isOK())
    return ffVal.sign();
  else
    return getExactSign(); 
}   

// you need force to approximate before call this!!
CORE_INLINE BigFloat ExprRep::getBigFloat()
{ return getAppValue().getBigFloat(); }

#endif
