/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: Real.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

//  infty and tiny
const long halfIntMax = INT_MAX /2;
const long halfIntMin = INT_MIN /2;

const long halfLongMax = LONG_MAX /2;
const long halfLongMin = LONG_MIN /2;

////////////////////////////////////////////////////////////////////
// class RealLong
template<>
CORE_INLINE Real RealLong::approx(const extLong& r, const extLong& a) const
{ BigFloat x; x.approx(BigInt(ker), r, a); return Real(x); }

//  unary minus
template<>
CORE_INLINE Real RealLong::operator -() const {
  if (ker < - LONG_MAX)
    return Real(- BigInt(ker));
  else
    return Real(- ker);
}

//  addition
template<>
CORE_INLINE Real RealLong::operator +(const Real& x) const
{ return x.getRep()->addLong(*this); }

template<>
CORE_INLINE Real RealLong::addLong(const RealLong& l) const {
  if (l.ker > halfLongMax && ker > halfLongMax
      || l.ker < halfLongMin && ker < halfLongMin)
    return Real(BigInt(l.ker) + BigInt(ker));
  else
    return Real(l.ker + ker);
}

template<>
CORE_INLINE Real RealLong::addDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealLong::addBigInt(const RealBigInt& I) const
{ return Real(I.ker + BigInt(ker)); }

template<>
CORE_INLINE Real RealLong::addBigFloat(const RealBigFloat& B) const
{ return Real(B.ker + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealLong::addBigRat(const RealBigRat& R) const
{ return Real(R.ker + BigRat(ker)); }

//  subtraction
template<>
CORE_INLINE Real RealLong::operator -(const Real& x) const
{ return x.getRep()->subLong(*this); }

template<>
CORE_INLINE Real RealLong::subLong(const RealLong& l) const {
  if (l.ker < halfLongMin && ker > halfLongMax
      || l.ker > halfLongMax && ker < halfLongMin)
    return Real(BigInt(l.ker) - BigInt(ker));
  else
    return Real(l.ker - ker);
}

template<>
CORE_INLINE Real RealLong::subDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealLong::subBigInt(const RealBigInt& I) const
{ return Real(I.ker - BigInt(ker)); }

template<>
CORE_INLINE Real RealLong::subBigFloat(const RealBigFloat& B) const
{ return Real(B.ker - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealLong::subBigRat(const RealBigRat& R) const
{ return Real(R.ker - BigRat(ker)); }

//  multiplication
template<>
CORE_INLINE Real RealLong::operator *(const Real& x) const
{ return x.getRep()->mulLong(*this); }

template<>
CORE_INLINE Real RealLong::mulLong(const RealLong& l) const {
  if (flrLg(l.ker) + flrLg(ker) < (int)(LONG_BIT - 2))
    return Real(l.ker * ker);
  else
    return Real(BigInt(l.ker) * BigInt(ker));
}

template<>
CORE_INLINE Real RealLong::mulDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealLong::mulBigInt(const RealBigInt& I) const
{ return Real(I.ker * BigInt(ker)); }

template<>
CORE_INLINE Real RealLong::mulBigFloat(const RealBigFloat& B) const
{ return Real(B.ker * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealLong::mulBigRat(const RealBigRat& R) const
{ return Real(BigRat(R.ker.numerator() * ker, R.ker.denominator())); }

//  division
template<>
CORE_INLINE Real RealLong::div(const Real& x, const extLong& r) const
{ return x.getRep()->divLong(*this, r); }

template<>
CORE_INLINE Real RealLong::divLong(const RealLong& l, const extLong&) const
{ return Real(BigRat(l.ker, ker)); }

template<>
CORE_INLINE Real RealLong::divDouble(const RealDouble& d, const extLong& r) const
{ return Real(BigFloat(d.ker).div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealLong::divBigInt(const RealBigInt& I, const extLong&) const
{ return Real(BigRat(I.ker, BigInt(ker))); }

template<>
CORE_INLINE Real RealLong::divBigFloat(const RealBigFloat& B, const extLong& r) const
{ return Real(B.ker.div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealLong::divBigRat(const RealBigRat& R, const extLong&) const
{ return Real(BigRat(R.ker.numerator(), R.ker.denominator() * ker)); }

//  squareroot
template<>
CORE_INLINE Real RealLong::sqrt(const extLong& a) const
{ return Real(BigFloat(ker).sqrt(a)); }

//  squareroot with initial approximation A
template<>
CORE_INLINE Real RealLong::sqrt(const extLong& a, const BigFloat& A) const
{ return Real(BigFloat(ker).sqrt(a, A) ); }

//  equality
template<>
CORE_INLINE bool RealLong::operator ==(const Real& x) const
{ return x.getRep()->eqlLong(*this); }

template<>
CORE_INLINE bool RealLong::eqlLong(const RealLong& l) const
{ return l.ker == ker; }

template<>
CORE_INLINE bool RealLong::eqlDouble(const RealDouble& d) const
{ return d.ker == double(ker); }

template<>
CORE_INLINE bool RealLong::eqlBigInt(const RealBigInt& I) const
{ return I.ker == BigInt(ker); }

template<>
CORE_INLINE bool RealLong::eqlBigFloat(const RealBigFloat& B) const
{ return B.ker == BigFloat(ker); }

template<>
CORE_INLINE bool RealLong::eqlBigRat(const RealBigRat& R) const
{ return R.ker == BigRat(ker); }

//  smaller-than
template<>
CORE_INLINE bool RealLong::operator <(const Real& x) const
{ return x.getRep()->grtLong(*this); }

template<>
CORE_INLINE bool RealLong::grtLong(const RealLong& l) const
{ return l.ker < ker; }

template<>
CORE_INLINE bool RealLong::grtDouble(const RealDouble& d) const
{ return d.ker < double(ker); }

template<>
CORE_INLINE bool RealLong::grtBigInt(const RealBigInt& I) const
{ return I.ker < BigInt(ker); }

template<>
CORE_INLINE bool RealLong::grtBigFloat(const RealBigFloat& B) const
{ return B.ker < BigFloat(ker); }

template<>
CORE_INLINE bool RealLong::grtBigRat(const RealBigRat& R) const
{ return R.ker < BigRat(ker); }

//  builtin functions
template<>
CORE_INLINE bool RealLong::isExact() const
{ return true; }

template<>
CORE_INLINE
void RealLong::ULV_E(extLong &up, extLong &lp, extLong &v2p,
                     extLong &v2m, extLong &v5p, extLong &v5m) const
{
  // TODO : extract the power of 5.
  up = lp = v2p = v2m = v5p = v5m = 0;
  if (ker == 0)
    return;

  // Extract the power of 2.
  unsigned long exp = 0;
  unsigned long tmp_ker = ker;
  while ((tmp_ker&1) != 0) {
    tmp_ker = tmp_ker/2;
    ++exp;
  }
  up = clLg(tmp_ker);
  lp = 0;
  v2p = exp;
}

template<>
CORE_INLINE extLong RealLong::flrLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE extLong RealLong::clLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE bool RealLong::isZeroIn() const
{ return ker == 0; }

template<>
CORE_INLINE unsigned long RealLong::degree() const
{ return 1; }

template<>
CORE_INLINE unsigned long RealLong::length() const
{ return clLg(1+ core_abs(ker)); }	// length is (log_2(1+ker^2)) /2.

template<>
CORE_INLINE unsigned long RealLong::height() const
{ return clLg(core_max(1L, core_abs(ker))); }	// height is max{1, |ker|}

template<>
CORE_INLINE int RealLong::sgn() const
{ return ker > 0 ? 1 : (!ker ? 0 : - 1); }

// 31/05/2003 Joaquin Grech
template<>
CORE_INLINE std::string RealLong::toString(long prec, bool sci) const
{   
	std::stringstream st;
	st << ker;
	return st.str();
}

//////////////////////////////////////////////////////////////////
// class RealDouble
template<>
CORE_INLINE Real RealDouble::approx(const extLong& r, const extLong& a) const
{ BigFloat x; x.approx(BigRat(ker), r, a); return Real(x); }

//  unary minus
template<>
CORE_INLINE Real RealDouble::operator -() const
{ return Real(- ker); }

//  addition
template<>
CORE_INLINE Real RealDouble::operator +(const Real& x) const
{ return x.getRep()->addDouble(*this); }

template<>
CORE_INLINE Real RealDouble::addLong(const RealLong& l) const
{ return Real(BigFloat(l.ker) + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::addDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::addBigInt(const RealBigInt& I) const
{ return Real(BigFloat(I.ker) + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::addBigFloat(const RealBigFloat& B) const
{ return Real(B.ker + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::addBigRat(const RealBigRat& R) const
{ return Real(R.ker + BigRat(ker)); }

//  subtraction
template<>
CORE_INLINE Real RealDouble::operator -(const Real& x) const
{ return x.getRep()->subDouble(*this); }

template<>
CORE_INLINE Real RealDouble::subLong(const RealLong& l) const
{ return Real(BigFloat(l.ker) - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::subDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::subBigInt(const RealBigInt& I) const
{ return Real(BigFloat(I.ker) - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::subBigFloat(const RealBigFloat& B) const
{ return Real(B.ker - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::subBigRat(const RealBigRat& R) const
{ return Real(R.ker - BigRat(ker)); }

//  multiplication
template<>
CORE_INLINE Real RealDouble::operator *(const Real& x) const
{ return x.getRep()->mulDouble(*this); }

template<>
CORE_INLINE Real RealDouble::mulLong(const RealLong& l) const
{ return Real(BigFloat(l.ker) * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::mulDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::mulBigInt(const RealBigInt& I) const
{ return Real(BigFloat(I.ker) * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::mulBigFloat(const RealBigFloat& B) const
{ return Real(B.ker * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealDouble::mulBigRat(const RealBigRat& R) const
{ return Real(R.ker * BigRat(ker)); }

//  division
template<>
CORE_INLINE Real RealDouble::div(const Real& x, const extLong& r) const
{ return x.getRep()->divDouble(*this, r); }

template<>
CORE_INLINE Real RealDouble::divLong(const RealLong& l, const extLong& r) const
{ return Real(BigFloat(l.ker).div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealDouble::divDouble(const RealDouble& d, const extLong& r) const
{ return Real(BigFloat(d.ker).div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealDouble::divBigInt(const RealBigInt& I, const extLong& r) const
{ return Real(BigFloat(I.ker).div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealDouble::divBigFloat(const RealBigFloat& B, const extLong& r) const
{ return Real(B.ker.div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealDouble::divBigRat(const RealBigRat& R, const extLong&) const
{ return Real(R.ker / BigRat(ker)); }

//  squareroot
template<>
CORE_INLINE Real RealDouble::sqrt(const extLong& a) const
{ return Real(BigFloat(ker).sqrt(a)); }

//  squareroot, with initial approximation A
template<>
CORE_INLINE Real RealDouble::sqrt(const extLong& a, const BigFloat & A) const
{ return Real(BigFloat(ker).sqrt(a, A)); }

//  equality
template<>
CORE_INLINE bool RealDouble::operator ==(const Real& x) const
{ return x.getRep()->eqlDouble(*this); }

template<>
CORE_INLINE bool RealDouble::eqlLong(const RealLong& l) const
{ return double(l.ker) == ker; }

template<>
CORE_INLINE bool RealDouble::eqlDouble(const RealDouble& d) const
{ return d.ker == ker; }

template<>
CORE_INLINE bool RealDouble::eqlBigInt(const RealBigInt& I) const
{ return BigFloat(I.ker) == BigFloat(ker); }

template<>
CORE_INLINE bool RealDouble::eqlBigFloat(const RealBigFloat& B) const
{ return B.ker == BigFloat(ker); }

template<>
CORE_INLINE bool RealDouble::eqlBigRat(const RealBigRat& R) const
{ return R.ker == BigRat(ker); }

//  smaller-than
template<>
CORE_INLINE bool RealDouble::operator <(const Real& x) const
{ return x.getRep()->grtDouble(*this); }

template<>
CORE_INLINE bool RealDouble::grtLong(const RealLong& l) const
{ return double(l.ker) < ker; }

template<>
CORE_INLINE bool RealDouble::grtDouble(const RealDouble& d) const
{ return d.ker < ker; }

template<>
CORE_INLINE bool RealDouble::grtBigInt(const RealBigInt& I) const
{ return BigFloat(I.ker) < BigFloat(ker); }

template<>
CORE_INLINE bool RealDouble::grtBigFloat(const RealBigFloat& B) const
{ return B.ker < BigFloat(ker); }

template<>
CORE_INLINE bool RealDouble::grtBigRat(const RealBigRat& R) const
{ return R.ker < BigRat(ker); }

//  builtin functions
template<>
CORE_INLINE bool RealDouble::isExact() const
{ return true; }

template<>
CORE_INLINE
void RealDouble::ULV_E(extLong &up, extLong &lp, extLong &v2p,
		       extLong &v2m, extLong &v5p, extLong &v5m) const
{
  // TODO : can probably be made faster using frexp() or such.
  // TODO : extract the power of 5.
  BigRat R = BigRat(ker);
  up  = ceilLg(R.numerator());
  v2m = ceilLg(R.denominator());
  lp = v2p = v5m = v5p = 0;
}

template<>
CORE_INLINE extLong RealDouble::flrLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE extLong RealDouble::clLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE bool RealDouble::isZeroIn() const
{ return ker == 0.0; }

template<>
CORE_INLINE unsigned long RealDouble::degree() const
{ return 1; }

template<>
CORE_INLINE unsigned long RealDouble::length() const {
  BigRat R  = BigRat(ker); 
  long ln = 1 + ceilLg(R.numerator());
  long ld = 1 + ceilLg(R.denominator());
  return (ln>ld) ? ln : ld; ///< an upper bound on log_2(sqrt(num^2+den^2))
}

template<>
CORE_INLINE unsigned long RealDouble::height() const {
  BigRat R  = BigRat(ker);
  long ln = ceilLg(R.numerator());
  long ld = ceilLg(R.denominator());
  return (ln>ld) ? ln : ld; ///< an upper bound on log_2(max(|num|, |den|))
}

template<>
CORE_INLINE int RealDouble::sgn() const
{ return ker > 0.0 ? 1 : (ker == 0.0 ? 0 : - 1); }

// 31/05/2003 Joaquin Grech
// Bug? shall we need to use "sci" to control output
// Zilin Du
template<>
CORE_INLINE std::string RealDouble::toString(long prec, bool sci) const
{ 
	std::stringstream st;
	st << ker;
	return st.str();
}

///////////////////////////////////////////////////////////////////
// class RealBigInt
template<>
CORE_INLINE Real RealBigInt::approx(const extLong& r, const extLong& a) const
{ BigFloat x; x.approx(ker, r, a); return Real(x); }

//  unary minus
template<>
CORE_INLINE Real RealBigInt::operator -() const
{ return Real(- ker); }

//  addition
template<>
CORE_INLINE Real RealBigInt::operator +(const Real& x) const
{ return x.getRep()->addBigInt(*this); }

template<>
CORE_INLINE Real RealBigInt::addLong(const RealLong& l) const
{ return Real(BigInt(l.ker) + ker); }

template<>
CORE_INLINE Real RealBigInt::addDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealBigInt::addBigInt(const RealBigInt& I) const
{ return Real(I.ker + ker); }

template<>
CORE_INLINE Real RealBigInt::addBigFloat(const RealBigFloat& B) const
{ return Real(B.ker + (BigFloat(ker))); }

template<>
CORE_INLINE Real RealBigInt::addBigRat(const RealBigRat& R) const
{ return Real(R.ker + BigRat(ker)); }

//  subtraction
template<>
CORE_INLINE Real RealBigInt::operator -(const Real& x) const
{ return x.getRep()->subBigInt(*this); }

template<>
CORE_INLINE Real RealBigInt::subLong(const RealLong& l) const
{ return Real(BigInt(l.ker) - ker); }

template<>
CORE_INLINE Real RealBigInt::subDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealBigInt::subBigInt(const RealBigInt& I) const
{ return Real(I.ker - ker); }

template<>
CORE_INLINE Real RealBigInt::subBigFloat(const RealBigFloat& B) const
{ return Real(B.ker - (BigFloat(ker))); }

template<>
CORE_INLINE Real RealBigInt::subBigRat(const RealBigRat& R) const
{ return Real(R.ker - BigRat(ker)); }

//  multiplication
template<>
CORE_INLINE Real RealBigInt::operator *(const Real& x) const
{ return x.getRep()->mulBigInt(*this); }

template<>
CORE_INLINE Real RealBigInt::mulLong(const RealLong& l) const
{ return Real(BigInt(l.ker) * ker); }

template<>
CORE_INLINE Real RealBigInt::mulDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealBigInt::mulBigInt(const RealBigInt& I) const
{ return Real(I.ker * ker); }

template<>
CORE_INLINE Real RealBigInt::mulBigFloat(const RealBigFloat& B) const
{ return Real(B.ker * (BigFloat(ker))); }

template<>
CORE_INLINE Real RealBigInt::mulBigRat(const RealBigRat& R) const
{ return Real(BigRat(R.ker.numerator() * ker, R.ker.denominator())); }

//  division
template<>
CORE_INLINE Real RealBigInt::div(const Real& x, const extLong& r) const
{ return x.getRep()->divBigInt(*this, r); }

template<>
CORE_INLINE Real RealBigInt::divLong(const RealLong& l, const extLong&) const
{ return Real(BigRat(BigInt(l.ker), ker)); }

template<>
CORE_INLINE Real RealBigInt::divDouble(const RealDouble& d, const extLong& r) const
{ return Real(BigFloat(d.ker).div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealBigInt::divBigInt(const RealBigInt& I, const extLong&) const
{ return Real(BigRat(I.ker, ker)); }

template<>
CORE_INLINE Real RealBigInt::divBigFloat(const RealBigFloat& B, const extLong& r) const
{ return Real(B.ker.div(BigFloat(ker), r)); }

template<>
CORE_INLINE Real RealBigInt::divBigRat(const RealBigRat& R, const extLong&) const
{ return Real(BigRat(R.ker.numerator(), R.ker.denominator() * ker)); }

//  squareroot
template<>
CORE_INLINE Real RealBigInt::sqrt(const extLong& a) const
{ return Real(BigFloat(ker).sqrt(a)); }

//  squareroot, with initial approximation A
template<>
CORE_INLINE Real RealBigInt::sqrt(const extLong& a, const BigFloat& A) const
{ return Real(BigFloat(ker).sqrt(a, A)); }

//  equality
template<>
CORE_INLINE bool RealBigInt::operator ==(const Real& x) const
{ return x.getRep()->eqlBigInt(*this); }

template<>
CORE_INLINE bool RealBigInt::eqlLong(const RealLong& l) const
{ return BigInt(l.ker) == ker; }

template<>
CORE_INLINE bool RealBigInt::eqlDouble(const RealDouble& d) const
{ return BigFloat(d.ker) == BigFloat(ker); }

template<>
CORE_INLINE bool RealBigInt::eqlBigInt(const RealBigInt& I) const
{ return I.ker == ker; }

template<>
CORE_INLINE bool RealBigInt::eqlBigFloat(const RealBigFloat& B) const
{ return B.ker == BigFloat(ker); } 

template<>
CORE_INLINE bool RealBigInt::eqlBigRat(const RealBigRat& R) const
{ return R.ker == BigRat(ker); }

//  smaller-than
template<>
CORE_INLINE bool RealBigInt::operator <(const Real& x) const
{ return x.getRep()->grtBigInt(*this); }

template<>
CORE_INLINE bool RealBigInt::grtLong(const RealLong& l) const
{ return BigInt(l.ker) < ker; }

template<>
CORE_INLINE bool RealBigInt::grtDouble(const RealDouble& d) const
{ return BigFloat(d.ker) < BigFloat(ker); }

template<>
CORE_INLINE bool RealBigInt::grtBigInt(const RealBigInt& I) const
{ return I.ker < ker; }

template<>
CORE_INLINE bool RealBigInt::grtBigFloat(const RealBigFloat& B) const
{ return B.ker < BigFloat(ker); }

template<>
CORE_INLINE bool RealBigInt::grtBigRat(const RealBigRat& R) const
{ return R.ker < BigRat(ker); }

//  builtin functions
template<>
CORE_INLINE bool RealBigInt::isExact() const
{ return true; }

template<>
CORE_INLINE
void RealBigInt::ULV_E(extLong &up, extLong &lp, extLong &v2p,
		       extLong &v2m, extLong &v5p, extLong &v5m) const
{
  up = lp = v2p = v2m = v5p = v5m = 0;
  if (ker == 0)
    return;

  // Extract power of 5.
  int exp5;
  BigInt remainder5;
  ker.getKaryExpo(remainder5, exp5, 5);
  v5p = exp5;
  // Extract power of 2.
  int exp2 = remainder5.getBinExpo();
  up = ceilLg(remainder5) - exp2;
  v2p = exp2;
}

template<>
CORE_INLINE extLong RealBigInt::flrLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE extLong RealBigInt::clLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE bool RealBigInt::isZeroIn() const
{ return sign(ker) == 0; }

template<>
CORE_INLINE unsigned long RealBigInt::degree() const
{ return 1; }

template<>
CORE_INLINE unsigned long RealBigInt::length() const
{ return ceilLg(1 + core_abs(ker)); }

template<>
CORE_INLINE unsigned long RealBigInt::height() const
{ return ceilLg( core_max(BigInt(1), core_abs(ker) )); }

template<>
CORE_INLINE int RealBigInt::sgn() const
{ return sign(ker); }

// 31/05/2003 Joaquin Grech
template<>
CORE_INLINE std::string RealBigInt::toString(long prec, bool sci) const
{ return ker.toString(); }

/////////////////////////////////////////////////////////////
// class RealBigFloat

// Chee, 8/2/01:  restored computation of x.approx(ker, r, a) 
template<>
CORE_INLINE Real RealBigFloat::approx(const extLong& r, const extLong& a) const
{ BigFloat x; x.approx(ker, r, a); return Real(x); }

//  unary minus
  
template<>
CORE_INLINE Real RealBigFloat::operator -() const
{ return Real(- ker); }

//  addition
template<>
CORE_INLINE Real RealBigFloat::operator +(const Real& x) const
{ return x.getRep()->addBigFloat(*this); }

template<>
CORE_INLINE Real RealBigFloat::addLong(const RealLong& l) const
{ return Real(BigFloat(l.ker) + (ker)); }

template<>
CORE_INLINE Real RealBigFloat::addDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) + (ker)); }

template<>
CORE_INLINE Real RealBigFloat::addBigInt(const RealBigInt& I) const
{ return Real(BigFloat(I.ker) + (ker)); }

template<>
CORE_INLINE Real RealBigFloat::addBigFloat(const RealBigFloat& B) const
{ return Real(B.ker + (ker)); }

template<>
CORE_INLINE Real RealBigFloat::addBigRat(const RealBigRat& R) const {
  if (ker.isExact()) {
    // Chen: the following cause un-wanted implicit conversion from
    // BigFloat to BigRat
    //    return Real(R.ker + BigRat(ker));
    return Real(R.ker + ker.BigRatValue());
  } else {
    BigFloat x;
    x.approx(R.ker, CORE_posInfty, - ker.flrLgErr());
    return Real(x + (ker));
  }
}

//  subtraction

template<>
CORE_INLINE Real RealBigFloat::operator -(const Real& x) const
{ return x.getRep()->subBigFloat(*this); }

template<>
CORE_INLINE Real RealBigFloat::subLong(const RealLong& l) const
{ return Real(BigFloat(l.ker) - (ker)); }

template<>
CORE_INLINE Real RealBigFloat::subDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) - (ker)); }

template<>
CORE_INLINE Real RealBigFloat::subBigInt(const RealBigInt& I) const
{ return Real(BigFloat(I.ker) - (ker)); }

template<>
CORE_INLINE Real RealBigFloat::subBigFloat(const RealBigFloat& B) const
{ return Real(B.ker - (ker)); }

template<>
CORE_INLINE Real RealBigFloat::subBigRat(const RealBigRat& R) const {
  if (ker.isExact())
    //    return Real(R.ker - BigRat(ker));
    return Real(R.ker - ker.BigRatValue());
  else {
      BigFloat x;
      x.approx(R.ker, CORE_posInfty, - ker.flrLgErr());
      return Real(x - (ker));
  }
}

//  multiplication
template<>
CORE_INLINE Real RealBigFloat::operator *(const Real& x) const
{ return x.getRep()->mulBigFloat(*this); }

template<>
CORE_INLINE Real RealBigFloat::mulLong(const RealLong& l) const
{ return Real(BigFloat(l.ker) * (ker)); }

template<>
CORE_INLINE Real RealBigFloat::mulDouble(const RealDouble& d) const
{ return Real(BigFloat(d.ker) * (ker)); }

template<>
CORE_INLINE Real RealBigFloat::mulBigInt(const RealBigInt& I) const
{ return Real(BigFloat(I.ker) *(ker)); }

template<>
CORE_INLINE Real RealBigFloat::mulBigFloat(const RealBigFloat& B) const
{ return Real(B.ker *(ker)); }

template<>
CORE_INLINE Real RealBigFloat::mulBigRat(const RealBigRat& R) const {
  if (ker.isExact())
    return Real(R.ker * ker.BigRatValue());
  else {
    BigFloat x;
    x.approx(R.ker, extLong(ker.MSB() - ker.flrLgErr()) + 1, CORE_posInfty);
    return Real(x * (ker));
  }
}

//  division
template<>
CORE_INLINE Real RealBigFloat::div(const Real& x, const extLong& r) const
{ return x.getRep()->divBigFloat(*this, r); }

template<>
CORE_INLINE Real RealBigFloat::divLong(const RealLong& l, const extLong& r) const
{ return Real(BigFloat(l.ker).div(ker, r)); }

template<>
CORE_INLINE Real RealBigFloat::divDouble(const RealDouble& d, const extLong& r) const
{ return Real(BigFloat(d.ker).div(ker, r)); }

template<>
CORE_INLINE Real RealBigFloat::divBigInt(const RealBigInt& I, const extLong& r) const
{ return Real(BigFloat(I.ker).div(ker, r)); }

template<>
CORE_INLINE Real RealBigFloat::divBigFloat(const RealBigFloat& B,
					const extLong& r) const
{ return Real(BigFloat(B.ker).div(ker, r)); }

template<>
CORE_INLINE Real RealBigFloat::divBigRat(const RealBigRat& R,
					const extLong& r) const
{
  if (ker.isExact())
    return Real(R.ker / ker.BigRatValue());
  else {
    BigFloat x;
    x.approx(R.ker, extLong(ker.MSB() - ker.flrLgErr()) + 1, CORE_posInfty);
    return Real(x.div(ker, r));
  }
}

//  squareroot, computed to absolute precision a
template<>
CORE_INLINE Real RealBigFloat::sqrt(const extLong& a) const
{ return Real(ker.sqrt(a)); }

//  squareroot, computed to absolute precision a, but with initial approximation A
template<>
CORE_INLINE Real RealBigFloat::sqrt(const extLong& a, const BigFloat & A) const
{ return Real(ker.sqrt(a, A)); }

//  equality
template<>
CORE_INLINE bool RealBigFloat::operator ==(const Real& x) const
{ return x.getRep()->eqlBigFloat(*this); }

template<>
CORE_INLINE bool RealBigFloat::eqlLong(const RealLong& l) const
{ return BigFloat(l.ker) == ker; }

template<>
CORE_INLINE bool RealBigFloat::eqlDouble(const RealDouble& d) const
{ return BigFloat(d.ker) == ker; }

template<>
CORE_INLINE bool RealBigFloat::eqlBigInt(const RealBigInt& I) const
{ return BigFloat(I.ker) == ker; }

template<>
CORE_INLINE bool RealBigFloat::eqlBigFloat(const RealBigFloat& B) const
{ return B.ker == ker; }

template<>
CORE_INLINE bool RealBigFloat::eqlBigRat(const RealBigRat& R) const
{ return R.ker == ker.BigRatValue(); }

//  smaller-than
template<>
CORE_INLINE bool RealBigFloat::operator <(const Real& x) const
{ return x.getRep()->grtBigFloat(*this); }

template<>
CORE_INLINE bool RealBigFloat::grtLong(const RealLong& l) const
{ return BigFloat(l.ker) < ker; }

template<>
CORE_INLINE bool RealBigFloat::grtDouble(const RealDouble& d) const
{ return BigFloat(d.ker) < ker; }

template<>
CORE_INLINE bool RealBigFloat::grtBigInt(const RealBigInt& I) const
{ return BigFloat(I.ker) < ker; }

template<>
CORE_INLINE bool RealBigFloat::grtBigFloat(const RealBigFloat& B) const
{ return B.ker < ker; }

template<>
CORE_INLINE bool RealBigFloat::grtBigRat(const RealBigRat& R) const
{ return R.ker < ker.BigRatValue(); }

//  builtin functions
template<>
CORE_INLINE bool RealBigFloat::isExact() const
{ return ker.isExact(); }

template<>
CORE_INLINE
void RealBigFloat::ULV_E(extLong &up, extLong &lp, extLong &v2p,
		         extLong &v2m, extLong &v5p, extLong &v5m) const
{
  // TODO : extract power of 5.
  up = lp = v2p = v2m = v5p = v5m = 0;
  BigRat R = ker.BigRatValue();
  up  = ceilLg(R.numerator());
  v2m = ceilLg(R.denominator());
}

template<>
CORE_INLINE extLong RealBigFloat::flrLgErr() const
{ return ker.flrLgErr(); }

template<>
CORE_INLINE extLong RealBigFloat::clLgErr() const
{ return ker.clLgErr(); }

template<>
CORE_INLINE bool RealBigFloat::isZeroIn() const
{ return ker.isZeroIn(); }

template<>
CORE_INLINE unsigned long RealBigFloat::degree() const
{ return 1; }

template<>
CORE_INLINE unsigned long RealBigFloat::length() const {
  // Chen Li: A bug fixed.
  // The statement in the older version with the bug was:
  //   BigRat R  = BigRat(ker);
  // The BigRat(BigFloat) actually is a
  // conversion operator (defined in BigFloat.h), _NOT_
  // an ordinary class constructor! The C++ language
  // specify that an intialization is not an assignment
  // but a constructor operation!
  // Considering that BigRat(BigFloat) is a conversion
  // operator not really a constructor. The programmer's
  // intent is obvious to do an assignment.
  // However, the g++ seems to be confused by the above
  // initialization.
  BigRat R  = ker.BigRatValue();
  long   ln = 1 + ceilLg(R.numerator());
  long   ld = 1 + ceilLg(R.denominator());
  return ( ln > ld ) ? ln : ld; 
}

template<>
CORE_INLINE unsigned long RealBigFloat::height() const {
  // Chen Li: A bug fixed. The old statement with the bug was:
  //   BigRat R  = BigRat(ker);
  // Detailed reasons see above (in RealBigFloat::length()!
  BigRat R  = ker.BigRatValue();
  long     ln = ceilLg(R.numerator());
  long     ld = ceilLg(R.denominator());
  return   ( ln > ld ) ? ln : ld; 
}

template<>
CORE_INLINE int RealBigFloat::sgn() const
{ return ker.sign(); }

// 31/05/2003 Joaquin Grech
template<>
CORE_INLINE std::string RealBigFloat::toString(long prec, bool sci) const
{ return ker.toString(prec,sci); }

//////////////////////////////////////////////////////////////////
// class RealBigRat

template<>
CORE_INLINE Real RealBigRat::approx(const extLong& r, const extLong& a) const
{ BigFloat x; x.approx(ker, r, a); return Real(x); }

//  unary minus
template<>
CORE_INLINE Real RealBigRat::operator -() const
{ return Real(- ker); }

//  addition
template<>
CORE_INLINE Real RealBigRat::operator +(const Real& x) const
{ return x.getRep()->addBigRat(*this); }

template<>
CORE_INLINE Real RealBigRat::addLong(const RealLong& l) const
{ return Real(BigRat(l.ker) + ker); }

template<>
CORE_INLINE Real RealBigRat::addDouble(const RealDouble& d) const
{ return Real(BigRat(d.ker) + ker); }

template<>
CORE_INLINE Real RealBigRat::addBigInt(const RealBigInt& I) const
{ return Real(BigRat(I.ker) + ker); }

template<>
CORE_INLINE Real RealBigRat::addBigFloat(const RealBigFloat& B) const {
  if (B.ker.isExact())
    return Real(B.ker.BigRatValue() + ker);
  else {
    BigFloat x;
    x.approx(ker, CORE_posInfty, - B.ker.flrLgErr());
    return Real(B.ker + (x));
  }
}

template<>
CORE_INLINE Real RealBigRat::addBigRat(const RealBigRat& R) const
{ return Real(R.ker + ker); }

//  subtraction
template<>
CORE_INLINE Real RealBigRat::operator -(const Real& x) const
{ return x.getRep()->subBigRat(*this); }

template<>
CORE_INLINE Real RealBigRat::subLong(const RealLong& l) const
{ return Real(BigRat(l.ker) - ker); }

template<>
CORE_INLINE Real RealBigRat::subDouble(const RealDouble& d) const
{ return Real(BigRat(d.ker) - ker); }

template<>
CORE_INLINE Real RealBigRat::subBigInt(const RealBigInt& I) const
{ return Real(BigRat(I.ker) - ker); }

template<>
CORE_INLINE Real RealBigRat::subBigFloat(const RealBigFloat& B) const {
  if (B.ker.isExact())
    return Real(B.ker.BigRatValue() - ker);
  else {
    BigFloat x;
    x.approx(ker, CORE_posInfty, - B.ker.flrLgErr());
    return Real(B.ker - (x));
  }
}

template<>
CORE_INLINE Real RealBigRat::subBigRat(const RealBigRat& R) const
{ return Real(R.ker - ker); }

//  multiplication
template<>
CORE_INLINE Real RealBigRat::operator *(const Real& x) const
{ return x.getRep()->mulBigRat(*this); }

template<>
CORE_INLINE Real RealBigRat::mulLong(const RealLong& l) const
{ return Real(BigRat(l.ker * ker.numerator(), ker.denominator())); }

template<>
CORE_INLINE Real RealBigRat::mulDouble(const RealDouble& d) const
{ return Real(BigRat(d.ker) * ker); }

template<>
CORE_INLINE Real RealBigRat::mulBigInt(const RealBigInt& I) const
{ return Real(BigRat(I.ker * ker.numerator(), ker.denominator())); }

template<>
CORE_INLINE Real RealBigRat::mulBigFloat(const RealBigFloat& B) const {
  if (B.ker.isExact())
    return Real(B.ker.BigRatValue() * ker);
  else {
    BigFloat x;
    x.approx(ker, extLong(B.ker.MSB() - B.ker.flrLgErr()) + 1, CORE_posInfty);
    return Real(B.ker *(x));
  }
}

template<>
CORE_INLINE Real RealBigRat::mulBigRat(const RealBigRat& R) const
{ return Real(R.ker * ker); }

//  division
template<>
CORE_INLINE Real RealBigRat::div(const Real& x, const extLong& r) const
{ return x.getRep()->divBigRat(*this, r); }

template<>
CORE_INLINE Real RealBigRat::divLong(const RealLong& l, const extLong&) const
{ return Real(BigRat(l.ker * ker.denominator(), ker.numerator())); }

template<>
CORE_INLINE Real RealBigRat::divDouble(const RealDouble& d, const extLong&) const
{ return Real(BigRat(d.ker) / ker); }

template<>
CORE_INLINE Real RealBigRat::divBigInt(const RealBigInt& I, const extLong&) const
{ return Real(BigRat(I.ker * ker.denominator(), ker.numerator())); }

template<>
CORE_INLINE Real RealBigRat::divBigFloat(const RealBigFloat& B,
					const extLong& r) const
{
  if (B.ker.isExact())
    return Real(B.ker.BigRatValue() / ker);
  else {
    BigFloat x;
    x.approx(ker, extLong(B.ker.MSB() - B.ker.flrLgErr()) + 1, CORE_posInfty);
    return Real(B.ker.div(x, r));
  }
}

template<>
CORE_INLINE Real RealBigRat::divBigRat(const RealBigRat& R,
				      const extLong&) const
{ return Real(R.ker / ker); }

//  squareroot
template<>
CORE_INLINE Real RealBigRat::sqrt(const extLong& a) const{
  BigFloat x; x.approx(ker, CORE_posInfty, 2*defAbsPrec+8);
  return Real(x.sqrt(a));
}
  
// squareroot, computed to absolute precision a, but with initial approximation A
template<>
CORE_INLINE Real RealBigRat::sqrt(const extLong& a, const BigFloat & A) const
{ return Real(BigFloat(ker).sqrt(a,A)); }


//  equality
template<>
CORE_INLINE bool RealBigRat::operator ==(const Real& x) const
{ return x.getRep()->eqlBigRat(*this); }

template<>
CORE_INLINE bool RealBigRat::eqlLong(const RealLong& l) const
{ return BigRat(l.ker) == ker; }

template<>
CORE_INLINE bool RealBigRat::eqlDouble(const RealDouble& d) const
{ return BigRat(d.ker) == ker; }

template<>
CORE_INLINE bool RealBigRat::eqlBigInt(const RealBigInt& I) const
{ return BigRat(I.ker) == ker; }

template<>
CORE_INLINE bool RealBigRat::eqlBigFloat(const RealBigFloat& B) const
{ return B.ker.BigRatValue() == ker; }

template<>
CORE_INLINE bool RealBigRat::eqlBigRat(const RealBigRat& R) const
{ return R.ker == ker; }

//  smaller-than
template<>
CORE_INLINE bool RealBigRat::operator <(const Real& x) const
{ return x.getRep()->grtBigRat(*this); }

template<>
CORE_INLINE bool RealBigRat::grtLong(const RealLong& l) const
{ return BigRat(l.ker) < ker; }

template<>
CORE_INLINE bool RealBigRat::grtDouble(const RealDouble& d) const
{ return BigRat(d.ker) < ker; }

template<>
CORE_INLINE bool RealBigRat::grtBigInt(const RealBigInt& I) const
{ return BigRat(I.ker) < ker; }

template<>
CORE_INLINE bool RealBigRat::grtBigFloat(const RealBigFloat& B) const
{ return B.ker.BigRatValue() < ker; }

template<>
CORE_INLINE bool RealBigRat::grtBigRat(const RealBigRat& R) const
{ return R.ker < ker; }

//  builtin functions
template<>
CORE_INLINE bool RealBigRat::isExact() const
{ return true; }

template<>
CORE_INLINE
void RealBigRat::ULV_E(extLong &up, extLong &lp, extLong &v2p,
		       extLong &v2m, extLong &v5p, extLong &v5m) const
{
  up = lp = v2p = v2m = v5p = v5m = 0;
  if (ker == 0)
    return;

  // Extract power of 5.
  int exp5;
  BigInt num5, den5;
  ker.numerator().getKaryExpo(num5, exp5, 5);
  if (exp5 != 0) {
    v5p = exp5;
    den5 = ker.denominator();
  } else {
    ker.denominator().getKaryExpo(den5, exp5, 5);
    v5m = exp5;
  }

  // Now we work with num5/den5.
  int exp2 = num5.getBinExpo();
  if (exp2 != 0) {
    v2p = exp2;
  } else {
    exp2 = den5.getBinExpo();
    v2m = exp2;
  }

  up = ceilLg(num5) - v2p;
  lp = ceilLg(den5) - v2m;
}

template<>
CORE_INLINE extLong RealBigRat::flrLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE extLong RealBigRat::clLgErr() const
{ return extLong(CORE_negInfty); }

template<>
CORE_INLINE bool RealBigRat::isZeroIn() const
{ return ker.sign() == 0; }

template<>
CORE_INLINE unsigned long RealBigRat::degree() const
{ return 1; }

template<>
CORE_INLINE unsigned long RealBigRat::length() const {
  long ln = 1 + ceilLg(ker.numerator());
  long ld = 1 + ceilLg(ker.denominator());
  return ( ln > ld ) ? ln : ld; 
}

template<>
CORE_INLINE unsigned long RealBigRat::height() const {
  long ln = ceilLg(ker.numerator());
  long ld = ceilLg(ker.denominator());
  return (ln > ld ) ? ln : ld; 
}

template<>
CORE_INLINE int RealBigRat::sgn() const
{ return ker.sign(); }

// 31/05/2003 Joaquin Grech
template<>
CORE_INLINE std::string RealBigRat::toString(long prec, bool sci) const
{ return ker.toString(); }

//  cast operators
template<>
CORE_INLINE RealLong::operator double() const
{ return (double) ker; }

template<>
CORE_INLINE RealLong::operator float() const
{ return (float) ker; }

template<>
CORE_INLINE RealLong::operator long() const
{ return (long) ker; }

template<>
CORE_INLINE RealLong::operator int() const
{ return (int) ker; }

template<>
CORE_INLINE BigInt RealLong::BigIntValue() const
{ return BigInt(ker); }  

template<>
CORE_INLINE BigRat RealLong::BigRatValue() const
{ return BigRat(ker); }  

template<>
CORE_INLINE BigFloat RealLong::BigFloatValue() const
{ return BigFloat(ker); }  

//  cast operators
template<>
CORE_INLINE RealDouble::operator double() const
{ return (double) ker; }

template<>
CORE_INLINE RealDouble::operator float() const
{ return (float) ker; }

template<>
CORE_INLINE RealDouble::operator long() const
{ return (long) ker; }

template<>
CORE_INLINE RealDouble::operator int() const
{ return (int) ker; }

template<>
CORE_INLINE BigInt RealDouble::BigIntValue() const
{ return BigInt(ker); }  

template<>
CORE_INLINE BigRat RealDouble::BigRatValue() const
{ return BigRat(ker); }  

template<>
CORE_INLINE BigFloat RealDouble::BigFloatValue() const
{ return BigFloat(ker); }  

//  cast operators
template<>
CORE_INLINE RealBigInt::operator double() const
{ return ker.doubleValue(); }

template<>
CORE_INLINE RealBigInt::operator float() const
{ return (float)ker.doubleValue(); }

template<>
CORE_INLINE RealBigInt::operator long() const
{ return ker.longValue(); }

template<>
CORE_INLINE RealBigInt::operator int() const
{ return (int)ker.longValue(); }

template<>
CORE_INLINE BigInt RealBigInt::BigIntValue() const
{ return ker; }  

template<>
CORE_INLINE BigRat RealBigInt::BigRatValue() const
{ return BigRat(ker); }  

template<>
CORE_INLINE BigFloat RealBigInt::BigFloatValue() const
{ return BigFloat(ker); }  

//  cast operators
template<>
CORE_INLINE RealBigFloat::operator double() const
{ return ker.doubleValue(); }

template<>
CORE_INLINE RealBigFloat::operator float() const
{ return ker.floatValue(); }

template<>
CORE_INLINE RealBigFloat::operator long() const
{ return ker.longValue(); }

template<>
CORE_INLINE RealBigFloat::operator int() const
{ return ker.intValue(); }

template<>
CORE_INLINE BigInt RealBigFloat::BigIntValue() const
{ return ker.BigIntValue(); }  

template<>
CORE_INLINE BigRat RealBigFloat::BigRatValue() const
{ return ker.BigRatValue(); }  

template<>
CORE_INLINE BigFloat RealBigFloat::BigFloatValue() const
{ return ker; }  

//  cast operators
template<>
CORE_INLINE RealBigRat::operator double() const
{ return ker.doubleValue(); }

template<>
CORE_INLINE RealBigRat::operator float() const
{ return static_cast<float>(ker.doubleValue()); }

template<>
CORE_INLINE RealBigRat::operator long() const
{ return static_cast<long>(ker.doubleValue()); }

template<>
CORE_INLINE RealBigRat::operator int() const 
{ return static_cast<int>(ker.doubleValue()); }

template<>
CORE_INLINE BigInt RealBigRat::BigIntValue() const
{ return ker.BigIntValue(); }  

template<>
CORE_INLINE BigRat RealBigRat::BigRatValue() const
{ return ker; }  

template<>
CORE_INLINE BigFloat RealBigRat::BigFloatValue() const
{ return BigFloat(ker); }  


#endif

