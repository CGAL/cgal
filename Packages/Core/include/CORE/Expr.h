/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: Expr.h
 *
 * Written by 
 *       Koji Ouchi <ouchi@simulation.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *       Igor Pechtchanski <pechtcha@cs.nyu.edu>
 *       Vijay Karamcheti <vijayk@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Sylvain Pion <pion@cs.nyu.edu> 
 *       Vikram Sharma<sharma@cs.nyu.edu>
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_EXPR_H
#define CORE_EXPR_H

#include <CORE/ExprRep.h>

CORE_BEGIN_NAMESPACE

/// \class Expr Expr.h
/// \brief Expr is a class of Expression in Level 3
class Expr {
public:   
  ExprRep* rep; ///< handle to the "real" representation
  
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  Expr() { rep = new ConstDoubleRep(); }
  /// copy constructor
  Expr(const Expr& e) { rep = e.rep;  rep->incRefCount(); } 

  /// constructor for <tt>int</tt>
  Expr(int i) { rep = new ConstDoubleRep(i); }
  /// constructor for <tt>unsigned int</tt>
  Expr(unsigned int ui) { rep = new ConstDoubleRep(ui); }
  
  /// constructor for <tt>long</tt>
  Expr(long l) { rep = new ConstRealRep(Real(l)); }
  /// constructor for <tt>unsigned long</tt>
  Expr(unsigned long ul) { rep = new ConstRealRep(Real(ul)); }
  
  /// constructor for <tt>float</tt>
  /** \note the results of this constructor can be somewhat unpredictable.
   *  One might assume that new Expr(.1) is exactly equal to .1, but it is
   *  actually equal to
   *      .1000000000000000055511151231257827021181583404541015625.
   *  This is so because .1 cannot be represented exactly as a double
   *  (or, for that matter, as a binary fraction of any finite length). 
   *  Thus, the long value that is being passed in to the constructor is not
   *  exactly equal to .1, appearances nonwithstanding. 
   */
  Expr(float f) { // check for valid numbers (i.e., not infinite and not NaN)
    if (!finite(f)) {
      std::cerr << " ERROR : constructed an invalid float! " << std::endl;
      if (AbortFlag) abort();
      InvalidFlag = -1;
    }
  	rep = new ConstDoubleRep(f);
  }
  /// constructor for <tt>double</tt>
  Expr(double d) { // check for valid numbers (i.e., not infinite and not NaN)
    if (!finite(d)) {
      std::cerr << " ERROR : constructed an invalid double! " << std::endl;
      if (AbortFlag) abort();
      InvalidFlag = -2;
    }
  	rep = new ConstDoubleRep(d);
  }
  
  /// constructor for <tt>BigInt</tt>
  Expr(const BigInt& I) { rep = new ConstRealRep(Real(I)); }
  /// constructor for <tt>BigRat</tt>
  Expr(const BigRat& R) { rep = new ConstRealRep(Real(R)); }
  /// constructor for <tt>BigFloat</tt>
  Expr(const BigFloat& F) { rep = new ConstRealRep(Real(F)); }
  
  /// constructor for <tt>const char*</tt>
  /** construct Expr from a string representation \a s 
   * with precision \a prec. It is perfectly predictable:
   * new Expr(".1") is exactly equal to .1, as one would expect. Therefore,
   * it is generally recommended that the (String) constructor be used in
   * preference to the (double) constructor.
   */
  Expr(const char *s, const extLong& prec = defInputDigits)
  { rep = new ConstRealRep(Real(s, prec)); }

  /// constructor for <tt>std::string</tt>
  Expr(const std::string& s, const extLong& prec = defInputDigits)
  { rep = new ConstRealRep(Real(s, prec)); }
  
  /// constructor for <tt>Real</tt>
  Expr(const Real &r) { rep = new ConstRealRep(r); }
  
  /// constructor for Polynomial node (n-th root)
  /** default value n=0 means the first positive root */
  template <class NT>
  Expr(const Polynomial<NT>& p, int n = 0)
  { rep = new ConstPolyRep<NT>(p, n); }

  /// constructor for Polynomial node (root in Interval <tt>I</tt>)
  template <class NT>
  Expr(const Polynomial<NT>& p, const BFInterval& I)
  { rep = new ConstPolyRep<NT>(p, I); }
  
  /// constructor for ExprRep
  Expr(ExprRep* p) : rep(p) {}

  /// destructor
  ~Expr() { rep->decRefCount(); }          
  //@}

  /// \name Assignment Operators
  //@{ 
  /// = operator
  Expr& operator=(const Expr& e) {
    if (this == &e) return *this;
    e.rep->incRefCount(); rep->decRefCount();
    rep = e.rep; return *this;
  }
  //@}
  
  /// \name Compound Assignment Operators
  //@{ 
  /// += operator
  Expr& operator+=(const Expr& e) {
    ExprRep *old = rep;
    rep = new AddRep(rep, e.rep);
    old->decRefCount();
    return *this;
  }
  /// -= operator
  Expr& operator-=(const Expr& e) {
    ExprRep *old = rep;
    rep = new SubRep(rep, e.rep);
    old->decRefCount();
    return *this;
  }
  /// *= operator
  Expr& operator*=(const Expr& e) {
    ExprRep *old = rep;
    rep = new MultRep(rep, e.rep);
    old->decRefCount();
    return *this;
  }
  /// /= operator
  Expr& operator/=(const Expr& e) {
  	if ((e.rep)->getSign() == 0) {
      std::cerr << " ERROR : division by zero ! " << std::endl;
      if (AbortFlag) abort();
      InvalidFlag = -3;
    }
    ExprRep *old = rep;
    rep = new DivRep(rep, e.rep);
    old->decRefCount();
    return *this;
  }
  //@}

  /// \name Unary Minus, Increment and Decrement Operators
  //@{
  /// unary minus
  Expr operator-() const { return Expr(new NegRep(rep)); }
  /// left increment operator (++i)
  Expr& operator++() { *this += 1; return *this; }
  /// right increment operator (i++)
  Expr operator++(int) { Expr t = *this; *this += 1; return t; }
  /// left decrement operator (--i)
  Expr& operator--() { *this -= 1; return *this; }
  /// right deccrement operator (i--)
  Expr operator--(int) { Expr t = *this; *this -= 1; return t; }
  //@}

  /// \name String Conversion Functions
  //@{
  /// set value from <tt>const char*</tt>
  void fromString(const char* s, const extLong& prec = defInputDigits)
  { *this = Expr(s, prec); }
  /// convert to <tt>std::string</tt>
  /** give decimal string representation */
  std::string toString(long prec=defOutputDigits, bool sci=false) const
  { return rep->toString(prec, sci); }		  
  //@}
  //
  
  /// \name Conversion Functions
  //@{
  /// convert to \c int
  int intValue() const { return (this->approx(64, 1024)).intValue(); }
  /// convert to \c long
  long longValue() const { return (this->approx(64, 1024)).longValue(); }
  /// convert to \c float
  float floatValue() const { return (this->approx(53, 1024)).floatValue(); }
  /// convert to \c double
  /** chen: - use equivalent precision (rel:53, abs: 1024)
    as in IEEE double. enforce an evaluation in case
    before this has been done before casting. */
  double doubleValue() const { return (this->approx(53, 1024)).doubleValue(); }
  /// convert to an interval defined by a pair of \c double
  /** If value is exact, the two \c double will coincide
   */
  void doubleInterval(double & lb, double & ub) const;
  /// convert to \c BigInt (approximate it first!)
  BigFloat BigIntValue() const { return rep->BigIntValue(); }  	
  /// convert to \c BigRat (approximate it first!)
  BigFloat BigRatValue() const { return rep->BigRatValue(); }  	
  /// convert to \c BigFloat (approximate it first!)
  /** Ought to allow BigFloatValue() take an optional precision argument */
  BigFloat BigFloatValue() const { return rep->BigFloatValue(); }  	
  //@}

  /// \name Approximation Function
  //@{
  /// Compute approximation to combined precision [\a r, \a a].
  /** Here is the definition of what this means: 
       If e is the exact value and ee is the approximate value,
       then  |e - ee| <= 2^{-a} or  |e - ee| <= 2^{-r} |e|. */
  const Real & approx(const extLong& relPrec = defRelPrec,
	              const extLong& absPrec = defAbsPrec) const
  { return rep->getAppValue(relPrec, absPrec); }
  //@}

  /// \name Helper Functions
  //@{
  /// get the sign
  int sign() const { return rep->getSign(); }
  /// is zero?
  bool isZero() const { return sign() == 0; }
  /// absolute value
  Expr abs() const { Expr x = (sign() >= 0) ? (*this) : -(*this); return x; }
  
  /// compare function
  int cmp(const Expr& e) const 
  { return rep == e.rep ? 0 : SubRep(rep, e.rep).getSign(); }
  
  /// get exponent of current approximate value
  long getExponent() const { return BigFloatValue().exp(); } 
  /// get mantissa of current approximate value
  BigInt getMantissa() const { return BigFloatValue().m(); }
  /// return rep pointer:
  ExprRep* getRep() const { return rep; }
  //@}

  /// return Expr(0)
  static const Expr& getZero();
public:
  /// \name Debug Helper Function
  //@{
  /// debug function
  void  debug(int mode = TREE_MODE, int level = DETAIL_LEVEL, 
              int depthLimit = INT_MAX) const;
  //@}
  /// debug information levels
  enum {LIST_MODE, TREE_MODE, SIMPLE_LEVEL, DETAIL_LEVEL};
};// class Expr

#define CORE_EXPR_ZERO Expr::getZero()

/// I/O Stream operator<<
inline std::ostream& operator<<(std::ostream& o, const Expr& e)
{ o << *(e.rep); return o; }
/// I/O Stream operator>>
inline std::istream& operator>>(std::istream& i, Expr& e) {
  Real rVal; i >> rVal; // precision is = defInputDigits
  if (i)  e = rVal;		// only assign when reading is successful.
  return i;
}

/// floor function
BigInt floor(const Expr&, Expr&);
/// power function
Expr pow(const Expr&, unsigned long);

/// addition
inline Expr operator+(const Expr& e1, const Expr& e2)
{ return Expr(new AddRep(e1.rep, e2.rep)); }
/// substraction
inline Expr operator-(const Expr& e1, const Expr& e2)
{ return Expr(new SubRep(e1.rep, e2.rep)); }
/// multiplication
inline Expr operator*(const Expr& e1, const Expr& e2)
{ return Expr(new MultRep(e1.rep, e2.rep)); }
/// division
inline Expr operator/(const Expr& e1, const Expr& e2) {
  if (e2.sign() == 0) {
    std::cerr << " ERROR : division by zero ! " << std::endl;
    if (AbortFlag) abort();
    InvalidFlag = -4;
  }
  return Expr(new DivRep(e1.rep, e2.rep));
}
/// modulo operator
inline Expr operator%(const Expr& e1, const Expr& e2)
{ Expr result; floor(e1/e2, result); return result; }

/// operator ==
/** this is inefficient if you compare to zero: 
 *  e.g., if (e != 0) {...} use e.isZero() instead */
inline bool operator==(const Expr& e1, const Expr& e2)
{ return e1.cmp(e2) == 0; }
/// operator !=
inline bool operator!=(const Expr& e1, const Expr& e2)
{ return e1.cmp(e2) != 0; }
/// operator <
inline bool operator< (const Expr& e1, const Expr& e2)
{ return e1.cmp(e2) < 0; }
/// operator <=
inline bool operator<=(const Expr& e1, const Expr& e2)
{ return e1.cmp(e2) <= 0; }
/// operator <
inline bool operator> (const Expr& e1, const Expr& e2)
{ return e1.cmp(e2) > 0; }
/// operator >=
inline bool operator>=(const Expr& e1, const Expr& e2)
{ return e1.cmp(e2) >= 0; }

/// return sign
inline int sign(const Expr& e) { return e.sign(); }
/// is zero?
inline bool isZero(const Expr& e) { return e.isZero(); }
/// compare
/** compare two Expr \a e1 and \a e2, return
 * \retval -1 if e1 < e2,
 * \retval 0 if e1 = e2,
 * \retval 1 if e1 > e2. */
inline int cmp(const Expr& e1, const Expr& e2) { return e1.cmp(e2); }
/// absolute value
inline Expr abs(const Expr& x) { return x.abs(); }
/// absolute value (same as abs)
inline Expr fabs(const Expr& x) { return abs(x); }
/// floor 
inline BigInt floor(const Expr& e) { Expr tmp; return floor(e, tmp); }
/// ceiling
inline BigInt ceil(const Expr& e) { return -floor(-e); }
/// power
inline Expr power(const Expr& e, unsigned long p) { return pow(e, p); }
/// divisibility predicate
inline bool isDivisible(const Expr& e1, const Expr& e2)
{ Expr result; floor(e1/e2, result); return (result.sign() == 0); }
/// square root
inline Expr sqrt(const Expr& e) {
  if (e.sign() < 0) {
    std::cerr << " ERROR : sqrt of negative value ! " << std::endl;
    if (AbortFlag) abort();
    InvalidFlag = -5;
  }
  return Expr(new SqrtRep(e.rep));
}

/// helper function for constructing Polynomial node (n-th node)
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, int n = 0)
{ return Expr(p, n); }

/// helper function for constructing Polynomial node
template <class NT>
inline Expr rootOf(const Polynomial<NT>& p, const BFInterval& I)
{ return Expr(p, I); }

/// constructor for Polynomial node of the form x^m - n (i.e., radicals)
template <class NT>
inline Expr radical(const NT& n, int m) {
  Polynomial<NT> Q(m);
  BFInterval I(0, n);
  Q.setCoeff(0, -n); Q.setCoeff(m, 1);
  return Expr(new ConstPolyRep<NT>(Q, I));
}

CORE_END_NAMESPACE
#endif
