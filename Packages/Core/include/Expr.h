/******************************************************************
 * Core Library Version 1.5, August 2002
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
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_EXPR_H
#define CORE_EXPR_H

#include "CoreImpl.h"
#include "CoreAux.h"
#include "Real.h"
#include "Filter.h"
#include "MemoryPool.h"

CORE_BEGIN_NAMESPACE

// These counters are incremented each time each bound is recognized as equal
// to the best one in computeBound().
extern unsigned int BFMSS_counter;
extern unsigned int Measure_counter;
extern unsigned int Cauchy_counter;
extern unsigned int LiYap_counter;
// These counters are incremented each time each bound is recognized as equal
// to the best one in computeBound(), and it's strictly the best.
extern unsigned int BFMSS_only_counter;
extern unsigned int Measure_only_counter;
extern unsigned int Cauchy_only_counter;
extern unsigned int LiYap_only_counter;
// This counter is incremented each time the precision needed matches the
// root bound.
extern unsigned int rootBoundHitCounter;

//  forward reference
class ExprRep;

/// \class Expr Expr.h
/// \brief Expr is a class of Expression in Level 3
class Expr {
public:   
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  Expr();
  /// constructor for \c int
  Expr(int);
  /// constructor for <tt>unsigned int</tt>
  Expr(unsigned int);
  /// constructor for \c long
  Expr(long);
  /// constructor for <tt>unsigned long</tt>
  Expr(unsigned long);
  /// constructor for \c float
  Expr(float);
  /// constructor for \c double
  Expr(double);
  /// constructor for \c BigInt
  Expr(const BigInt &);
  /// constructor for \c BigRat
  Expr(const BigRat &);
  /// constructor for \c BigFloat
  Expr(const BigFloat &);
  /// constructor for \c string
  /** construct Expr from a string representation \a s 
   * with precision \a prec */
  Expr(const char *s, const extLong& prec = defInputDigits); 
  /// constructor for \c Real
  Expr(const Real &);
  /// copy constructor
  Expr(const Expr&); 
  /// destructor
  ~Expr();          
  //@}

  /// \name Aprroximation Function
  //@{
  /// Compute approximation to combined precision [\a r, \a a].
  /** Here is the definition of what this means: 
       If e is the exact value and ee is the approximate value,
       then  |e - ee| <= 2^{-a} or  |e - ee| <= 2^{-r} |e|. */
  const Real & approx(const extLong& relPrec = defRelPrec,
	              const extLong& absPrec = defAbsPrec) const;
  //@}

  /// \name Helper Functions
  //@{
  /// get current approximate value
  BigFloat getBigFloat() const;  	
  /// get exponent of current approximate value
  long getExponent() const; 
  /// get mantissa of current approximate value
  BigInt getMantissa() const;
  /// get the sign
  int sign() const;
  //@}

  /// \name Conversion Functions
  //@{
  /// convert to \c int
  int toInt() const ;
  /// convert to \c long
  long toLong() const ;
  /// convert to \c float
  float toFloat() const;
  /// convert to \c double
  double toDouble() const;
  /// convert to \c string
  /** give decimal string representation */
  const char *toString() const;
  //@}

  /// \name Assignment Operators
  //@{ 
  /// = operator
  Expr& operator=(const Expr&);
  /// += operator
  Expr& operator+=(const Expr&);
  /// -= operator
  Expr& operator-=(const Expr&);
  /// *= operator
  Expr& operator*=(const Expr&);
  /// /= operator
  Expr& operator/=(const Expr&);
  //@}

  /// \name Increment, Decrement, and Unary Minus Operators
  //@{
  /// left increment operator (++i)
  Expr& operator++();
  /// right increment operator (i++)
  Expr operator++(int);
  /// left decrement operator (--i)
  Expr& operator--();
  /// right deccrement operator (i--)
  Expr operator--(int);
  /// unary minus
  Expr operator-() const;
  //@}
  
  /// \name Arithematic Operators
  //@{
  /// addition
  friend Expr operator+(const Expr&, const Expr&);
  /// substraction
  friend Expr operator-(const Expr&, const Expr&);
  /// multiplication
  friend Expr operator*(const Expr&, const Expr&);
  /// division
  friend Expr operator/(const Expr&, const Expr&);
  /// square root
  friend Expr sqrt(const Expr&);
  //@}

  /// \name Comparison Operators
  //@{
  /// operator ==
  friend bool operator==(const Expr&, const Expr&);
  /// operator !=
  friend bool operator!=(const Expr&, const Expr&);
  /// operator <
  friend bool operator< (const Expr&, const Expr&);
  /// operator <=
  friend bool operator<=(const Expr&, const Expr&);
  /// operator <
  friend bool operator> (const Expr&, const Expr&);
  /// operator >=
  friend bool operator>=(const Expr&, const Expr&);
  //@}
  
  /// \name Builtin Functions
  //@{
  /// sign function
  friend int sign(const Expr&);
  /// isZero function
  friend bool isZero(const Expr&);
  /// compare function
  /** compare two Expr \a e1 and \a e2, return
   * \retval -1 if e1 < e2,
   * \retval 0 if e1 = e2,
   * \retval 1 if e1 > e2. */
  friend int compare(const Expr& e1, const Expr& e2);
  /// floor function
  friend BigInt floor(const Expr&);
  /// ceil function
  friend BigInt ceil(const Expr&);
  /// power function
  friend Expr pow(const Expr&, unsigned long);
  /// power function (same as pow())
  friend Expr power(const Expr&, unsigned long n);
  /// absolute value function
  friend Expr abs(const Expr&);
  /// absolute value function (same as abs())
  friend Expr fabs(const Expr&);
  //@}

  /// \name I/O Stream
  //@{
  /// write to ostream 
  friend std::ostream& operator<<(std::ostream&, const Expr&);
  /// read from istream
  friend std::istream& operator>>(std::istream&, Expr&);
  //@}

public:
  /// \name Debug Helper Function
  //@{
  /// debug function
  void  debug(int mode = TREE_MODE, int level = DETAIL_LEVEL, 
              int depthLimit = INT_MAX) const;
  //@}
  
  /// debug information levels
  enum {LIST_MODE, TREE_MODE, SIMPLE_LEVEL, DETAIL_LEVEL};

  /// \name Deprecated Functions
  //@{
  /// sign function 
  /** get sign (forces computation of exact flags if necessary) 
      \deprecated using sign() instead of */
  int getSign() const;
  /// convert to \c int
  /** \deprecated using toInt() instead of */
  int intValue() const ;
  /// convert to \c long
  /** \deprecated using toLong() instead of */
  long longValue() const ;
  /// convert to \c float
  /** \deprecated using toFloat() instead of */
  float floatValue() const;
  /// convert to \c double
  /** \deprecated using toDouble() instead of */
  double doubleValue() const;
  //@}
  
  /// return Expr(0)
  static const Expr& getZero();
  ExprRep* getRep() const { return rep; }
private:
  Expr(ExprRep* p) : rep(p) {}
protected:
  ExprRep* rep; ///< handle to the "real" representation

};// class Expr

/// \struct NodeInfo
/// \brief store information of a node
struct NodeInfo {
  Real     appValue;		///< current approximate value
  bool     appComputed;  	///< true if the approx value been computed
  bool     flagsComputed;  	///< true if rootBound parameters have been computed
  extLong  knownPrecision; 	///< Precision achieved by current approx value

#ifdef DEBUG
  extLong relPrecision;
  extLong absPrecision;
  unsigned long numNodes; 
#endif

  /// d_e bounds the degree of the minimal polynomial of a DAG expression
  /** Basically, d_e is equal to 2^k where k is the number of square-root nodes
   *   in the DAG.  If there are other kinds of non-linear nodes, this is
   *   generalized accordingly.   */
  extLong d_e; 

  bool visited;   	///< flag in counting # of sqrts
  int sign; 		///< sign of the value being represented.

  extLong  uMSB; ///< upper bound of the position of Most Significant Bit
  extLong  lMSB; ///< lower bound of the position of Most Significant Bit

  // For the degree-length method mentioned in Chee's book.
  /* the degree of defining polynomial P(X) obtained from Resultant calculus 
   * (deprecated now) */
  // extLong degree;

  extLong length; ///< || P(X) ||
  extLong measure; ///< height

  // For our new bound.
  /// 2^{high(E)} is an UPPER bound for the moduli of ALL conjugates of E.
  /** In our papers, high is equal to log_2(\overline{\mu(E)}). */
  extLong high;
  /// 2^{-low(E)} is an LOWER bound for the moduli of ALL NON_ZERO conjugate of E.
  /** BE CAREFUL!  NOTE THAT UNLIKE "high", the sign of low is negated here! 
      In our papers, low is equal to -log_2(\underline{\nu(E)})  */
  extLong low;   
  
  /// \brief upper bound of the leading coefficient of minimal defining
  ///        polynomial of $E$.
  extLong lc;
  /// \brief upper bound of the last non-zero coefficient of minimal defining 
  ///        polynomial of $E$.
  extLong tc;  

  // For the 2-ary BFMSS bound.
  extLong v2p, v2m;
  // For the 5-ary BFMSS bound.
  extLong v5p, v5m;
  /// 2^u25 is an upper bound for the moduli of all the conjugates of U(E)
  /** where E = 2^v2*5^v5*U(E)/L(E), U(E) and L(E) are division-free. */
  extLong u25;
  /// 2^l25 is an upper bound for the moduli of all the conjugates of L(E)
  /** where E = 2^v2*5^v5*U(E)/L(E), U(E) and L(E) are division-free. */
  extLong l25;

  int ratFlag;		///< rational flag
  BigRat* ratValue;	///< rational value

  /// default constructor
  NodeInfo();
};

/// \class ExprRep
/// \brief The sharable, internal representation of expressions
class ExprRep {
public:
  /// \name Constructor and Destructor
  //@{
  /// default constructor
  ExprRep();
  /// virtual destructor for this base class
  virtual ~ExprRep() {
    if (nodeInfo != NULL) // This check is only for optimization.
      delete nodeInfo;
  }
  //@}
  
  /// \name Reference Counting
  //@{
  /// increase reference counter
  void incRefCount() { ++refCount; }
  /// decrease reference counter
  void decRefCount() { if ((--refCount) == 0) delete this; }
  /// check whether reference counter == 1
  int isUnique() const { return refCount == 1; }
  //@}
 
  /// \name Helper Functions
  //@{
  /// Get the approximate value
  const Real & getAppValue(const extLong& relPrec = defRelPrec,
		  const extLong& absPrec = defAbsPrec);
  /// Get the sign.
  int getSign(); 
  int getExactSign();

  const Real& appValue() const { return nodeInfo->appValue; }
  Real& appValue() { return nodeInfo->appValue; }
  
  const bool& appComputed() const { return nodeInfo->appComputed; }
  bool& appComputed() { return nodeInfo->appComputed; }

  const bool& flagsComputed() const { return nodeInfo->flagsComputed; }
  bool& flagsComputed() { return nodeInfo->flagsComputed; }

  const extLong& knownPrecision() const { return nodeInfo->knownPrecision; }
  extLong& knownPrecision() { return nodeInfo->knownPrecision; }

#ifdef DEBUG  
  const extLong& relPrecision() const { return nodeInfo->relPrecision; }
  extLong& relPrecision() { return nodeInfo->relPrecision; }

  const extLong& absPrecision() const { return nodeInfo->absPrecision; }
  extLong& absPrecision() { return nodeInfo->absPrecision; }

  const unsigned long& numNodes() const { return nodeInfo->numNodes; }
  unsigned long& numNodes() { return nodeInfo->numNodes; }
#endif

  const extLong& d_e() const { return nodeInfo->d_e; }
  extLong& d_e() { return nodeInfo->d_e; }

  const bool& visited() const { return nodeInfo->visited; }
  bool& visited() { return nodeInfo->visited; }

  const int& sign() const { return nodeInfo->sign; }
  int& sign() { return nodeInfo->sign; }

  const extLong& uMSB() const { return nodeInfo->uMSB; }
  extLong& uMSB() { return nodeInfo->uMSB; }

  const extLong& lMSB() const { return nodeInfo->lMSB; }
  extLong& lMSB() { return nodeInfo->lMSB; }

  const extLong& length() const { return nodeInfo->length; }
  extLong& length() { return nodeInfo->length; }

  const extLong& measure() const { return nodeInfo->measure; }
  extLong& measure() { return nodeInfo->measure; }

  const extLong& high() const { return nodeInfo->high; }
  extLong& high() { return nodeInfo->high; }

  const extLong& low() const { return nodeInfo->low; }
  extLong& low() { return nodeInfo->low; }

  const extLong& lc() const { return nodeInfo->lc; }
  extLong& lc() { return nodeInfo->lc; }

  const extLong& tc() const { return nodeInfo->tc; }
  extLong& tc() { return nodeInfo->tc; }

  const extLong& v2p() const { return nodeInfo->v2p; }
  extLong& v2p() { return nodeInfo->v2p; }

  const extLong& v2m() const { return nodeInfo->v2m; }
  extLong& v2m() { return nodeInfo->v2m; }

  extLong v2() const { return v2p()-v2m(); }

  const extLong& v5p() const { return nodeInfo->v5p; }
  extLong& v5p() { return nodeInfo->v5p; }

  const extLong& v5m() const { return nodeInfo->v5m; }
  extLong& v5m() { return nodeInfo->v5m; }

  extLong v5() const { return v5p()-v5m(); }

  const extLong& u25() const { return nodeInfo->u25; }
  extLong& u25() { return nodeInfo->u25; }
  
  const extLong& l25() const { return nodeInfo->l25; }
  extLong& l25() { return nodeInfo->l25; }
  
  const int& ratFlag() const { return nodeInfo->ratFlag; }
  int& ratFlag() { return nodeInfo->ratFlag; }

  const BigRat* ratValue() const { return nodeInfo->ratValue; }
  BigRat*& ratValue() { return nodeInfo->ratValue; }

  /// Get BigFloat
  BigFloat getBigFloat();
  /// represent as a string (not implemented yet)
  const char *toString() const { return "to be implemented"; }
  //@}

  /// \name Debug functions
  //@{
  /// dump the contents in this DAG node
  const std::string dump(int = OPERATOR_VALUE) const; 
  /// print debug information in list mode
  virtual void debugList(int level, int depthLimit) const = 0;
  /// print debug information in tree mode
  virtual void debugTree(int level, int indent, int depthLimit) const = 0; 
  //@}

  /// \name I/O Stream
  //@{
  friend std::ostream& operator<<(std::ostream&, ExprRep&);
  //@}

  CORE_MEMORY(ExprRep)

private:
  unsigned refCount;    // reference count

public:
  enum {OPERATOR_ONLY, VALUE_ONLY, OPERATOR_VALUE, FULL_DUMP};

  NodeInfo* nodeInfo;	///< node information
  filteredFp ffVal;	///< filtered value

  /// \name Approximation Functions
  //@{
  /// initialize nodeInfo
  virtual void initNodeInfo() = 0;
  /// compute the sign, uMSB, lMSB, etc.
  virtual void computeExactFlags() = 0;
  /// compute the minimal root bound
  extLong computeBound();
  /// driver function to approximate
  void approx(const extLong& relPrec, const extLong& absPrec); 
  /// compute an approximate value satifying the specified precisions
  virtual void computeApproxValue(const extLong&, const extLong&) = 0;
  /// Test whether the current approx. value satisfies [relPrec, absPrec]
  bool withinKnownPrecision(const extLong&, const extLong&);
  //@}

  /// \name Misc Functions
  //@{
  /// reduce current node
  void reduceToBigRat(const BigRat&);
  /// reduce current node
  void reduceTo(const ExprRep*);
  /// reduce current node to zero
  void reduceToZero();
  /// return operator string
  virtual const std::string op() const { return "UNKNOWN"; }
  //@}

  /// \name Degree Bound Functions
  //@{
  /// compute "d_e" based on # of sqrts
  extLong degreeBound();
  /// count # of sqrts in the DAG
  virtual unsigned long count() = 0; 
  /// reset the flag "visited"
  virtual void clearFlag() = 0;
  //@}
#ifdef DEBUG
  virtual unsigned long dagSize() = 0;
  virtual void fullClearFlag() = 0;
#endif 
};//ExprRep

/// \class ConstRep
/// \brief constant node
class ConstRep : public ExprRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  ConstRep() {}
  /// destructor
  ~ConstRep() {}
  //@}

  /// \name Debug Functions
  //@{
  /// print debug information in list mode
  void debugList(int level, int depthLimit) const;
  /// print debug information in tree mode
  void debugTree(int level, int indent, int depthLimit) const; 
  //@}
  CORE_MEMORY(ConstRep)
protected:
  /// initialize nodeInfo
  void initNodeInfo();
  /// return operator in string
  const std::string op() const { return "C"; }
  /// count # of square roots
  unsigned long count() { return 0; }
  /// clear visited flag
  void clearFlag() { visited() = false; }
#ifdef DEBUG
  unsigned long dagSize();
  void fullClearFlag();
#endif 
};

/// \class ConstDoubleRep
/// \brief constant node
class ConstDoubleRep : public ConstRep{
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  ConstDoubleRep() { }
  /// constructor for all \c double type
  ConstDoubleRep(double d) { ffVal = d; }
  /// destructor
  ~ConstDoubleRep() {}
  //@}
  CORE_MEMORY(ConstDoubleRep)
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
};

/// \class ConstRealRep
/// \brief constant node
class ConstRealRep : public ConstRep{
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  ConstRealRep() : value(CORE_REAL_ZERO) { }
  /// constructor for all \c Real type
  ConstRealRep(const Real &);
  /// destructor
  ~ConstRealRep() {}
  //@}
  CORE_MEMORY(ConstRealRep)
private:
  Real value; ///< internal representation of node 
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
};

/// \class UnaryOpRep
/// \brief unary operator node
class UnaryOpRep : public ExprRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  UnaryOpRep(ExprRep* c) : child(c)  { child->incRefCount(); }
  /// destructor
  virtual ~UnaryOpRep() { child->decRefCount(); }
  //@}

  /// \name Debug Functions
  //@{
  /// print debug information in list mode
  void debugList(int level, int depthLimit) const;
  /// print debug information in tree mode
  void debugTree(int level, int indent, int depthLimit) const; 
  //@}

  CORE_MEMORY(UnaryOpRep)
protected:
  ExprRep* child; ///< pointer to its child node
  /// initialize nodeInfo
  virtual void initNodeInfo();
  /// clear visited flag
  void clearFlag();
#ifdef DEBUG
  unsigned long dagSize();
  void fullClearFlag();
#endif 
};

/// \class NegRep
/// \brief unary minus operator node
class NegRep : public UnaryOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  NegRep(ExprRep* c) : UnaryOpRep(c)
  {
    ffVal = - child->ffVal;
  }
  /// destructor
  ~NegRep() {}
  //@}

  CORE_MEMORY(NegRep)
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const { return "Neg"; }
  /// count # of square roots
  unsigned long count();
};

/// \class SqrtRep
/// \brief squartroot operator node
class SqrtRep : public UnaryOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  SqrtRep(ExprRep* c) : UnaryOpRep(c)
  {
    ffVal = (child->ffVal).sqrt();
  }
  /// destructor
  ~SqrtRep() {}
  //@}

  CORE_MEMORY(SqrtRep)
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const { return "Sqrt"; }
  /// count # of square roots
  unsigned long count();
};

/// \class BinOpRep
/// \brief binary operator node
class BinOpRep : public ExprRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  BinOpRep(ExprRep* f, ExprRep* s) : first(f), second(s)
  { first->incRefCount(); second->incRefCount(); }
  /// destructor
  virtual ~BinOpRep()
  { first->decRefCount(); second->decRefCount(); }
  //@}

  /// \name Debug Functions
  //@{
  /// print debug information in list mode
  void debugList(int level, int depthLimit) const;
  /// print debug information in tree mode
  void debugTree(int level, int indent, int depthLimit) const; 
  //@}

  CORE_MEMORY(BinOpRep)
protected:
  ExprRep* first;  ///< first operand
  ExprRep* second; ///< second operand
  
  /// initialize nodeInfo
  virtual void initNodeInfo();
  /// clear visited flags 
  void clearFlag();
  /// count # of square roots
  unsigned long count();
#ifdef DEBUG
  unsigned long dagSize();
  void fullClearFlag();
#endif 
};

/// \struct Add
/// \brief "functor" class used as parameter to AddSubRep<>
struct Add {
  /// name
  static const char* name;

  /// unary operator
  template <class T>
  const T& operator()(const T& t) const
  { return t; }

  /// binary operator
  template <class T>
  T operator()(const T& a, const T& b) const
  { return a+b; }
};

/// \struct Sub
/// \brief "functor" class used as parameter to AddSubRep<>
struct Sub {
  /// name
  static const char* name;

  /// unary operator
  template <class T>
  T operator()(const T& t) const
  { return -t; }

  /// binary operator
  template <class T>
  T operator()(const T& a, const T& b) const
  { return a-b; }
};

/// \class AddSubRep
/// \brief template class where operator is supposed to be Add or Sub
template <class Operator>
class AddSubRep : public BinOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  AddSubRep(ExprRep* f, ExprRep* s) : BinOpRep(f, s)
  {
    ffVal = Op(first->ffVal, second->ffVal);
  }
  /// destructor
  ~AddSubRep() {}
  //@}

  CORE_MEMORY(AddSubRep)
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const { return Operator::name; }
private:
  static Operator Op;
};

template <class Operator>
Operator AddSubRep<Operator>::Op;

/// \typedef AddRep
/// \brief AddRep for easy of use
typedef AddSubRep<Add> AddRep;

/// \typedef SubRep
/// \brief SuRep for easy of use
typedef AddSubRep<Sub> SubRep;

/// \class MultRep
/// \brief multiplication operator node
class MultRep : public BinOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  MultRep(ExprRep* f, ExprRep* s) : BinOpRep(f, s)
  {
    ffVal = first->ffVal * second->ffVal;
  }
  /// destructor
  ~MultRep() {}
  //@}

  CORE_MEMORY(MultRep)
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const { return "*"; }
};

/// \class DivRep
/// \brief division operator node
class DivRep : public BinOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  DivRep(ExprRep* f, ExprRep* s) : BinOpRep(f, s)
  {
    ffVal = first->ffVal / second->ffVal;
  }
  /// destructor
  ~DivRep() {}
  //@}

  CORE_MEMORY(DivRep)
protected:
  /// compute sign and MSB
  void computeExactFlags();
  /// compute approximation value
  void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const { return "/"; }
};

BigInt floor(const Expr&, Expr&);
Expr pow(const Expr&, unsigned long);

#ifdef CORE_ENABLE_INLINES
  #include "Expr.inl"
#else
  // friend functions for Expr class
  // (need declarations in case they are not inlined)
  Expr operator+(const Expr&, const Expr&);
  Expr operator-(const Expr&, const Expr&);
  Expr operator*(const Expr&, const Expr&);
  Expr operator/(const Expr&, const Expr&);
  Expr sqrt(const Expr&);
  bool operator==(const Expr&, const Expr&);
  bool operator!=(const Expr&, const Expr&);
  bool operator< (const Expr&, const Expr&);
  bool operator<=(const Expr&, const Expr&);
  bool operator> (const Expr&, const Expr&);
  bool operator>=(const Expr&, const Expr&);
  int sign(const Expr&);
  bool isZero(const Expr&);
  int compare(const Expr& e1, const Expr& e2);
  BigInt floor(const Expr&);
  BigInt ceil(const Expr&);
  Expr power(const Expr&, unsigned long n);
  Expr abs(const Expr&);
  Expr fabs(const Expr&);
  std::ostream& operator<<(std::ostream&, const Expr&);
  std::istream& operator>>(std::istream&, Expr&);
#endif

#define CORE_EXPR_ZERO Expr::getZero()
  
CORE_END_NAMESPACE
#endif
