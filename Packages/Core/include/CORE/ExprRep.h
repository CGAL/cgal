/******************************************************************
 * Core Library Version 1.6, June 2003
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: ExprRep.h
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

#ifndef CORE_EXPRREP_H
#define CORE_EXPRREP_H

#include "CoreImpl.h"
#include "CoreAux.h"
#include "Real.h"
#include "Filter.h"
#include "MemoryPool.h"
#include "poly/Sturm.h"

CORE_BEGIN_NAMESPACE

#ifdef DEBUG_BOUND
// These counters are incremented each time each bound is recognized as equal
// to the best one in computeBound().
extern unsigned int BFMSS_counter;
extern unsigned int Measure_counter;
// extern unsigned int Cauchy_counter;
extern unsigned int LiYap_counter;
// These counters are incremented each time each bound is recognized as equal
// to the best one in computeBound(), and it's strictly the best.
extern unsigned int BFMSS_only_counter;
extern unsigned int Measure_only_counter;
// extern unsigned int Cauchy_only_counter;
extern unsigned int LiYap_only_counter;
// This counter is incremented each time the precision needed matches the
// root bound.
extern unsigned int rootBoundHitCounter;
#endif

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

  // extLong length; ///< length is really lg(|| P(X) ||)
  extLong measure; ///< measure is really lg(Measure)

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

//  forward reference
// class Expr;

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

//  const extLong& length() const { return nodeInfo->length; }
//  extLong& length() { return nodeInfo->length; }

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
  BigInt BigIntValue();
  BigRat BigRatValue();
  BigFloat BigFloatValue();
  /// represent as a string in decimal value
  // toString() Joaquin Grech 31/5/2003
  std::string toString(long prec=defOutputDigits, bool sci=false) const
  {	return nodeInfo->appValue.toString(prec,sci); }
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
  /// count actually computes the degree bound of current node.
  virtual extLong count() = 0; 
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
  virtual void initNodeInfo();
  /// return operator in string
  const std::string op() const { return "C"; }
  /// count returns the degree of current node
  extLong count() { return d_e(); }
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

/// \class Constant Polynomial Node
/// \brief template class where NT is supposed to be some number type
template <class NT>
class ConstPolyRep : public ConstRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  ConstPolyRep() { }

  /// constructor for Polynomial
  ConstPolyRep(const Polynomial<NT>& p, int n) : ss(p) {
    // isolate roots using Sturm Sequences
    I = ss.isolateRoot(n);
    // check whether n-th root exists
    if (I.first == 1 && I.second == 0) {
      std::cerr << "error! root index out of bound." << std::endl;
      abort();
    }
    // test if the root isolated in I is 0:
    if ((I.first == 0)&&(I.second == 0))
        ffVal = 0;
    else
    	ffVal = computeFilteredValue();
  }

  /// constructor for Polynomial
  ConstPolyRep(const Polynomial<NT>& p, const BFInterval& II) : ss(p), I(II) {
    if (ss.numberOfRoots(I.first, I.second) != 1) {
      std::cerr << "error! non-isolating interval." << std::endl;
      abort();
    }
    ffVal = computeFilteredValue(); 
  }

  /// destructor
  ~ConstPolyRep() {}
  //@}
  CORE_MEMORY(ConstPolyRep)
private:
  Sturm<NT> ss; ///< internal Sturm sequences 
  BFInterval I; ///< current interval contains the real value
  filteredFp computeFilteredValue() {
    // refine initial interval to absolute error of 2^(lMSB(k)-54)
    // 	  where k is a lower bound on the root (use Cauchy Lower Bound here).
    //    Hence, the precision we pass to refine should be 54-lMSB(k).
    
    // refine with newton (new method)
    I = ss.newtonRefine(I, 54-(ss.seq[0].CauchyLowerBound()).lMSB().asLong());

    //return I.first.doubleValue(); // NOTE:  This is not quite right!
    			// It should be "centralize" which should set
			// the error bit correctly.
			// E.g., radical(4,2) will print wrongly.
    if ((I.first == 0) && (I.second == 0))	// Checkfor zero value
	    return filteredFp(0);
    BigFloat x = centerize(I.first, I.second);
    double val = x.doubleValue();
    long ee = x.exp()*CHUNK_BIT;
    unsigned long err = ee > 0 ? (x.err() << ee) : (x.err() >> (-ee));
    double max = core_abs(val) + err;
    int ind = ((BigInt(x.err()) << 53) / (x.m() + x.err())).longValue(); 
    return filteredFp(val, max, ind);
  }

protected:
  void initNodeInfo() {
    nodeInfo = new NodeInfo();
    d_e() = ss.seq[0].getTrueDegree(); // return degree of the polynomial
  }
  /// compute sign and MSB
  void computeExactFlags() {

    if ((I.first == 0) && (I.second == 0)) {
	reduceToZero(); 
	return;
    } else if (I.second > 0) {
       uMSB() = I.second.uMSB();
       lMSB() = I.first.lMSB();
       sign() = 1;
    } else { // we know that I.first < 0
       lMSB() = I.second.lMSB();
       uMSB() = I.first.uMSB();
       sign() = -1;
    }
    // length() = 1+ ss.seq[0].length().uMSB();
    measure() = 1+ ss.seq[0].length().uMSB();	// since measure<= length
  
    // compute u25, l25, v2p, v2m, v5p, v5m
    v2p() = v2m() = v5p() = v5m() = 0;
    u25() = 1+ss.seq[0].CauchyUpperBound().uMSB();
    l25() = ceilLg(ss.seq[0].getLeadCoeff());  // assumed coeff is integer!!
  
    // compute high, low, lc, tc
    high() = u25();
    low() = - (ss.seq[0].CauchyLowerBound().lMSB()); // note the use of negative
    lc() = l25();
    tc() = ceilLg(ss.seq[0].getTailCoeff());
  
    // no rational reduction
    if (rationalReduceFlag) ratFlag() = -1;

    flagsComputed() = true;
    appValue()=centerize(I.first, I.second);// set an initial value for appValue
  }
  /// compute approximation value
  void computeApproxValue(const extLong& relPrec, const extLong& absPrec) {
    extLong pr = -lMSB() + relPrec;
    extLong p = pr < absPrec ? pr : absPrec;
    
    // bisection sturm (old method)
    //I = ss.refine(I, p.asLong()+1);
    
    // refine with newton (new method)
    I = ss.newtonRefine(I, p.asLong()+1);
    appValue() = centerize(I.first, I.second);
  }
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
  /// count computes the degree of current node, i.e., d_e().
  /** This is now a misnomer, but historically accurate.
   */
  extLong count();
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
  /// count computes the degree of current node, i.e., d_e().
  /** This is now a misnomer, but historically accurate.
   */
  extLong count();
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
  /// count computes the degree of current node, i.e., d_e().
  /** This is now a misnomer, but historically accurate.
   */
  extLong count();
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

// inline functions
inline int ExprRep::getExactSign() { 
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
inline int ExprRep::getSign() { 
  if (ffVal.isOK())
    return ffVal.sign();
  else
    return getExactSign(); 
}   

// you need force to approximate before call these functions!!
inline BigInt ExprRep::BigIntValue()
{ return getAppValue().BigIntValue(); }

inline BigRat ExprRep::BigRatValue()
{ return getAppValue().BigRatValue(); }

inline BigFloat ExprRep::BigFloatValue()
{ return getAppValue().BigFloatValue(); }
CORE_END_NAMESPACE
#endif
