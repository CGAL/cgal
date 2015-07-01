/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 * You can redistribute it and/or modify it under the terms of the GNU
 * Lesser General Public License as published by the Free Software Foundation,
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
 * File: ExprRep.h
 * Synopsis: Internal Representation of Expr.
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
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#ifndef _CORE_EXPRREP_H_
#define _CORE_EXPRREP_H_

#include <CGAL/CORE/Real.h>
#include <CGAL/CORE/Filter.h>
#include <CGAL/CORE/poly/Sturm.h>

namespace CORE { 

#ifdef CORE_DEBUG_BOUND
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

const extLong EXTLONG_BIG = (1L << 30);
const extLong EXTLONG_SMALL = -(1L << 30);

const double log_5 = log(double(5))/log(double(2));

// Returns the ceil of log_2(5^a).
inline extLong ceilLg5(const extLong & a) {
#if defined (_MSC_VER)
  return (int) ::ceil(log_5 * a.toLong());
#else
  return (int) std::ceil(log_5 * a.toLong());
#endif
}

/// \struct NodeInfo
/// \brief store information of a node
struct NodeInfo {
  Real     appValue;		///< current approximate value
  bool     appComputed;  	///< true if the approx value been computed
  bool     flagsComputed;  	///< true if rootBound parameters have been computed
  extLong  knownPrecision; 	///< Precision achieved by current approx value

#ifdef CORE_DEBUG
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
  /// 2^{-low(E)} is LOWER bound on the moduli of ALL NON_ZERO conjugate of E.
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
  CGAL_CORE_EXPORT NodeInfo();
};//NodeInfo struct

//  forward reference
// class Expr;

/// \class ExprRep
/// \brief The sharable, internal representation of expressions
//  Members: private: int refCount,
//            public:  NodeInfo* nodeInfo,
//                     filteredFp ffVal.
class ExprRep {
public:
  /// \name Constructor and Destructor
  //@{
  /// default constructor
  CGAL_CORE_EXPORT ExprRep();
  /// virtual destructor for this base class
  virtual ~ExprRep() {
    if (nodeInfo != NULL) // This check is only for optimization.
      delete nodeInfo;
  }
  //@}

  /// \name Reference Counting
  //@{
  /// increase reference counter
  void incRef() {
    ++refCount;
  }
  /// decrease reference counter
  void decRef() {
    if (--refCount == 0)
      delete this;
  }
  /// return reference counter
  int getRefCount() const {
    return refCount;
  }
  /// check whether reference counter == 1
  int isUnique() const {
    return refCount == 1;
  }
  //@}

  /// \name Helper Functions
  //@{
  /// Get the approximate value
  CGAL_CORE_EXPORT const Real & getAppValue(const extLong& relPrec = defRelPrec,
                                            const extLong& absPrec = defAbsPrec);
  /// Get the sign.
  int getSign();
  int getExactSign();

  const Real& appValue() const {
    return nodeInfo->appValue;
  }
  Real& appValue() {
    return nodeInfo->appValue;
  }

  const bool& appComputed() const {
    return nodeInfo->appComputed;
  }
  bool& appComputed() {
    return nodeInfo->appComputed;
  }

  const bool& flagsComputed() const {
    return nodeInfo->flagsComputed;
  }
  bool& flagsComputed() {
    return nodeInfo->flagsComputed;
  }

  const extLong& knownPrecision() const {
    return nodeInfo->knownPrecision;
  }
  extLong& knownPrecision() {
    return nodeInfo->knownPrecision;
  }

#ifdef CORE_DEBUG
  const extLong& relPrecision() const {
    return nodeInfo->relPrecision;
  }
  extLong& relPrecision() {
    return nodeInfo->relPrecision;
  }

  const extLong& absPrecision() const {
    return nodeInfo->absPrecision;
  }
  extLong& absPrecision() {
    return nodeInfo->absPrecision;
  }

  const unsigned long& numNodes() const {
    return nodeInfo->numNodes;
  }
  unsigned long& numNodes() {
    return nodeInfo->numNodes;
  }
#endif

  const extLong& d_e() const {
    return nodeInfo->d_e;
  }
  extLong& d_e() {
    return nodeInfo->d_e;
  }

  const bool& visited() const {
    return nodeInfo->visited;
  }
  bool& visited() {
    return nodeInfo->visited;
  }

  const int& sign() const {
    return nodeInfo->sign;
  }
  int& sign() {
    return nodeInfo->sign;
  }

  const extLong& uMSB() const {
    return nodeInfo->uMSB;
  }
  extLong& uMSB() {
    return nodeInfo->uMSB;
  }

  const extLong& lMSB() const {
    return nodeInfo->lMSB;
  }
  extLong& lMSB() {
    return nodeInfo->lMSB;
  }

  //  const extLong& length() const { return nodeInfo->length; }
  //  extLong& length() { return nodeInfo->length; }

  const extLong& measure() const {
    return nodeInfo->measure;
  }
  extLong& measure() {
    return nodeInfo->measure;
  }

  const extLong& high() const {
    return nodeInfo->high;
  }
  extLong& high() {
    return nodeInfo->high;
  }

  const extLong& low() const {
    return nodeInfo->low;
  }
  extLong& low() {
    return nodeInfo->low;
  }

  const extLong& lc() const {
    return nodeInfo->lc;
  }
  extLong& lc() {
    return nodeInfo->lc;
  }

  const extLong& tc() const {
    return nodeInfo->tc;
  }
  extLong& tc() {
    return nodeInfo->tc;
  }

  const extLong& v2p() const {
    return nodeInfo->v2p;
  }
  extLong& v2p() {
    return nodeInfo->v2p;
  }

  const extLong& v2m() const {
    return nodeInfo->v2m;
  }
  extLong& v2m() {
    return nodeInfo->v2m;
  }

  extLong v2() const {
    return v2p()-v2m();
  }

  const extLong& v5p() const {
    return nodeInfo->v5p;
  }
  extLong& v5p() {
    return nodeInfo->v5p;
  }

  const extLong& v5m() const {
    return nodeInfo->v5m;
  }
  extLong& v5m() {
    return nodeInfo->v5m;
  }

  extLong v5() const {
    return v5p()-v5m();
  }

  const extLong& u25() const {
    return nodeInfo->u25;
  }
  extLong& u25() {
    return nodeInfo->u25;
  }

  const extLong& l25() const {
    return nodeInfo->l25;
  }
  extLong& l25() {
    return nodeInfo->l25;
  }

  const int& ratFlag() const {
    return nodeInfo->ratFlag;
  }
  int& ratFlag() {
    return nodeInfo->ratFlag;
  }

  const BigRat* ratValue() const {
    return nodeInfo->ratValue;
  }
  BigRat*& ratValue() {
    return nodeInfo->ratValue;
  }

  /// Get BigFloat
  BigInt BigIntValue();
  BigRat BigRatValue();
  BigFloat BigFloatValue();
  /// represent as a string in decimal value
  // toString() Joaquin Grech 31/5/2003
  std::string toString(long prec=defOutputDigits, bool sci=false) {
    return (getAppValue(defRelPrec, defAbsPrec)).toString(prec,sci);
  }
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
  CGAL_CORE_EXPORT friend std::ostream& operator<<(std::ostream&, ExprRep&);
  //@}

private:
  int refCount;    // reference count

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
  CGAL_CORE_EXPORT extLong computeBound();
  /// driver function to approximate
  CGAL_CORE_EXPORT void approx(const extLong& relPrec, const extLong& absPrec);
  /// compute an approximate value satifying the specified precisions
  virtual void computeApproxValue(const extLong&, const extLong&) = 0;
  /// Test whether the current approx. value satisfies [relPrec, absPrec]
  CGAL_CORE_EXPORT bool withinKnownPrecision(const extLong&, const extLong&);
  //@}

  /// \name Misc Functions
  //@{
  /// reduce current node
  CGAL_CORE_EXPORT void reduceToBigRat(const BigRat&);
  /// reduce current node
  CGAL_CORE_EXPORT void reduceTo(const ExprRep*);
  /// reduce current node to zero
  CGAL_CORE_EXPORT void reduceToZero();
  /// return operator string
  virtual const std::string op() const {
    return "UNKNOWN";
  }
  //@}

  /// \name Degree Bound Functions
  //@{
  /// compute "d_e" based on # of sqrts
  CGAL_CORE_EXPORT extLong degreeBound();
  /// count actually computes the degree bound of current node.
  virtual extLong count() = 0;
  /// reset the flag "visited"
  virtual void clearFlag() = 0;
  //@}
#ifdef CORE_DEBUG
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
  virtual ~ConstRep() {}
  //@}

  /// \name Debug Functions
  //@{
  /// print debug information in list mode
  CGAL_CORE_EXPORT void debugList(int level, int depthLimit) const;
  /// print debug information in tree mode
  CGAL_CORE_EXPORT void debugTree(int level, int indent, int depthLimit) const;
  //@}
protected:
  /// initialize nodeInfo
  CGAL_CORE_EXPORT virtual void initNodeInfo();
  /// return operator in string
  const std::string op() const {
    return "C";
  }
  /// count returns the degree of current node
  //extLong count() { return d_e(); }
  CGAL_CORE_EXPORT extLong count();
  /// clear visited flag
  void clearFlag() {
    visited() = false;
  }
#ifdef CORE_DEBUG
  unsigned long dagSize();
  void fullClearFlag();
#endif
};

/// \class ConstDoubleRep
/// \brief constant node
class ConstDoubleRep : public ConstRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  ConstDoubleRep() {}
  /// constructor for all \c double type
  ConstDoubleRep(double d) {
    ffVal = d;
  }
  /// destructor
  ~ConstDoubleRep() {}
  //@}
  CORE_MEMORY(ConstDoubleRep)
protected:
  /// compute sign and MSB
  CGAL_CORE_EXPORT void computeExactFlags();
  /// compute approximation value
  CGAL_CORE_EXPORT void computeApproxValue(const extLong&, const extLong&);
};

/// \class ConstRealRep
/// \brief constant node
class ConstRealRep : public ConstRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// default constructor
  ConstRealRep() : value(CORE_REAL_ZERO) { }
  /// constructor for all \c Real type
  CGAL_CORE_EXPORT ConstRealRep(const Real &);
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
      core_error("CORE ERROR! root index out of bound",
		      __FILE__, __LINE__, true);
      abort();
    }
    // test if the root isolated in I is 0:
    if ((I.first == 0)&&(I.second == 0))
      ffVal = 0;
    else
      ffVal = computeFilteredValue();	// silly to use a filter here!
    					// since sign is known.
  }

  /// constructor for Polynomial
  ConstPolyRep(const Polynomial<NT>& p, const BFInterval& II)
	  : ss(p), I(II) {
    BFVecInterval v;
    ss.isolateRoots(I.first, I.second, v);
    I = v.front();
    if (v.size() != 1) {
      core_error("CORE ERROR! non-isolating interval",
		      __FILE__, __LINE__, true);
      abort();
    }
    ffVal = computeFilteredValue(); // Chee: this line seems unnecessary
  }

  /// destructor
  ~ConstPolyRep() {}
  //@}
  CORE_MEMORY(ConstPolyRep)

private:
  Sturm<NT> ss; ///< internal Sturm sequences
  BFInterval I; ///< current interval contains the real value
  		// IMPORTANT: I.first and I.second are exact BigFloats
  filteredFp computeFilteredValue() {
    // refine initial interval to absolute error of 2^(lMSB(k)-54) where
    // 	  k is a lower bound on the root (use Cauchy Lower Bound).
    //    Hence, the precision we pass to refine should be 54-lMSB(k).

    // refine with newton (new method)
    // ss.seq[0] could be zero!!
    // I=ss.newtonRefine(I,
    //            54-(ss.seq[0].CauchyLowerBound()).lMSB().asLong());
    extLong lbd = ss.seq[0].CauchyLowerBound().lMSB();

    if (lbd.isTiny())
      I = ss.newtonRefine(I, 54);
    else
      I = ss.newtonRefine(I, 54-lbd.asLong()); // is this necessary?

    //return I.first.doubleValue(); // NOTE:  This is not quite right!
    // It should be "centralized" to set
    // the error bit correctly.
    // E.g., otherwise, radical(4,2) will print wrongly.
    if ((I.first == 0) && (I.second == 0))	// Checkfor zero value
      return filteredFp(0);
    BigFloat x = centerize(I.first, I.second);
    double val = x.doubleValue();
    double max = core_max(core_abs(I.first), core_abs(I.second)).doubleValue();
    int ind = 1;
    /*
    long ee = x.exp()*CHUNK_BIT;
    unsigned long err = ee > 0 ? (x.err() << ee) : (x.err() >> (-ee));
    double max = core_abs(val) + err;
    int ind = longValue((BigInt(x.err()) << 53) / (x.m() + x.err())); 
    */
    return filteredFp(val, max, ind); // Aug 8, 2004, Comment from Chee:
       // I think we should get rid of filters here!  Given the interval I,
       // we either know the sign (I.first >=0) or (I.second <=0)
       // or we don't.  We don't need to compute all the index stuff.
       // In fact, you have lost the sign in the above computation...
       // ALSO, why bother to use filter?
  }//computeFilteredValue

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
    		// ceilLg(BigInt) and ceilLg(Expr) are defined. But if
    		// NT=int, ceilLg(int) is ambiguous!  Added ceilLg(int)
		// under BigInt.h

    // compute high, low, lc, tc
    high() = u25();
    low() = - (ss.seq[0].CauchyLowerBound().lMSB()); // note the use of negative
    lc() = l25();
    tc() = ceilLg(ss.seq[0].getTailCoeff());

    // no rational reduction
    if (rationalReduceFlag)
      ratFlag() = -1;

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
  UnaryOpRep(ExprRep* c) : child(c)  {
    child->incRef();
  }
  /// destructor
  virtual ~UnaryOpRep() {
    child->decRef();
  }
  //@}

  /// \name Debug Functions
  //@{
  /// print debug information in list mode
  CGAL_CORE_EXPORT void debugList(int level, int depthLimit) const;
  /// print debug information in tree mode
  CGAL_CORE_EXPORT void debugTree(int level, int indent, int depthLimit) const;
  //@}
protected:
  ExprRep* child; ///< pointer to its child node
  /// initialize nodeInfo
  CGAL_CORE_EXPORT virtual void initNodeInfo();
  /// clear visited flag
  CGAL_CORE_EXPORT void clearFlag();
#ifdef CORE_DEBUG
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
  NegRep(ExprRep* c) : UnaryOpRep(c) {
    ffVal = - child->ffVal;
  }
  /// destructor
  ~NegRep() {}
  //@}

  CORE_MEMORY(NegRep)
protected:
  /// compute sign and MSB
  CGAL_CORE_EXPORT void computeExactFlags();
  /// compute approximation value
  CGAL_CORE_EXPORT void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const {
    return "Neg";
  }
  /// count computes the degree of current node, i.e., d_e().
  /** This is now a misnomer, but historically accurate.
   */
  CGAL_CORE_EXPORT extLong count();
};

/// \class SqrtRep
/// \brief squartroot operator node
class SqrtRep : public UnaryOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  SqrtRep(ExprRep* c) : UnaryOpRep(c) {
    ffVal = (child->ffVal).sqrt();
  }
  /// destructor
  ~SqrtRep() {}
  //@}

  CORE_MEMORY(SqrtRep)
protected:
  /// compute sign and MSB
  CGAL_CORE_EXPORT void computeExactFlags();
  /// compute approximation value
  CGAL_CORE_EXPORT void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const {
    return "Sqrt";
  }
  /// count computes the degree of current node, i.e., d_e().
  /** This is now a misnomer, but historically accurate.
   */
  CGAL_CORE_EXPORT extLong count();
};

/// \class BinOpRep
/// \brief binary operator node
class BinOpRep : public ExprRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  BinOpRep(ExprRep* f, ExprRep* s) : first(f), second(s) {
    first->incRef();
    second->incRef();
  }
  /// destructor
  virtual ~BinOpRep() {
    first->decRef();
    second->decRef();
  }
  //@}

  /// \name Debug Functions
  //@{
  /// print debug information in list mode
  CGAL_CORE_EXPORT void debugList(int level, int depthLimit) const;
  /// print debug information in tree mode
  CGAL_CORE_EXPORT void debugTree(int level, int indent, int depthLimit) const;
  //@}
protected:
  ExprRep* first;  ///< first operand
  ExprRep* second; ///< second operand

  /// initialize nodeInfo
  CGAL_CORE_EXPORT virtual void initNodeInfo();
  /// clear visited flags
  CGAL_CORE_EXPORT void clearFlag();
  /// count computes the degree of current node, i.e., d_e().
  /** This is now a misnomer, but historically accurate.
   */
  CGAL_CORE_EXPORT extLong count();
#ifdef CORE_DEBUG
  CGAL_CORE_EXPORT unsigned long dagSize();
  CGAL_CORE_EXPORT void fullClearFlag();
#endif
};

/// \struct Add
/// \brief "functor" class used as parameter to AddSubRep<>
struct Add {
  /// name
  CGAL_CORE_EXPORT static const char* name;

  /// unary operator
  template <class T>
  const T& operator()(const T& t) const {
    return t;
  }

  /// binary operator
  template <class T>
  T operator()(const T& a, const T& b) const {
    return a+b;
  }
};

/// \struct Sub
/// \brief "functor" class used as parameter to AddSubRep<>
struct Sub {
  /// name
  CGAL_CORE_EXPORT static const char* name;

  /// unary operator
  template <class T>
  T operator()(const T& t) const {
    return -t;
  }

  /// binary operator
  template <class T>
  T operator()(const T& a, const T& b) const {
    return a-b;
  }
};

/// \class AddSubRep
/// \brief template class where operator is supposed to be Add or Sub
template <class Operator>
class AddSubRep : public BinOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  AddSubRep(ExprRep* f, ExprRep* s) : BinOpRep(f, s) {
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
  const std::string op() const {
    return Operator::name;
  }
private:
  static Operator Op;
};//AddSubRep class

template <class Operator>
Operator AddSubRep<Operator>::Op;

  /// AddSubRep<Op>::computeExactFlags()
  ///     This function is the heart of Expr class,
  ///     and hence the heart of Core Library!
  /// Here is where we use the root bounds.
template <class Operator>
void AddSubRep<Operator>::computeExactFlags() {
  if (!first->flagsComputed())
    first->computeExactFlags();
  if (!second->flagsComputed())
    second->computeExactFlags();

  int sf = first->sign();
  int ss = second->sign();

  if ((sf == 0) && (ss == 0)) { // the node is zero
    reduceToZero();
    return;
  } else if (sf == 0) { // first operand is zero
    reduceTo(second);
    sign() = Op(ss);
    appValue() = Op(appValue());
    if (rationalReduceFlag && ratFlag() > 0)
      *(ratValue()) = Op(*(ratValue()));
    return;
  } else if (ss == 0) { // second operand is zero
    reduceTo(first);
    return;
  }
  // rational node
  if (rationalReduceFlag) {
    if (first->ratFlag() > 0 && second->ratFlag() > 0) {
      BigRat val=Op(*(first->ratValue()), *(second->ratValue()));
      reduceToBigRat(val);
      ratFlag() = first->ratFlag() + second->ratFlag();
      return;
    } else
      ratFlag() = -1;
  }

  // neither operand is zero
  extLong df    = first->d_e();
  extLong ds    = second->d_e();
  // extLong md    = df < ds ? df : ds;
  // extLong l1    = first->length();
  // extLong l2    = second->length();
  extLong m1    = first->measure();
  extLong m2    = second->measure();

  // length() = df * l2 + ds * l1 + d_e() + md;
  measure() = m1 * ds + m2 * df + d_e();

  // BFMSS[2,5] bound.
  v2p() = core_min(first->v2p() + second->v2m(),
		  first->v2m() + second->v2p());
  v2m() = first->v2m() + second->v2m();
  v5p() = core_min(first->v5p() + second->v5m(),
		  first->v5m() + second->v5p());
  v5m() = first->v5m() + second->v5m();

  if (v2p().isInfty() || v5p().isInfty())
    u25() = CORE_INFTY;
  else
    u25() = EXTLONG_ONE + core_max(first->v2p() + second->v2m()
	    - v2p() + ceilLg5(first->v5p() + second->v5m() - v5p())
            + first->u25() + second->l25(),
                                   first->v2m() + second->v2p() - v2p()
            + ceilLg5(first->v5m() + second->v5p() - v5p())
            + first->l25() + second->u25());
  l25() = first->l25() + second->l25();

  lc() = ds * first->lc() + df * second->lc();
  tc() = measure();

  high() = core_max(first->high(),second->high())+EXTLONG_ONE;
  // The following is a subset of the minimization in computeBound().
  low() = core_min(measure(), (d_e()-EXTLONG_ONE)*high() + lc());

  extLong lf = first->lMSB();
  extLong ls = second->lMSB();
  extLong uf = first->uMSB();
  extLong us = second->uMSB();

  extLong l  = core_max(lf, ls);
  extLong u  = core_max(uf, us);

#ifdef CORE_TRACE
  std::cout << "INSIDE Add/sub Rep: " << std::endl;
#endif

  if (Op(sf, ss) != 0) {     // can't possibly cancel out
#ifdef CORE_TRACE
    std::cout << "Add/sub Rep:  Op(sf, ss) non-zero" << std::endl;
#endif

    uMSB() = u + EXTLONG_ONE;
    lMSB() = l;            // lMSB = core_min(lf, ls)+1 better
    sign() = sf;
  } else {               // might cancel out
#ifdef CORE_TRACE
    std::cout << "Add/sub Rep:  Op(sf, ss) zero" << std::endl;
#endif

    uMSB() = u + EXTLONG_ONE;
    uMSB() = u;
    if (lf >= us + EXTLONG_TWO) {// one is at least 1 order of magnitude larger
#ifdef CORE_TRACE
      std::cout << "Add/sub Rep:  Can't cancel" << std::endl;
#endif

      lMSB() = lf - EXTLONG_ONE;     // can't possibly cancel out
      sign() = sf;
    } else if (ls >= uf + EXTLONG_TWO) {
#ifdef CORE_TRACE
      std::cout << "Add/sub Rep:  Can't cancel" << std::endl;
#endif

      lMSB() = ls - EXTLONG_ONE;
      sign() = Op(ss);
    } else if (ffVal.isOK()) {// begin filter computation
#ifdef CORE_TRACE
      std::cout << "Add/sub Rep:  filter used" << std::endl;
#endif
#ifdef CORE_DEBUG_FILTER
      std::cout << "call filter in " << op() << "Rep" << std::endl;
#endif
      sign() = ffVal.sign();
      lMSB() = ffVal.lMSB();
      uMSB() = ffVal.uMSB();
    } else {			// about the same size, might cancel out
#ifdef CORE_TRACE
      std::cout << "Add/sub Rep:  iteration start" << std::endl;
#endif

      extLong lowBound = computeBound();
      /* Zilin 06/11/2003
       * as BFMSS[2] might be a negative number, lowBound can be negative.
       * In this case, we just set it to 1 since we need at least one bit
       * to get the sign.  In the future, we may need to improve this.
       */
      if (lowBound <= EXTLONG_ZERO)
        lowBound = EXTLONG_ONE;

      if (!progressiveEvalFlag) {
        // convert the absolute error requirement "lowBound" to
        // a relative error requirement "ur", s.t.
        //    |x|*2^(-ur) <= 2^(-lowBound).
        // ==> r >= a + lg(x) >= a + (uMSB + 1);
        //	    extLong  rf = lowBound + (uf + 1);
        //	    extLong  rs = lowBound + (us + 1);
        //	    first->approx(rf, CORE_INFTY);
        //	    second->approx(rs, CORE_INFTY);
        // Chen: considering the uMSB is also an approximate bound.
        // we choose to use absolute precision up-front.
        Real newValue = Op(first->getAppValue(CORE_INFTY,
				lowBound + EXTLONG_ONE),
                           second->getAppValue(CORE_INFTY,
				   lowBound + EXTLONG_ONE));
        if (!newValue.isZeroIn()) { // Op(first, second) != 0
          lMSB() = newValue.lMSB();
          uMSB() = newValue.uMSB();   // chen: to get tighers value.
          sign() = newValue.sign();
        } else if (lowBound.isInfty()) {//check if rootbound is too big
          core_error("AddSubRep:root bound has exceeded the maximum size\n \
            but we still cannot decide zero.\n", __FILE__, __LINE__, false);
        } else {               // Op(first, second) == 0
          lMSB() = CORE_negInfty;
          sign() = 0;
        }
      } else {  // else do progressive evaluation
#ifdef CORE_TRACE
        std::cout << "Add/sub Rep:  progressive eval" << std::endl;
#endif
        // Oct 30, 2002: fixed a bug here!  Old versions used relative
        // precision bounds, but one should absolute precision for addition!
        // Moreover, this is much more efficient.

        // ua is the upper bound on the absolute precision in our iteration
	// Chee, Aug 8, 2004: it is important that ua be strictly
	//     larger than lowBound AND the defaultInitialProgressivePrec,
	//     so that we do at least one iteration of the for-loop. So:
	// i is the variable for iteration.
        extLong i = core_min(defInitialProgressivePrec, lowBound.asLong());
        extLong ua = lowBound.asLong() + EXTLONG_ONE;
        //   NOTE: ua is allowed to be CORE_INFTY
	
#ifdef CORE_DEBUG_BOUND
        std::cout << "DebugBound:" << "ua = " << ua << std::endl;
#endif
        // We initially set the lMSB and sign as if the value is zero:
        lMSB() = CORE_negInfty;
        sign() = 0;

        EscapePrecFlag = 0;	// reset the Escape Flag

        // Now we try to determine the real lMSB and sign,
        // in case it is not really zero:
#ifdef CORE_TRACE
        std::cout << "Upper bound (ua) for iteration is " << ua << std::endl;
	std::cout << "Starting iteration at i = " << i << std::endl;
#endif

        for ( ; i<ua; i*=EXTLONG_TWO) {
          // relative bits = i
	  //
	  // PROBLEM WITH NEXT LINE: you must ensure that
	  // first and second are represented by BigFloats...
	  //
          Real newValue = Op(first->getAppValue(CORE_INFTY, i),
                             second->getAppValue(CORE_INFTY, i));

#ifdef CORE_TRACE
	  if (newValue.getRep().ID() == REAL_BIGFLOAT) 
	  std::cout << "BigFloat! newValue->rep->ID() = "
		  << newValue.getRep().ID() << std::endl;
	  else 
	  std::cout << "ERROR, Not BigFloat! newValue->rep->ID() ="
		  << newValue.getRep().ID() << std::endl;
	  std::cout << "newValue = Op(first,second) = "
		  << newValue << std::endl;
	  std::cout << "first:appVal, appComputed, knownPrec, sign ="
		  << first->appValue() << ","
		  << first->appComputed() << ","
		  << first->knownPrecision() << ","
		  << first->sign() << std::endl;
	  std::cout << "second:appVal, appComputed, knownPrec, sign ="
		  << second->appValue() << ","
		  << second->appComputed() << ","
		  << second->knownPrecision() << ","
		  << second->sign() << std::endl;
#endif
          if (!newValue.isZeroIn()) {   // Op(first, second) != 0
            lMSB() = newValue.lMSB();
            uMSB() = newValue.uMSB();
            sign() = newValue.sign();
#ifdef CORE_DEBUG_BOUND
            std::cout << "DebugBound(Exit Loop): " << "i=" << i << std::endl;
#endif
#ifdef CORE_TRACE
	    std::cout << "Zero is not in, lMSB() = " << lMSB()
		    << ", uMSB() = " << uMSB()
		    << ", sign() = " << sign() << std::endl;
	    std::cout << "newValue = " << newValue << std::endl;
#endif

            break; // assert -- this must happen in the loop if nonzero!
          }
          //8/9/01, Chee: implement escape precision here:
          if (i> EscapePrec) {
            EscapePrecFlag = -i.asLong();//negative means EscapePrec is used
	    core_error("Escape precision triggered at",
            		 __FILE__, __LINE__, false);
            if (EscapePrecWarning)
              std::cout<< "Escape Precision triggered at "
		      << EscapePrec << " bits" << std::endl;
#ifdef CORE_DEBUG
            std::cout << "EscapePrecFlags=" << EscapePrecFlag << std::endl;
            std::cout << "ua =" << ua  << ",lowBound=" << lowBound << std::endl;
#endif
            break;
          }// if
        }// for (long i=1...)

#ifdef CORE_DEBUG_BOUND
        rootBoundHitCounter++;
#endif

        if (sign() == 0 && ua .isInfty()) {
          core_error("AddSubRep: root bound has exceeded the maximum size\n \
            but we still cannot decide zero.\n", __FILE__, __LINE__, true);
        } // if (sign == 0 && ua .isInfty())
      }// else do progressive
    }
  }
  flagsComputed() = true;
}// AddSubRep::computeExactFlags

template <class Operator>
void AddSubRep<Operator>::computeApproxValue(const extLong& relPrec,
    const extLong& absPrec) {
  // Nov 13, 2002: added the analog of "reduceTo(first)" and "reduceTo(second)"
  //  that is found in computeExactFlags.  This is more efficient, but
  //  it also removes a NaN warning in subsequent logic!
  //  E.g., if first=0, then first->uMSB and first->lMSB are -infty, and
  //  subtracting them creates NaN.  Chee and Zilin.
  if (first->sign() == 0) {
    appValue() = Op(second->getAppValue(relPrec, absPrec));
    return;
  }
  if (second->sign() == 0) {
    appValue() = first->getAppValue(relPrec, absPrec);
    return;
  }
  if (lMSB() < EXTLONG_BIG && lMSB() > EXTLONG_SMALL) {
    extLong rf = first->uMSB()-lMSB()+relPrec+EXTLONG_FOUR;  // 2 better
    if (rf < EXTLONG_ZERO)
      rf = EXTLONG_ZERO;  // from Koji's thesis P63: Proposition 26
    extLong rs = second->uMSB()-lMSB()+relPrec+EXTLONG_FOUR; // 2 better
    if (rs < EXTLONG_ZERO)
      rs = EXTLONG_ZERO;  // from Koji's thesis P63: Proposition 26
    extLong  a  = absPrec + EXTLONG_THREE;                      // 1 better
    appValue() = Op(first->getAppValue(rf, a), second->getAppValue(rs, a));
  } else {
    std::cerr << "lMSB = " << lMSB() << std::endl; // should be in core_error
    core_error("CORE WARNING: a huge lMSB in AddSubRep",
	 	__FILE__, __LINE__, false);
  }
}

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
  MultRep(ExprRep* f, ExprRep* s) : BinOpRep(f, s) {
    ffVal = first->ffVal * second->ffVal;
  }
  /// destructor
  ~MultRep() {}
  //@}
  
  CORE_MEMORY(MultRep)
  protected:
  /// compute sign and MSB
  CGAL_CORE_EXPORT void computeExactFlags();
  /// compute approximation value
  CGAL_CORE_EXPORT void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const {
    return "*";
  }
};

/// \class DivRep
/// \brief division operator node
class DivRep : public BinOpRep {
public:
  /// \name Constructors and Destructor
  //@{
  /// constructor
  DivRep(ExprRep* f, ExprRep* s) : BinOpRep(f, s) {
    ffVal = first->ffVal / second->ffVal;
  }
  /// destructor
  ~DivRep() {}
  //@}

  CORE_MEMORY(DivRep)
protected:
  /// compute sign and MSB
  CGAL_CORE_EXPORT void computeExactFlags();
  /// compute approximation value
  CGAL_CORE_EXPORT void computeApproxValue(const extLong&, const extLong&);
  /// return operator in string
  const std::string op() const {
    return "/";
  }
};

// inline functions
inline int ExprRep::getExactSign() {
  if (!nodeInfo)
    initNodeInfo();

  if (!flagsComputed()) {
    degreeBound();
#ifdef CORE_DEBUG
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
inline BigInt ExprRep::BigIntValue() {
  return getAppValue().BigIntValue();
}

inline BigRat ExprRep::BigRatValue() {
  return getAppValue().BigRatValue();
}

inline BigFloat ExprRep::BigFloatValue() {
  return getAppValue().BigFloatValue();
}

} //namespace CORE
#endif // _CORE_EXPRREP_H_
