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
 * File: Expr.cpp
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
 * $URL$
 * $Id$
 ***************************************************************************/

#include <CGAL/CORE/Expr.h>
#include <cmath>

namespace CORE { 

#ifdef CORE_DEBUG_BOUND
unsigned int BFMSS_counter = 0;
unsigned int BFMSS_only_counter = 0;
unsigned int Measure_counter = 0;
unsigned int Measure_only_counter = 0;
unsigned int Cauchy_counter = 0;
unsigned int Cauchy_only_counter = 0;
unsigned int LiYap_counter = 0;
unsigned int LiYap_only_counter = 0;
unsigned int rootBoundHitCounter = 0;
unsigned int computeBoundCallsCounter = 0;
#endif

const char* Add::name = "+";
const char* Sub::name = "-";

/********************************************************
 *  class Expr
 ********************************************************/
const Expr& Expr::getZero() {
  static Expr Zero(0);
  return Zero;
}
const Expr& Expr::getOne() {
  static Expr One(1);
  return One;
}

// computes an interval comprising a pair of doubles
// Note:
//
// This function returns are two consecutive representable binary
// IEEE double values whichs contain the real value, but when you print out
// them, you might be confused by the decimal represention due to round.
//
void Expr::doubleInterval(double & lb, double & ub) const {
  double d = doubleValue();
  if (! CGAL_CORE_finite(d)) {	// if overflow, underflow or NaN
    lb = ub = d;
    return;
  }
  int sign = ((* this) -Expr(d)).sign();
  // Seems like doubleValue() always give a lower bound,
  // 	so sign = 0 or 1 (never -1).
  //std::cout << "Sign = " << sign << std::endl;
  if (sign == 0) {
    lb = ub = d;
    return;
  }
  int exp;
  frexp(d, & exp);  	// get the exponent of d
  exp--;		// the exp from frexp satisfies
  //     2^{exp-1} <= d < 2^{exp}
  // But, we want exp to satisfy
  //     2^{exp} <= d < 2^{exp+1}
  if (sign > 0) {
    lb = d;
    ub = d + ldexp(1.0, -52+exp);
    return;
  } else {
    ub = d;
    lb = d - ldexp(1.0, -52+exp);
    return;
  }
}

// floor(e, sub) returns the floor(e), and puts the
//      remainder into sub.
BigInt floor(const Expr& e, Expr &sub) {
  if (e==0) {
	  return 0;
  }
  BigInt f = e.approx(CORE_INFTY, 2).BigIntValue();
  sub = e-f;
  // Adjustment
  if (sub<0)
    ++sub, --f;
  if (sub>=1)
    --sub, ++f;
  assert(sub >=0 && sub<1); // got an assertion error? (Chee 3/24/04)
  return f;
}

// Chenli: implemented algorithm from Goldberg's article.
// 7/01: Thanks to Serge Pashkov for fixing infinite loop when n=0.
Expr pow(const Expr& e, unsigned long n) {
  if (n == 0)
    return Expr(1);
  else if (n == 1)
    return e;
  else {
    Expr x = e;
    while ((n % 2) == 0) { // n is even
      x *= x;
      n >>= 1;
    }
    Expr u = x;
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

NodeInfo::NodeInfo() : appValue(CORE_REAL_ZERO), appComputed(false),
    flagsComputed(false), knownPrecision(CORE_negInfty),
#ifdef CORE_DEBUG
    relPrecision(EXTLONG_ZERO), absPrecision(CORE_negInfty), numNodes(0),
#endif
    // Most of the following data members don't need to be
    // initialized here.
    d_e(EXTLONG_ZERO), visited(false), sign(0),
    uMSB(CORE_negInfty), lMSB(CORE_negInfty),
    // length(0),
    measure(EXTLONG_ZERO), high(EXTLONG_ZERO), low(EXTLONG_ONE),
    lc(EXTLONG_ZERO), tc(EXTLONG_ZERO),
    v2p(EXTLONG_ZERO), v2m(EXTLONG_ZERO),
    v5p(EXTLONG_ZERO), v5m(EXTLONG_ZERO),
    u25(EXTLONG_ZERO), l25(EXTLONG_ZERO),
    ratFlag(0), ratValue(NULL) { }

/********************************************************
 *  class ExprRep
 ********************************************************/
//  constructor
ExprRep::ExprRep() : refCount(1), nodeInfo(NULL), ffVal(0.0) { }

// Computes the root bit bound of the expression.
// In effect, computeBound() returns the current value of low.
extLong ExprRep::computeBound() {
  extLong measureBd = measure();
  // extLong cauchyBd = length();
  extLong ourBd = (d_e() - EXTLONG_ONE) * high() + lc();
  // BFMSS[2,5] bound.
  extLong bfmsskBd;
  if (v2p().isInfty() || v2m().isInfty())
    bfmsskBd = CORE_INFTY;
  else
    bfmsskBd = l25() + u25() * (d_e() - EXTLONG_ONE) - v2() - ceilLg5(v5());

  // since we might compute \infty - \infty for this bound
  if (bfmsskBd.isNaN())
    bfmsskBd = CORE_INFTY;

  extLong bd = core_min(measureBd,
                        // core_min(cauchyBd,
                        core_min(bfmsskBd, ourBd));
#ifdef CORE_SHOW_BOUNDS
    std::cout << "Bounds (" << measureBd <<
                    "," << bfmsskBd << ", " << ourBd << "),  ";
    std::cout << "MIN = " << bd << std::endl;
    std::cout << "d_e= " << d_e() << std::endl;
#endif

#ifdef CORE_DEBUG_BOUND
  // Some statistics about which one is/are the winner[s].
  computeBoundCallsCounter++;
  int number_of_winners = 0;
  std::cerr << " New contest " << std::endl;
  if (bd == bfmsskBd) {
    BFMSS_counter++;
    number_of_winners++;
    std::cerr << " BFMSS is the winner " << std::endl;
  }
  if (bd == measureBd) {
    Measure_counter++;
    number_of_winners++;
    std::cerr << " measureBd is the winner " << std::endl;
  }
  /*  if (bd == cauchyBd) {
      Cauchy_counter++;
      number_of_winners++;
      std::cerr << " cauchyBd is the winner " << std::endl;
    }
   */
  if (bd == ourBd) {
    LiYap_counter++;
    number_of_winners++;
    std::cerr << " ourBd is the winner " << std::endl;
  }

  assert(number_of_winners >= 1);

  if (number_of_winners == 1) {
    if (bd == bfmsskBd) {
      BFMSS_only_counter++;
      std::cerr << " BFMSSBd is the only winner " << std::endl;
    } else if (bd == measureBd) {
      Measure_only_counter++;
      std::cerr << " measureBd is the only winner " << std::endl;
    }
    /* else if (bd == cauchyBd) {
      Cauchy_only_counter++;
      std::cerr << " cauchyBd is the only winner " << std::endl;
    } */
    else if (bd == ourBd) {
      LiYap_only_counter++;
      std::cerr << " ourBd is the only winner " << std::endl;
    }
  }
#endif

  return bd;
}//computeBound()

void ExprRep::reduceToBigRat(const BigRat& rat) {
  Real value(rat);

  //appValue() = value;
  appComputed() = false; // since appValue is not assigned until approx() is called
  flagsComputed() = true;
  knownPrecision() = CORE_negInfty;

#ifdef CORE_DEBUG
  relPrecision() = EXTLONG_ZERO;
  absPrecision() = CORE_negInfty;
  //numNodes() = numNodes();
#endif

  d_e() = EXTLONG_ONE;
  //visited() = e->visited();
  sign() = value.sign();
  uMSB() = value.MSB();
  lMSB() = value.MSB();
  // length() = value.length(); 	// fixed? original = 1
  measure() = value.height();		// measure <= height for rational value

  // BFMSS[2,5] bound.
  value.ULV_E(u25(), l25(), v2p(), v2m(), v5p(), v5m());

  extLong u_e = u25() + v2p();
  extLong l_e = l25() + v2m();

  u_e = u_e + ceilLg5(v5p());
  l_e = l_e + ceilLg5(v5m());

  if (l_e == EXTLONG_ZERO) {           // no divisions introduced
    high() = u_e;
    low() = EXTLONG_ONE - u_e; // - (u_e - 1)
  } else {
    high() = u_e - l_e + EXTLONG_ONE;
    low() = 2 - high();
  }

  lc() = l_e;
  tc() = u_e;

  if (ratValue() == NULL)
    ratValue() = new BigRat(rat);
  else
    *(ratValue()) = rat;
}

// This only copies the current information of the argument e to
// 	*this ExprRep.
void ExprRep::reduceTo(const ExprRep *e) {
  if (e->appComputed()) {
    appValue() = e->appValue();
    appComputed() = true;
    flagsComputed() = true;
    knownPrecision() = e->knownPrecision();
#ifdef CORE_DEBUG
    relPrecision() = e->relPrecision();
    absPrecision() = e->absPrecision();
    numNodes() = e->numNodes();
#endif

  }
  d_e() = e->d_e();
  //visited() = e->visited();
  sign() = e->sign();
  uMSB() = e->uMSB();
  lMSB() = e->lMSB();
  // length() = e->length(); 	// fixed? original = 1
  measure() = e->measure();

  // BFMSS[2,5] bound.
  u25() = e->u25();
  l25() = e->l25();
  v2p() = e->v2p();
  v2m() = e->v2m();
  v5p() = e->v5p();
  v5m() = e->v5m();

  high() = e->high();
  low() = e->low();		// fixed? original = 0
  lc() = e->lc();
  tc() = e->tc();

  // Chee (Mar 23, 2004), Notes on ratFlag():
  // ===============================================================
  // For more information on the use of this flag, see progs/pentagon.
  // This is an integer valued member of the NodeInfo class.
  // Its value is used to determine whether
  // we can ``reduce'' an Expression to a single node containing
  // a BigRat value.  This reduction is done if the global variable
  // rationalReduceFlag=true.  The default value is false.
  // This is the intepretation of ratFlag:
  //	ratFlag < 0 means irrational
  //	ratFlag = 0 means not initialized
  //	ratFlag > 0 means rational
  // Currently, ratFlag>0 is an upper bound on the size of the expression,
  // since we recursively compute
  // 		ratFlag(v) = ratFlag(v.lchild)+ratFlag(v.rchild) + 1.
  // PROPOSAL: if ratFlag() > RAT_REDUCE_THRESHHOLD
  // 	then we automatically do a reduction.  We must determine
  // 	an empirical value for RAT_REDUCE_THRESHOLD

  if (rationalReduceFlag) {
    ratFlag() = e->ratFlag();

    if (e->ratFlag() > 0 && e->ratValue() != NULL) {
      ratFlag() ++;
      if (ratValue() == NULL)
        ratValue() = new BigRat(*(e->ratValue()));
      else
        *(ratValue()) = *(e->ratValue());
    } else
      ratFlag() = -1;
  }
}

void ExprRep::reduceToZero() {
  appValue() = CORE_REAL_ZERO;
  appComputed() = true;
  flagsComputed() = true;
  knownPrecision() = CORE_negInfty;
#ifdef CORE_DEBUG
  relPrecision() = EXTLONG_ZERO;
  absPrecision() = CORE_negInfty;
  //  numNodes() = 0;
#endif

  d_e() = EXTLONG_ONE;
  visited() = false;
  sign() = 0;
  uMSB() = CORE_negInfty;
  lMSB() = CORE_negInfty;
  // length() = 0; 	// fixed? original = 1
  measure() = EXTLONG_ZERO;

  // BFMSS[2,5] bound.
  u25() = l25() = v2p() = v2m() = v5p() = v5m() = EXTLONG_ZERO;

  low() = EXTLONG_ONE;		// fixed? original = 0
  high() = lc() = tc() = EXTLONG_ZERO;

  if (rationalReduceFlag) {
    if (ratFlag() > 0) {
      ratFlag() ++;
      if (ratValue() == NULL)
        ratValue() = new BigRat(0);
      else
        *(ratValue()) = 0;
    } else
      ratFlag() = 1;
  }
}

////////////////////////////////////////////////////////////
//  Test whether the current approximate value satisfies
//  the composite precision requirements [relPrec, absPrec].
////////////////////////////////////////////////////////////

bool ExprRep::withinKnownPrecision(const extLong& relPrec,
                                   const extLong& absPrec) {
  if (appComputed()) { // an approximate value has been evaluated.
    if (appValue().isExact()) {
      return true;
    } else { // approximation has certain error.
      // decide to which position it is required to compute correctly.
      extLong required = core_max(-absPrec, appValue().lMSB()-relPrec);
      // see whether the existing error is smaller than the requirement.
      return (knownPrecision() <= required);
    }
  } else
    return false;
}//withinKnownPrecision(a, r)

// approximate the expression to certain precisions when
// necessary (either no approximate value available or
// the existing one is not accurate enough).
void ExprRep::approx(const extLong& relPrec = defRelPrec,
                     const extLong& absPrec = defAbsPrec) {
  if (!getSign())
    return; // if it is exactly zero...

  // NOTE: The Filter might give a precise enough approximation already.
  if (!getExactSign())
    return;

  if (!appComputed() || (!withinKnownPrecision(relPrec, absPrec))) {
    // it's necessary to re-evaluate.
    // to avoid huge lMSB which would cause long time and problems.

    // if it is a rational node
    if (rationalReduceFlag && ratFlag() > 0 && ratValue() != NULL)
      appValue() = Real(*(ratValue())).approx(relPrec, absPrec); //< shouldn't
                         // this case be done by computeApproxValue()?
    else
      computeApproxValue(relPrec, absPrec);

    // update flags
    appComputed() = true;
    knownPrecision() = appValue().clLgErr();
#ifdef CORE_DEBUG
    if (relPrecision() < relPrec)
      relPrecision() = relPrec;
    if (absPrecision() < absPrec)
      absPrecision() = absPrec;
#endif

  }
}

// return an approximate value to certain precision.
const Real& ExprRep::getAppValue(const extLong& relPrec,
		const extLong& absPrec) {
  if (getSign()) {
    approx(relPrec, absPrec);
    return appValue();
  } else
    return CORE_REAL_ZERO;
}

std::ostream& operator<<(std::ostream& o, ExprRep& rep) {
  if (rep.getSign()) {
    rep.approx(defRelPrec, defAbsPrec);
    o << rep.appValue();
  } else {
    o << "0";
  }
  return o;
}

// Chee, Zilin: July 17, 2002
//  Original algorithm is wrongly implemented, and can take time
//	 exponential in the size of the dag.
//
//  METHOD:
//	Inductively assume that all "visited" flags are false.
//	This calls for a reimplementation of "count()" and "clearFlag()".
//	Actually, we did not have to fix the count() function.
//
//  (1) First recursively compute d_e for each node
//		by calling the count() function.
//	Important thing is count() will turn the "visited" flags
//		to be true, so that there is no double counting.
//	Furthermore, if d_e had already been computed, the
//		arithmetic for d_e can be avoided (in this case,
//		it is only the setting of "visited" flag that we
//		are interested in!
//  (2) At the end of count(), we have set all reachable nodes
//		to "visited", and their d_e have been computed.
//  (3) Now, call clearFlag() to recursively clear all reachable
//		nodes.  NOTE THAT PREVIOUSLY, clearFlag() was called
//		first!  This obvious is wrong

extLong ExprRep::degreeBound() {
  if (d_e() == EXTLONG_ONE) // no radical nodes below
    return EXTLONG_ONE;
  count();
  clearFlag();
  return d_e();
}
//  class ConstRealRep
//  constructor
ConstRealRep::ConstRealRep(const Real & r) : value(r) {
  if (!value.isExact()) {
    // clone the BigFloat and set its error to zero.
    value = value.BigFloatValue().makeExact();
  }
  ffVal = filteredFp(value);
}

// initialize nodeInfo
void ConstRep::initNodeInfo() {
  nodeInfo = new NodeInfo();
  d_e() = EXTLONG_ONE;
}
void UnaryOpRep::initNodeInfo() {
  if (child->nodeInfo == NULL)
    child->initNodeInfo();
  nodeInfo = new NodeInfo();
}
void BinOpRep::initNodeInfo() {
  if (first->nodeInfo == NULL)
    first->initNodeInfo();
  if (second->nodeInfo == NULL)
    second->initNodeInfo();
  nodeInfo = new NodeInfo();
}

#ifdef CORE_DEBUG
unsigned long ConstRep::dagSize() {
  if (!visited()) {
    visited() = true;
    numNodes() = 1;
  } else
    numNodes() = 0;
  return numNodes();
}

unsigned long UnaryOpRep::dagSize() {
  if (!visited()) {
    visited() = true;
    numNodes() = child->dagSize() + 1;
  } else
    numNodes() = 0;
  return numNodes();
}

unsigned long BinOpRep::dagSize() {
  if (!visited()) {
    visited() = true;
    numNodes() = first->dagSize() + second->dagSize() + 1;
  } else
    numNodes() = 0;
  return numNodes();
}

void ConstRep::fullClearFlag() {
  if (visited())
    visited() = false;
}

void UnaryOpRep::fullClearFlag() {
  if (visited()) {
    child->fullClearFlag();
    visited() = false;
  }
}

void BinOpRep::fullClearFlag() {
  if (visited()) {
    first->fullClearFlag();
    second->fullClearFlag();
    visited() = false;
  }
}
#endif

//
// clear visited flag
//
/* see Expr.h
  void ConstRep::clearFlag()
  { visited = false; }
*/
void UnaryOpRep::clearFlag() {
  if (d_e() == EXTLONG_ONE)
    return; // no radicals below.
  if (visited()) {
    visited() = false;
    child->clearFlag();
  }
}
//  class BinOpRep
void BinOpRep::clearFlag() {
  if (d_e() == EXTLONG_ONE)
    return; // rational below
  if (visited()) {
    visited() = false;
    first->clearFlag();
    second->clearFlag();
  }
}

//
// count # of squareroot
//
extLong ConstRep::count() {
  if (visited())
    return EXTLONG_ONE;
  visited() = true;
  return d_e();
}

extLong NegRep::count() {
  if (d_e() == EXTLONG_ONE)
    return EXTLONG_ONE;
  if (visited())
    return EXTLONG_ONE;
  visited() = true;
  d_e() = child->count();
  return d_e();
}

extLong SqrtRep::count() {
  if (d_e() == EXTLONG_ONE)
    return EXTLONG_ONE;
  if (visited())
    return EXTLONG_ONE;
  visited() = true;
  d_e() = child->count() * EXTLONG_TWO;
  return d_e();
}

extLong BinOpRep::count() {
  if (d_e() == EXTLONG_ONE)
    return EXTLONG_ONE;
  if (visited())
    return EXTLONG_ONE;
  visited() = true;
  d_e() = first->count() * second->count();
  return d_e();
}

//
// compute exact flags functions
//
//  exact value

void computeExactFlags_temp(ConstRep* t, const Real &value) {
  // Chen Li: the following is incorrect:
  //    uMSB = lMSB = value.MSB();
  // because the value could be a BigFloat which is an interval.
  if (value.isExact()) {
    t->uMSB() = t->lMSB() = value.MSB();
  } else {
    t->uMSB() = value.uMSB();
    t->lMSB() = value.lMSB();
    core_error("Leaves in DAG is not exact!", __FILE__, __LINE__, true);
  }

  t->sign() = value.sign();
  // t->length() = value.length();
  t->measure() = value.height(); // for rationals and integers,
  // measure = height.

  // BFMSS[2,5] bound.
  value.ULV_E(t->u25(), t->l25(), t->v2p(), t->v2m(), t->v5p(), t->v5m());

  // The original BFMSS parameters can be set from the BFMSS[2,5] parameters.
  // Here we just need them locally.
  extLong u_e = t->u25() + t->v2p() + ceilLg5(t->v5p());
  extLong l_e = t->l25() + t->v2m() + ceilLg5(t->v5m());

#ifdef ORIGINAL_BFMSS
  // To go back to the original BFMSS :
  t->u25() = u_e;
  t->l25() = l_e;
  t->v2p() = t->v2m() = t->v5p() = t->v5m() = EXTLONG_ZERO;
#elif defined BFMSS_2_ONLY
  // To go back to BFMSS[2] only :
  t->u25() = t->u25() + ceilLg5(t->v5p());
  t->l25() = t->l25() + ceilLg5(t->v5m());
  t->v5p() = t->v5m() = EXTLONG_ZERO;
#endif

  if (l_e == EXTLONG_ZERO) {           // no divisions introduced
    t->high() = u_e;
    t->low() = EXTLONG_ONE - u_e; // - (u_e - 1)
  } else {
    t->high() = u_e - l_e + EXTLONG_ONE;
    t->low() = EXTLONG_TWO - t->high();
  }

  t->lc() = l_e;
  t->tc() = u_e;

  // set BigRat value
  if (rationalReduceFlag) {
    t->ratFlag() = 1;
    t->ratValue() = new BigRat(value.BigRatValue());
  }

  t->flagsComputed() = true;
}

void ConstDoubleRep::computeExactFlags() {// can be made more efficient
  computeExactFlags_temp(this, Real(ffVal.getValue()));
}

void ConstRealRep::computeExactFlags() {
  computeExactFlags_temp(this, value);
}

void NegRep::computeExactFlags() {
  if (!child->flagsComputed())
    child->computeExactFlags();

  if (child->sign() == 0) {
    reduceToZero();
    return;
  }

  if (rationalReduceFlag) {
    if (child->ratFlag()>0 && child->ratValue() != NULL) {
      BigRat val = -(*(child->ratValue()));
      reduceToBigRat(val);
      ratFlag() = child->ratFlag()+1;
      return;
    } else
      ratFlag() = -1;
  }

  sign() = -child->sign();
  uMSB() = child->uMSB();
  lMSB() = child->lMSB();

  // length() = child->length();
  measure() = child->measure();
  u25() = child->u25();
  l25() = child->l25();
  v2p() = child->v2p();
  v2m() = child->v2m();
  v5p() = child->v5p();
  v5m() = child->v5m();
  high() = child->high();
  low() = child->low();
  lc()  = child->lc();
  tc()  = child->tc();
  flagsComputed() = true;
}//NegRep::computeExactFlags

void SqrtRep::computeExactFlags() {
  if (!child->flagsComputed())
    child->computeExactFlags();

  if (rationalReduceFlag)
    ratFlag() = -1;

  sign() = child->sign();
  if (sign() < 0)
    core_error("squareroot is called with negative operand.",
               __FILE__, __LINE__, true);

  uMSB() = child->uMSB() / EXTLONG_TWO;
  lMSB() = child->lMSB() / EXTLONG_TWO;

  // length() = child->length();
  measure() = child->measure();

  // BFMSS[2,5] bound.
  if (child->v2p() + ceilLg5(child->v5p()) + child->u25() >=
      child->v2m() + ceilLg5(child->v5m()) + child->l25()) {
    extLong vtilda2 = child->v2p() + child->v2m();
    v2p() = vtilda2 / EXTLONG_TWO;
    v2m() = child->v2m();
    extLong vmod2;
    if (v2p().isInfty())
      vmod2 = CORE_INFTY;
    else
      vmod2 = vtilda2 - EXTLONG_TWO*v2p(); // == vtilda2 % 2
    extLong vtilda5 = child->v5p() + child->v5m();
    v5p() = vtilda5 / EXTLONG_TWO;
    v5m() = child->v5m();
    extLong vmod5;
    if (v5p().isInfty())
      vmod5 = CORE_INFTY;
    else
      vmod5 = vtilda5 - EXTLONG_TWO*v5p(); // == vtilda5 % 2
    u25() = (child->u25() + child->l25() + vmod2 + ceilLg5(vmod5) + EXTLONG_ONE) / EXTLONG_TWO;
    l25() = child->l25();
  } else {
    extLong vtilda2 = child->v2p() + child->v2m();
    v2p() = child->v2p();
    v2m() = vtilda2 / EXTLONG_TWO;
    extLong vmod2;
    if (v2m().isInfty())
      vmod2 = CORE_INFTY;
    else
      vmod2 = vtilda2 - EXTLONG_TWO*v2m(); // == vtilda2 % 2
    extLong vtilda5 = child->v5p() + child->v5m();
    v5p() = child->v5p();
    v5m() = vtilda5 / EXTLONG_TWO;
    u25() = child->u25();
    extLong vmod5;
    if (v5m().isInfty())
      vmod5 = CORE_INFTY;
    else
      vmod5 = vtilda5 - EXTLONG_TWO*v5m(); // == vtilda5 % 2
    l25() = (child->u25() + child->l25() + vmod2 + ceilLg5(vmod5) + EXTLONG_ONE) / EXTLONG_TWO;
  }

  high() = (child->high() +EXTLONG_ONE)/EXTLONG_TWO;
  low() = child->low() / EXTLONG_TWO;
  lc() = child->lc();
  tc() = child->tc();
  flagsComputed() = true;
}// SqrtRep::computeExactFlags

void MultRep::computeExactFlags() {
  if (!first->flagsComputed())
    first->computeExactFlags();
  if (!second->flagsComputed())
    second->computeExactFlags();

  if ((!first->sign()) || (!second->sign())) {
    // value must be exactly zero.
    reduceToZero();
    return;
  }
  // rational node
  if (rationalReduceFlag) {
    if (first->ratFlag() > 0 && second->ratFlag() > 0) {
      BigRat val = (*(first->ratValue()))*(*(second->ratValue()));
      reduceToBigRat(val);
      ratFlag() = first->ratFlag() + second->ratFlag();
      return;
    } else
      ratFlag() = -1;
  }

  // value is irrational.
  uMSB() = first->uMSB() + second->uMSB() + EXTLONG_ONE;
  lMSB() = first->lMSB() + second->lMSB();
  sign() = first->sign() * second->sign();

  extLong df = first->d_e();
  extLong ds = second->d_e();
  // extLong lf = first->length();
  // extLong ls = second->length();

  // length() = df * ls + ds * lf;
  measure() = (first->measure()) * ds+(second->measure()) * df;

  // BFMSS[2,5] bound.
  v2p() = first->v2p() + second->v2p();
  v2m() = first->v2m() + second->v2m();
  v5p() = first->v5p() + second->v5p();
  v5m() = first->v5m() + second->v5m();
  u25() = first->u25() + second->u25();
  l25() = first->l25() + second->l25();

  high() = first->high() + second->high();
  low() = first->low() + second->low();

  lc() = ds * first->lc() + df * second->lc();
  tc() = core_min(ds * first->tc() + df * second->tc(), measure());

  flagsComputed() = true;
}// MultRep::computeExactFlags

void DivRep::computeExactFlags() {
  if (!first->flagsComputed())
    first->computeExactFlags();
  if (!second->flagsComputed())
    second->computeExactFlags();

  if (!second->sign())
    core_error("zero divisor.", __FILE__, __LINE__, true);

  if (!first->sign()) {// value must be exactly zero.
    reduceToZero();
    return;
  }

  // rational node
  if (rationalReduceFlag) {
    if (first->ratFlag() > 0 && second->ratFlag() > 0) {
      BigRat val = (*(first->ratValue()))/(*(second->ratValue()));
      reduceToBigRat(val);
      ratFlag() = first->ratFlag() + second->ratFlag();
      return;
    } else
      ratFlag() = -1;
  }

  // value is irrational.
  uMSB() = first->uMSB() - second->lMSB();
  lMSB() = first->lMSB() - second->uMSB() - EXTLONG_ONE;
  sign() = first->sign() * second->sign();

  extLong df = first->d_e();
  extLong ds = second->d_e();
  // extLong lf = first->length();
  // extLong ls = second->length();

  // length() = df * ls + ds * lf;
  measure() = (first->measure())*ds + (second->measure())*df;

  // BFMSS[2,5] bound.
  v2p() = first->v2p() + second->v2m();
  v2m() = first->v2m() + second->v2p();
  v5p() = first->v5p() + second->v5m();
  v5m() = first->v5m() + second->v5p();
  u25() = first->u25() + second->l25();
  l25() = first->l25() + second->u25();

  high() = first->high() + second->low();
  low() = first->low() + second->high();

  lc() = ds * first->lc() + df * second->tc();
  tc() = core_min(ds * first->tc() + df * second->lc(), measure());

  flagsComputed() = true;
}

//
//  approximation functions
//
void ConstDoubleRep::computeApproxValue(const extLong& /*relPrec*/,
                                        const extLong& /*absPrec*/)
// can ignore precision bounds since ffVal.getValue() returns exact value
{
  appValue() = Real(ffVal.getValue());
}

void ConstRealRep::computeApproxValue(const extLong& relPrec,
                                      const extLong& absPrec) {
  appValue() = value.approx(relPrec, absPrec);
}

void NegRep::computeApproxValue(const extLong& relPrec,
                                const extLong& absPrec) {
  appValue() = -child->getAppValue(relPrec, absPrec);
}

void SqrtRep::computeApproxValue(const extLong& relPrec,
                                 const extLong& absPrec) {
  extLong r = relPrec + relPrec + EXTLONG_EIGHT; // chenli: ???
  extLong a = absPrec + absPrec + EXTLONG_EIGHT;
  extLong pr = - lMSB() + r;
  extLong p  = pr < a ? pr : a;

  Real val = child->getAppValue(r, a);
  if (incrementalEvalFlag) {
    if (appValue() == CORE_REAL_ZERO)
      appValue() = val;
    appValue() = val.sqrt(p, appValue().BigFloatValue());
  } else
    appValue() = val.sqrt(p);
}

void MultRep::computeApproxValue(const extLong& relPrec,
                                 const extLong& absPrec) {
  if (lMSB() < EXTLONG_BIG && lMSB() > EXTLONG_SMALL) {
    extLong r   = relPrec + EXTLONG_FOUR;
    extLong  afr = - first->lMSB() + EXTLONG_ONE;
    extLong  afa = second->uMSB() + absPrec + EXTLONG_FIVE;
    extLong  af  = afr > afa ? afr : afa;
    extLong  asr = - second->lMSB() + EXTLONG_ONE;
    extLong  asa = first->uMSB() + absPrec + EXTLONG_FIVE;
    extLong  as  = asr > asa ? asr : asa;
    appValue() = first->getAppValue(r, af)*second->getAppValue(r, as);
  } else {
    std::cerr << "lMSB = " << lMSB() << std::endl;
    core_error("a huge lMSB in MulRep", __FILE__, __LINE__, false);
  }
}

void DivRep::computeApproxValue(const extLong& relPrec,
                                const extLong& absPrec) {
  if (lMSB() < EXTLONG_BIG && lMSB() > EXTLONG_SMALL) {
    extLong rr  = relPrec + EXTLONG_SEVEN;		// These rules come from
    extLong ra  = uMSB() + absPrec + EXTLONG_EIGHT;	// Koji's Master Thesis, page 65
    extLong ra2 = core_max(ra, EXTLONG_TWO);
    extLong r   = core_min(rr, ra2);
    extLong  af  = - first->lMSB() + r;
    extLong  as  = - second->lMSB() + r;

    extLong pr = relPrec + EXTLONG_SIX;
    extLong pa = uMSB() + absPrec + EXTLONG_SEVEN;
    extLong p  = core_min(pr, pa);	// Seems to be an error:
    // p can be negative here!
    // Also, this does not conform to
    // Koji's thesis which has a default
    // relative precision (p.65).

    appValue() = first->getAppValue(r, af).div(second->getAppValue(r, as), p);
  } else {
    std::cerr << "lMSB = " << lMSB() << std::endl;
    core_error("a huge lMSB in DivRep", __FILE__, __LINE__, false);
  }
}

//
// Debug Help Functions
//
void Expr::debug(int mode, int level, int depthLimit) const {
  std::cout << "-------- Expr debug() -----------" << std::endl;
  std::cout << "rep = " << rep << std::endl;
  if (mode == Expr::LIST_MODE)
    rep->debugList(level, depthLimit);
  else if (mode == Expr::TREE_MODE)
    rep->debugTree(level, 0, depthLimit);
  else
    core_error("unknown debugging mode", __FILE__, __LINE__, false);
  std::cout << "---- End Expr debug(): " << std::endl;
}


const std::string ExprRep::dump(int level) const {
  std::ostringstream ost;
  if (level == OPERATOR_ONLY) {
    ost << op();
  } else if (level == VALUE_ONLY) {
    ost << appValue();
  } else if (level == OPERATOR_VALUE) {
    ost << op() << "[val: " << appValue() << "]";
  } else if (level == FULL_DUMP) {
    ost << op()
    << "[val: "  << appValue() << "; "
    << "kp: " << knownPrecision() << "; "
#ifdef CORE_DEBUG
    << "r: " << relPrecision() << "; "
    << "a: " << absPrecision() << "; "
#endif
    << "lMSB: " << lMSB() << "; "
    << "uMSB: " << uMSB() << "; "
    << "sign: " << sign() << "; "
    // << "length: " << length() << "; "
    << "measure: " << measure() << "; "
    << "d_e: " << d_e() << "; "
    << "u25: " << u25() << "; "
    << "l25: " << l25() << "; "
    << "v2p: " << v2p() << "; "
    << "v2m: " << v2m() << "; "
    << "v5p: " << v5p() << "; "
    << "v5m: " << v5m() << "; "
    << "high: " << high() << "; "
    << "low: " << low() << "; "
    << "lc: " << lc() << "; "
    << "tc: " << tc()
    << "]";
  }
  return std::string(ost.str());
  // note that str() return an array not properly terminated!
}


void UnaryOpRep::debugList(int level, int depthLimit) const {
  if (depthLimit <= 0)
    return;
  if (level == Expr::SIMPLE_LEVEL) {
    std::cout << "(" << dump(OPERATOR_VALUE);
    child->debugList(level, depthLimit - 1);
    std::cout << ")";
  } else if (level == Expr::DETAIL_LEVEL) {
    std::cout << "(" << dump(FULL_DUMP);
    child->debugList(level, depthLimit - 1);
    std::cout << ")";
  }
}

void UnaryOpRep::debugTree(int level, int indent, int depthLimit) const {
  if (depthLimit <= 0)
    return;
  for (int i = 0; i<indent; i++ )
    std::cout << "  ";
  std::cout << "|_";
  if (level == Expr::SIMPLE_LEVEL)
    std::cout << dump(OPERATOR_VALUE);
  else if (level == Expr::DETAIL_LEVEL)
    std::cout << dump(FULL_DUMP);
  std::cout << std::endl;
  child->debugTree(level, indent + 2, depthLimit - 1);
}

void ConstRep::debugList(int level, int depthLimit) const {
  if (depthLimit <= 0)
    return;
  if (level == Expr::SIMPLE_LEVEL) {
    std::cout << "(" << dump(OPERATOR_VALUE) << ")";
  } else if (level == Expr::DETAIL_LEVEL) {
    std::cout << "(" << dump(FULL_DUMP) << ")";
  }
}

void ConstRep::debugTree(int level, int indent, int depthLimit) const {
  if (depthLimit <= 0)
    return;
  for (int i=0; i<indent; i++)
    std::cout << "  ";
  std::cout << "|_";
  if (level == Expr::SIMPLE_LEVEL)
    std::cout << dump(OPERATOR_VALUE);
  else if (level == Expr::DETAIL_LEVEL)
    std::cout << dump(FULL_DUMP);
  std::cout << std::endl;
}

void BinOpRep::debugList(int level, int depthLimit) const {
  if (depthLimit <= 0 )
    return;
  std::cout << "(";
  if (level == Expr::SIMPLE_LEVEL) {
    std::cout << dump(OPERATOR_VALUE);
  } else if (level == Expr::DETAIL_LEVEL) {
    std::cout << dump(FULL_DUMP);
  }
  first->debugList(level, depthLimit - 1);
  std::cout << ", ";
  second->debugList(level, depthLimit - 1);
  std::cout << ")" ;
}

void BinOpRep::debugTree(int level, int indent, int depthLimit) const {
  if (depthLimit <= 0)
    return;
  for (int i=0; i<indent; i++)
    std::cout << "  ";
  std::cout << "|_";
  if (level == Expr::SIMPLE_LEVEL) {
    std::cout << dump(OPERATOR_VALUE);
  } else if (level == Expr::DETAIL_LEVEL) {
    std::cout << dump(FULL_DUMP);
  }
  std::cout << std::endl;
  first->debugTree(level, indent + 2, depthLimit - 1);
  second->debugTree(level, indent + 2, depthLimit - 1);
}

} //namespace CORE
