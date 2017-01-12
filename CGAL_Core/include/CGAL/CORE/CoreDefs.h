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
 * File: CoreDefs.h
 * Synopsis:
 *       This contains useful Core Library global parameters which
 *       users may modify at runtime or compile time
 *       For each parameter, we provide corresponding methods to
 *       modify or examine the values.
 *
 * Written by 
 *       Chee Yap <yap@cs.nyu.edu>
 *       Chen Li <chenli@cs.nyu.edu>
 *       Zilin Du <zilin@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $URL$
 * $Id$
 ***************************************************************************/

#ifndef _CORE_COREDEFS_H_
#define _CORE_COREDEFS_H_

#include <CGAL/CORE/extLong.h>
#include <CGAL/atomic.h>

#ifdef CGAL_HEADER_ONLY
  
  #define CGAL_GLOBAL_STATE_VAR(TYPE, NAME, VALUE)  \
    inline TYPE & get_static_##NAME()               \
    {                                               \
      static TYPE NAME(VALUE);                      \
      return NAME;                                  \
    }
  
#else // CGAL_HEADER_ONLY

  #define CGAL_GLOBAL_STATE_VAR(TYPE, NAME, VALUE)  \
    CGAL_CORE_EXPORT extern TYPE NAME;                   \
    inline TYPE& get_static_##NAME()                \
    {                                               \
      return NAME;                                  \
    }

#endif // CGAL_HEADER_ONLY

namespace CORE { 

//////////////////////////////////////////////////////////////
// defined constants
//////////////////////////////////////////////////////////////

/// default accuracy level
#define DEFAULT_CORE_LEVEL 3

/// short hand for positive infinity
#define CORE_INFTY  (CORE_posInfty)

//////////////////////////////////////////////////////////////
// global precision parameters
//////////////////////////////////////////////////////////////

/// Abort Flag -- default value is true
/** The normal behavior is to abort when an invalid expression
 * is constructed.  This flag can be used to turn off this abort.
 * In any case, an error message will be printed */
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(bool, AbortFlag, true)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<bool>, AbortFlag, true)
#endif

/// Invalid Flag -- initiallly value is non-negative
/** If the Abort Flag is false, then the Invalid flag will be set to
 *  a negative value whenever an invalid expression is constructed.
 *  It is the user's responsibility to check this flag and to make
 *  it non-negative again. */
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(int, InvalidFlag, 0)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<int>, InvalidFlag, 0)
#endif

/// Escape Precision in bits
CGAL_GLOBAL_STATE_VAR(extLong, EscapePrec, CORE_posInfty)

/// current ur when EscapePrec triggered
/** this flag becomes negative when default EscapePrec is applied */
CGAL_GLOBAL_STATE_VAR(long, EscapePrecFlag, 0)

/// Escape Precision Warning Flag
/** this flag is true by default, and will cause a warning to be printed
    when EscapePrec is reached */
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(bool, EscapePrecWarning, true)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<bool>, EscapePrecWarning, true)
#endif

// These following two values determine the precision of computing
// approximations in Expr.

/// default Relative Precision in bits
CGAL_GLOBAL_STATE_VAR(extLong, defRelPrec, 60)
/// default Absolute Precision in bits
CGAL_GLOBAL_STATE_VAR(extLong, defAbsPrec, CORE_posInfty)

/// default # of decimal digits for conversion from a BF to string.
/** This value cannot be CORE_INFTY.
    See also defOutputDigits. 
    */
/*  QUESTION: the following comment seems to contradict the above comment:
	"controls the printout precision of std::cout for BigFloat"
    Perhaps, we should merge defOutputDigits and defBigFloatOutputDigits?
    */
#ifdef CGAL_NO_ATOMIC
    CGAL_GLOBAL_STATE_VAR(long, defBigFloatOutputDigits, 10)
#else
    CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<long>, defBigFloatOutputDigits, 10)
#endif

/// default input precision in digits for converting a string to a Real or Expr
/** This value can be CORE_INFTY */
CGAL_GLOBAL_STATE_VAR(extLong, defInputDigits, CORE_posInfty)

/// controls the printout precision of std::cout for Real and Expr
/** This value cannot be CORE_INFTY
    See also defBigFloatOutputDigits. 
    (it really should be an int, as in std::cout.setprecision(int)). */
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(long, defOutputDigits, 10) // == get_static_defBigFloatOutputDigits()
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<long>, defOutputDigits, 10) // == get_static_defBigFloatOutputDigits()
#endif

/// default input precision in digits for converting a string to a BigFloat
/** This value cannot be CORE_INFTY. */
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(long, defBigFloatInputDigits, 16)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<long>, defBigFloatInputDigits, 16)
#endif

inline
long
getBigFloatInputDigits()
{
  return get_static_defBigFloatInputDigits();
}

/// default BigFloat Division Relative Precision
CGAL_GLOBAL_STATE_VAR(extLong, defBFdivRelPrec, 54)
/// default BigFloat Sqrt Absolute Precision
CGAL_GLOBAL_STATE_VAR(extLong, defBFsqrtAbsPrec, 54)

//////////////////////////////////////////////////////////////
// Mode parameters: incremental, progressive, filters
//////////////////////////////////////////////////////////////

/// floating point filter flag
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(bool, fpFilterFlag, true)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<bool>, fpFilterFlag, true)
#endif


/// if true, evaluation of expressions would be incremental
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(bool, incrementalEvalFlag, true)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<bool>, incrementalEvalFlag, true)
#endif


/// progressive evaluation flag
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(bool, progressiveEvalFlag, true)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<bool>, progressiveEvalFlag, true)
#endif


/// rational reduction flag
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(bool, rationalReduceFlag, false)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<bool>, rationalReduceFlag, false)
#endif

/// default initial (bit) precision for AddSub Progressive Evaluation
#ifdef CGAL_NO_ATOMIC
CGAL_GLOBAL_STATE_VAR(long, defInitialProgressivePrec, 64)
#else
CGAL_GLOBAL_STATE_VAR(CGAL::cpp11::atomic<long>, defInitialProgressivePrec, 64)
#endif

//////////////////////////////////////////////////////////////
// methods for setting global precision parameters
// 	including: scientific vs. positional format
//	All the set methods return the previous global value if any
//////////////////////////////////////////////////////////////

/// set default composite precision [defAbsPrec, defRelPrec]
/** It determines the precision to which an Expr evaluates its
    (exact, implicit) constant value. */
inline void setDefaultPrecision(const extLong &r, const extLong &a) {
  get_static_defRelPrec() = r;
  get_static_defAbsPrec() = a;
}

/// set default relative precision
inline extLong setDefaultRelPrecision(const extLong &r) {
  extLong old = get_static_defRelPrec();
  get_static_defRelPrec() = r;
  return old;
}

/// set default absolute precision
inline extLong setDefaultAbsPrecision(const extLong &a) {
  extLong old = get_static_defAbsPrec();
  get_static_defAbsPrec() = a;
  return old;
}

/// set default input digits (for Expr, Real)
/** it controls the absolute error */
inline extLong setDefaultInputDigits(const extLong &d) {
  extLong old = get_static_defInputDigits();
  get_static_defInputDigits() = d;
  return old;
}

/// set default output digits (for Expr, Real)
inline long setDefaultOutputDigits(long d = get_static_defOutputDigits(),
                                   std::ostream& o = std::cout) {
  long old = get_static_defOutputDigits();
  get_static_defOutputDigits() = d;
  o.precision(d);
  return old;
}

/// set default input digits for BigFloat
inline long setDefaultBFInputDigits(long d) {
  long old = get_static_defBigFloatInputDigits();
  get_static_defBigFloatInputDigits() = d;
  return old;
}

/// set default output digits for BigFloat
inline long setDefaultBFOutputDigits(long d) {
  long old = get_static_defBigFloatOutputDigits();
  get_static_defBigFloatOutputDigits() = d;
  return old;
}

/// turn floating-point filter on/off
inline bool setFpFilterFlag(bool f) {
  bool oldf = get_static_fpFilterFlag();
  get_static_fpFilterFlag() = f;
  return oldf;
}

/// turn incremental evaluation flag on/off
inline bool setIncrementalEvalFlag(bool f) {
  bool oldf = get_static_incrementalEvalFlag();
  get_static_incrementalEvalFlag() = f;
  return oldf;
}

/// turn progressive evaluation flag on/off
inline bool setProgressiveEvalFlag(bool f) {
  bool oldf = get_static_progressiveEvalFlag();
  get_static_progressiveEvalFlag() = f;
  return oldf;
}

  inline long getInitialProgressivePrec()
  {
    return get_static_defInitialProgressivePrec();
  }

/// set initial bit precision for progressive evaluation:
inline long setDefInitialProgressivePrec(long n) {
  long oldn = get_static_defInitialProgressivePrec();
  get_static_defInitialProgressivePrec() = n;
  return oldn;
}

/// turn rational reduction flag on/off
inline bool setRationalReduceFlag(bool f) {
  bool oldf = get_static_rationalReduceFlag();
  get_static_rationalReduceFlag() = f;
  return oldf;
}

/// CORE_init(..) is the CORE initialization function.
/** We recommend calling it before anything else.  Originally motivated
    by need to get around gnu's compiler bug in which the variable 
    "defAbsPrec" was not properly initialized.  But it has other uses,
    e.g., overriding the default std::cout precision (most systems 
    initializes this value to 6) to our own */
inline void CORE_init(long d) {
  get_static_defAbsPrec() = CORE_posInfty;
  get_static_defOutputDigits() = d;
  std::setprecision(get_static_defOutputDigits());
}

/// change to scientific output format
inline void setScientificFormat(std::ostream& o = std::cout) {
  o.setf(std::ios::scientific, std::ios::floatfield);
}

/// change to positional output format
inline void setPositionalFormat(std::ostream& o = std::cout) {
  o.setf(std::ios::fixed, std::ios::floatfield);
}

} //namespace CORE

#ifdef CGAL_HEADER_ONLY
#include <CGAL/CORE/CoreDefs_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // _CORE_COREDEFS_H_
