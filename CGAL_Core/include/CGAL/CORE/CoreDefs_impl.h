/****************************************************************************
 * Core Library Version 1.7, August 2004
 * Copyright (c) 1995-2004 Exact Computation Project
 * All rights reserved.
 *
 * This file is part of CGAL (www.cgal.org).
 *
 * File: CoreDefs.cpp
 * Synopsis:
 *         Useful parameters for Core Library which users may change
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
 * SPDX-License-Identifier: LGPL-3.0-or-later
 ***************************************************************************/

#include "CGAL/CORE/CoreDefs.h"

#include <atomic>

namespace CORE {

//  Default Values

/* ************************************************************
 * ERROR FLAGS
 * ************************************************************ */

#ifndef CGAL_HEADER_ONLY

/** I/O error flag (Default value 0, indicating no error)
 *  User's responsibility to check and reset value to 0. */
// This is currently used in geom2d/points2d.cpp for I/O of points

// Note from 2014: does not seem to be used anywhere, and it is not declared
// in CoreDefs.h so it is not accessible
// Left here for compatibilty when CGAL_HEADER_ONLY is not defined

int IOErrorFlag = 0;

/**
 * If AbortFlag is true when invalid expression is constructed, system will abort
 */
#ifdef CGAL_NO_ATOMIC
bool AbortFlag = true;
#else
std::atomic<bool> AbortFlag(true);
#endif

/**
 * InvalidFlag is set to negative whenever an invalid expression is constructed.
 * The user has the responsibility to reset to non-negative value.
 */

#ifdef CGAL_NO_ATOMIC
int InvalidFlag = 0;
#else
std::atomic<int> InvalidFlag(0);
#endif

/* ************************************************************
 * PRECISION PARAMETERS
 * ************************************************************ */

/**
 *  Default BigFloat Division Relative Precision
 *  -- this is used by BigFloat division when the arguments are error-free.
 */

extLong defBFdivRelPrec = 54;

/**
 *  Default BigFloat Sqrt Absolute Precision
 *  -- this is used by BigFloat sqrt when the argument is error-free.
 */

extLong defBFsqrtAbsPrec = 54;

/**
 * Escape Precision
 *   -- we will not compare a number to precision higher than this
 *   -- if this is infinity, there there is no escape precision */
extLong EscapePrec  = CORE_posInfty;

/** this flag becomes negative if the EscapePrec is used. */
long EscapePrecFlag = 0;

/// Escape Precision Warning Flag
/** this flag is true by default, and will cause a warning to be printed
    when EscapePrec is reached */
#ifdef CGAL_NO_ATOMIC
bool EscapePrecWarning = true;
#else
std::atomic<bool> EscapePrecWarning(true);
#endif

/** The Composite Precision [defAbsPrec, defRelPrec]
 *  determines the precision to which an Expr evaluates its
 *  (exact, implicit) constant value. */

/**  absolute precision  = 2^31 - 1 */
extLong defAbsPrec = CORE_posInfty;
/** default relative precision is 60 relative bits.
 *  Why 60?  We would really like this to be 54, so that the default
 *  conversion duplicates the IEEE double precision.  But it turns out
 *  (see README file under BUGS), we need 59 to ensure this.
 *  Chee Yap (7/1/01) */
extLong defRelPrec = 60;

/**  number of BigFloat digits to print out */
#ifdef CGAL_NO_ATOMIC
long defBigFloatOutputDigits = 10;
#else
std::atomic<long> defBigFloatOutputDigits(10);
#endif

/**  NORMALLY, we like to make this equal to defBigFloatOutputDigits
  *  8/3/01, Chee: re-introduced this parameter */
#ifdef CGAL_NO_ATOMIC
long defOutputDigits = 10;
#else
std::atomic<long> defOutputDigits(10); // == defBigFloatOutputDigits;
#endif

/** String Input Precision */

/** Set this to 16 if you want machine precision. This controls the
 *  absolute error in string decimal inputs to Real or Expr.
 *  If defInputDigits is finite, then the absolute error will be
 *  at most 10^{-defInputDigits}.  Otherwise, the input is exactly
 *  represented by some BigFloat or BigRat value. */
extLong defInputDigits = CORE_posInfty;

/** This controls the absolute error in converting from string to BigFloat
 *  The absolute error will be at most 10^{-defInputDigits} */
#ifdef CGAL_NO_ATOMIC
long defBigFloatInputDigits = 16;
#else
std::atomic<long> defBigFloatInputDigits(16);
#endif

/* ************************************************************
 * EVALUATION FLAGS
 * ************************************************************ */

/** Floating Point filter
 *  true = turn on floating point filter */
#ifdef CGAL_NO_ATOMIC
bool fpFilterFlag = true;
#else
std::atomic<bool> fpFilterFlag(true);
#endif

/** IncrementaL evaluation flag
 *  incremental evaluation is requested, This means, we try to use previous
 *  approximate values to improve an approximation */
#ifdef CGAL_NO_ATOMIC
bool incrementalEvalFlag = true;
#else
std::atomic<bool> incrementalEvalFlag(true);
#endif

/** Progressive evaluation flag
 *  true = turn on progressive evaluation flag */
#ifdef CGAL_NO_ATOMIC
bool progressiveEvalFlag = true;
#else
std::atomic<bool> progressiveEvalFlag(true);
#endif

/** Initial progressive evaluation precision
 *  Used by AddSubRep */
#ifdef CGAL_NO_ATOMIC
long defInitialProgressivePrec = 64;
#else
std::atomic<long> defInitialProgressivePrec(64);
#endif

/** RATIONAL REDUCTION FLAG
 *  true = turn on rational reduction */
#ifdef CGAL_NO_ATOMIC
bool rationalReduceFlag = false;
#else
std::atomic<bool> rationalReduceFlag(false);
#endif
#endif // CGAL_HEADER_ONLY

} //namespace CORE

