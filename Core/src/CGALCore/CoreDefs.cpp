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
 * File: CoreDefs.cpp
 * Synopsis:
 *	 Useful parameters for Core Library which users may change
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

#include "CGAL/CORE/CoreDefs.h"

namespace CORE { 

//  Default Values

/* ************************************************************
 * ERROR FLAGS
 * ************************************************************ */

/** I/O error flag (Default value 0, indicating no error)
 *  User's responsibility to check and reset value to 0. */
// This is currently used in geom2d/points2d.cpp for I/O of points

int IOErrorFlag = 0;

/**
 * If AbortFlag is true when invalid expression is constructed, system will abort
 */

bool AbortFlag = true;

/**
 * InvalidFlag is set to negative whenever an invalid expression is constructed.
 * The user has the responsibility to reset to non-negative value.
 */

int InvalidFlag = 0;

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
bool EscapePrecWarning = true;

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
long defBigFloatOutputDigits = 10;

/**  NORMALLY, we like to make this equal to defBigFloatOutputDigits
  *  8/3/01, Chee: re-introduced this parameter */
long defOutputDigits = defBigFloatOutputDigits;

/** String Input Precision */

/** Set this to 16 if you want machine precision. This controls the
 *  absolute error in string decimal inputs to Real or Expr.
 *  If defInputDigits is finite, then the absolute error will be 
 *  at most 10^{-defInputDigits}.  Otherwise, the input is exactly 
 *  represented by some BigFloat or BigRat value. */
extLong defInputDigits = CORE_posInfty;

/** This controls the absolute error in converting from string to BigFloat
 *  The absolute error will be at most 10^{-defInputDigits} */
long defBigFloatInputDigits = 16;

/* ************************************************************
 * EVALUATION FLAGS
 * ************************************************************ */

/** Floating Point filter
 *  true = turn on floating point filter */
bool fpFilterFlag = true;

/** IncrementaL evaluation flag
 *  incremental evaluation is requested, This means, we try to use previous
 *  approximate values to improve an approximation */
bool incrementalEvalFlag = true;

/** Progressive evaluation flag
 *  true = turn on progressive evaluation flag */
bool progressiveEvalFlag = true;

/** Initial progressive evaluation precision
 *  Used by AddSubRep */
long defInitialProgressivePrec = 64;

/** RATIONAL REDUCTION FLAG
 *  true = turn on rational reduction */
bool rationalReduceFlag = false;

} //namespace CORE

