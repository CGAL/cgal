/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
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
 * $Id$
 *****************************************************************/

#ifndef CORE_DEFS_H
#define CORE_DEFS_H

#include "CoreImpl.h"
#include "extLong.h"

CORE_BEGIN_NAMESPACE

//////////////////////////////////////////////////////////////
// defined constants
//////////////////////////////////////////////////////////////

/// default accuracy level
#define DEFAULT_LEVEL 3		

/// short hand for positive infinity
#define CORE_INFTY  (CORE_posInfty) 

//////////////////////////////////////////////////////////////
// global precision parameters 
//////////////////////////////////////////////////////////////

/// default BigFloat Division Relative Precision 
extern extLong defBFdivRelPrec;

/// default BigFloat Sqrt Absolute Precision 
extern extLong defBFsqrtAbsPrec;

/// Escape Precision in bits
extern extLong EscapePrec;
/// current ur when EscapePrec triggered
/** this flag becomes negative when default EscapePrec is applied */
extern long EscapePrecFlag;

// These following two values determine the precision of computing 
// approximations in Expr.

/// default Relative Precision in bits
extern extLong defRelPrec;
/// default Absolute Precision in bits
extern extLong defAbsPrec;

/// default # of decimal digits for conversion from a BF to string.
extern long defBigFloatOutputDigits;

/// default input precision in digits for converting a string to a Real or Expr
/** This value can be CORE_INFTY */
extern extLong defInputDigits; 
					
/// controls the printout precision of std::cout for Real and Expr
/** This value cannot be CORE_INFTY 
 *  (it really should be an int, as in std::cout.setprecision(int)). */
extern long defOutputDigits;

/// default input precision in digits for converting a string to a BigFloat
/** This value cannot be CORE_INFTY */
extern long defBigFloatInputDigits;

/// controls the printout precision of std::cout for BigFloat
/** This value cannot be CORE_INFTY */
extern long defBigFloatOutputDigits;

//////////////////////////////////////////////////////////////
// Mode parameters: incremental, progressive, filters
//////////////////////////////////////////////////////////////

/// floating point filter flag
extern bool fpFilterFlag; 
/// if true, evaluation of expressions would be incremental
extern bool incrementalEvalFlag;
/// progressive evaluation flag
extern bool progressiveEvalFlag; 
/// rational reduction flag
extern bool rationalReduceFlag; 

//////////////////////////////////////////////////////////////
// methods for setting global precision parameters 
// 	including: scientific vs. positional format
//	All the set methods return the previous global value if any
//////////////////////////////////////////////////////////////

/// set default composite precision [defAbsPrec, defRelPrec] 
/** It determines the precision to which an Expr evaluates its 
    (exact, implicit) constant value. */
CORE_INLINE void setDefaultPrecision(const extLong &r, const extLong &a);

/// set default relative precision
CORE_INLINE extLong setDefaultRelPrecision(const extLong &r);

/// set default absolute precision
CORE_INLINE extLong setDefaultAbsPrecision(const extLong &a);


/// set default input digits (for Expr, Real)
/** it controls the absolute error */
CORE_INLINE extLong setDefaultInputDigits(const extLong &d);

/// set default output digits (for Expr, Real)
CORE_INLINE long setDefaultOutputDigits(long d = defOutputDigits,
		std::ostream& o = std::cout);

/// set default input digits for BigFloat
CORE_INLINE long setDefaultBFInputDigits(long d);

/// set default output digits for BigFloat
CORE_INLINE long setDefaultBFOutputDigits(long d);

/// turn floating-point filter on/off
CORE_INLINE bool setFpFilterFlag(bool f);

/// turn incremental evaluation flag on/off
CORE_INLINE bool setIncrementalEvalFlag(bool f);

/// turn progressive evaluation flag on/off
CORE_INLINE bool setProgressiveEvalFlag(bool f);

/// turn rational reduction flag on/off
CORE_INLINE bool setRationalReduceFlag(bool f);

/// CORE_init(..) is the CORE initialization function.
/** We recommend calling it before anything else.  Originally motivated
    by need to get around gnu's compiler bug in which the variable 
    "defAbsPrec" was not properly initialized.  But it has other uses,
    e.g., overriding the default std::cout precision to our own */
CORE_INLINE void CORE_init(long d);

/// change to scientific output format
CORE_INLINE void setScientificFormat(std::ostream& o = std::cout);

/// change to positional output format
CORE_INLINE void setPositionalFormat(std::ostream& o = std::cout);

#ifdef CORE_ENABLE_INLINES
#include "CoreDefs.inl"
#endif

CORE_END_NAMESPACE
#endif
