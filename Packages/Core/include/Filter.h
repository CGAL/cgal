/******************************************************************
 * Core Library Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: Filter.h
 *
 * Synopsis:
 *      This is a simple filtered floating point number,
 *      represented by the main class, FilterFp.
 *      based on the Burnikel-Funke-Schirra (BFS) filter scheme.
 *      We do not use IEEE exception mechanism here.
 *      It is used by the Expr class.
 *
 * Written by 
 *       Zilin Du <zilin@cs.nyu.edu>
 *       Chee Yap <yap@cs.nyu.edu>
 *
 * WWW URL: http://cs.nyu.edu/exact/
 * Email: exact@cs.nyu.edu
 *
 * $Id$
 *****************************************************************/

#ifndef CORE_FILTER_H
#define CORE_FILTER_H

#include "CoreImpl.h"
#include "Real.h"

#if defined (_MSC_VER) || defined (__MINGW32__) // add support for MinGW
  #define finite(x)	_finite(x)
  #define ilogb(x)	(int)_logb(x)
#else
  extern "C" int finite(double);	// since SunOS defined "finite" in 
  extern "C" int ilogb(double);		// <ieeefp.h>, but gnu defined it
					// in <math.h>, so we declared it
					// here explictly.
					// Zilin Du: 07/18/2002
#endif

CORE_BEGIN_NAMESPACE

const int POWTWO_26 = (1 << 26);  ///< constant 2^26

#ifdef _MSC_VER
#pragma warning(disable: 4056) // warning message "overflow in floating-point constant arithmetic"
#endif
const double DBL_INFTY = 2*DBL_MAX;  ///< constant infty double
#ifdef _MSC_VER
#pragma warning(default: 4056)
#endif

/// \class filteredFp Filter.h
/// \brief filteredFp represents filtered floating point
///        numbers, based on BFS filter
class filteredFp {
  double fpVal;         // approximate double value for some "real value"
  double maxAbs;        // if (|fpVal| > maxAbs * ind * 2^{-53}) then
  int ind;              // sign of value is sign(fpVal).  Else, don't know.
                        // REFERENCE: Burnikel, Funke, Schirra (BFS) filter
public:
  /// \name Constructors
  //@{
  /// constructor
  filteredFp (double val = 0.0) 
    : fpVal(val), maxAbs(core_abs(val)), ind(0) {}
  /// constructor
  filteredFp (double val, double m, int i) 
    : fpVal(val), maxAbs(m), ind(i) {}
  
  /// construct a filteredFp from Real v.
  /** if v causes an overflow, fpVal = +/- Infty
      if v causes an underflow, fpVal = ...? */
  filteredFp (const Real & value) : fpVal(0.0), maxAbs(0.0), ind(0) {
    if (value != CORE_REAL_ZERO) {
      ind = 1;
      fpVal = value.toDouble();
      maxAbs = core_abs(fpVal); // NaN are propagated correctly by core_abs.
    }
  }
  //@}

  /// \name Help Functions
  //@{
  /// return filtered value (for debug)
  double getValue() const { return fpVal; }
  /// check whether filtered value is OK
  bool isOK() const 
  { return (fpFilterFlag  && // To disable filter
	    finite(fpVal) && // Test for infinite and NaNs
	    (core_abs(fpVal) >= maxAbs*ind*CORE_EPS)); }
  /// return the sign of fitered value.
  /** (note: call isOK() to check whether the sign is ok 
      before call this function.) */
  int sign() const {
#ifdef DEBUG
    assert(isOK());
#endif
    if (fpVal == 0.0)
      return 0;
    else
      return fpVal > 0.0 ? 1: -1;
  }
  /// lower bound on MSB
  /** defined to be cel(lg(real value));
      ilogb(x) is floor(log_2(|x|)). 
      Also, ilogb(0) = -INT_MAX.
   	    ilogb(NaN) = ilogb(+/-Inf) = INT_MAX */
  extLong lMSB() const 
  { return extLong(ilogb(core_abs(fpVal)-maxAbs*ind*CORE_EPS));}
  /// upper bound on MSB 
  extLong uMSB() const 
  { return extLong(ilogb(core_abs(fpVal)+maxAbs*ind*CORE_EPS));}
  //@}

  /// \name Operators
  //@{
  /// unary minus
  filteredFp operator -() const
  { return filteredFp(-fpVal, maxAbs, ind); }
  /// addition
  filteredFp operator+ (const filteredFp& x) const
  { return filteredFp(fpVal+x.fpVal, maxAbs+x.maxAbs, 1+core_max(ind, x.ind)); }
  /// subtraction 
  filteredFp operator- (const filteredFp& x) const
  { return filteredFp(fpVal-x.fpVal, maxAbs+x.maxAbs, 1+core_max(ind, x.ind)); }
  /// multiplication
  filteredFp operator* (const filteredFp& x) const
  { return filteredFp(fpVal*x.fpVal, maxAbs*x.maxAbs+DBL_MIN, 1+ind+x.ind); }
  /// division
  filteredFp operator/ (const filteredFp& x) const {
    if (x.fpVal == 0.0)
      core_error("possible zero divisor!", __FILE__, __LINE__, false);
    double xxx = core_abs(x.fpVal) / x.maxAbs - (x.ind+1)*CORE_EPS + DBL_MIN;
    if (xxx > 0) {
      double val =  fpVal / x.fpVal;
      double maxVal = ( core_abs(val) + maxAbs / x.maxAbs) / xxx + DBL_MIN;
      return filteredFp(val, maxVal, 1 + core_max(ind, x.ind + 1));
    } else
      return filteredFp(DBL_INFTY, 0.0, 0);
  }
  /// square root 
  filteredFp sqrt () const {
    if (fpVal < 0.0)
      core_error("possible negative sqrt!", __FILE__, __LINE__, false);    
    if (fpVal > 0.0) {
      double val = ::sqrt(fpVal);
      return filteredFp(val,  ( maxAbs / fpVal ) * val, 1 + ind);
    } else 
      return filteredFp(0.0, ::sqrt(maxAbs) * POWTWO_26, 1 + ind);
  }

  void dump (std::ostream&os) const {
    os << " Filter = [fpVal = " << fpVal << " , maxAbs = " << maxAbs
	    << " , ind = " << ind << " ]"; 
  }
  //@}
}; //filteredFp class

inline
std::ostream & operator<< (std::ostream & os, const filteredFp& fp)
{
  fp.dump(os);
  return os;
}

CORE_END_NAMESPACE

#endif

