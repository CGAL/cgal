/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: CoreDefs.inl
 * $Id$ 
******************************************************************/

#ifdef CORE_INLINE

//////////////////////////////////////////////////////////////
// The Composite Precision [defAbsPrec, defRelPrec] 
//    determines the precision to which an Expr evaluates its 
//    (exact, implicit) constant value.
//////////////////////////////////////////////////////////////

CORE_INLINE extLong setDefaultRelPrecision(const extLong &r) {
  extLong old = defRelPrec;
  defRelPrec = r;
  return (old);
}

CORE_INLINE extLong setDefaultAbsPrecision(const extLong &a) {
  extLong old = defAbsPrec;
  defAbsPrec = a;
  return (old);
}

CORE_INLINE void setDefaultPrecision(const extLong &r, const extLong &a) {
  defRelPrec = r;
  defAbsPrec = a;
  return;
}

//////////////////////////////////////////////////////////////
// Output precision
//////////////////////////////////////////////////////////////

CORE_INLINE long setDefaultBFOutputDigits(long d) {
  long old = defBigFloatOutputDigits;
  defBigFloatOutputDigits = d;
  return (old);
}

CORE_INLINE long setDefaultOutputDigits(long d, std::ostream& o) {
  long old = defOutputDigits;
  defOutputDigits = d;
  o.precision(d);  
  return (old);
}

//////////////////////////////////////////////////////////////
// String Input Precision
//////////////////////////////////////////////////////////////

CORE_INLINE long setDefaultBFInputDigits(long d) {
  long old = defBigFloatInputDigits;
  defBigFloatInputDigits = d;
  return (old);
}

CORE_INLINE extLong setDefaultInputDigits(const extLong &d) {
  extLong old = defInputDigits;
  defInputDigits = d;		//  this controls the absolute error
  return (old);
}

//////////////////////////////////////////////////////////////
// FLOATING POINT FILTER
//////////////////////////////////////////////////////////////

CORE_INLINE bool setFpFilterFlag(bool f) {  // Set filter flag, return old flag
  bool oldf = fpFilterFlag;
  fpFilterFlag = f;
  return oldf;
}

//////////////////////////////////////////////////////////////
// INCREMENTAL EVALUATION FLAGS
//////////////////////////////////////////////////////////////

CORE_INLINE bool setIncrementalEvalFlag(bool f) {
  bool oldf = incrementalEvalFlag;
  incrementalEvalFlag = f;
  return oldf;
}

//////////////////////////////////////////////////////////////
// PROGRESSIVE EVALUATION FLAGS
//////////////////////////////////////////////////////////////

CORE_INLINE bool setProgressiveEvalFlag(bool f) {
  bool oldf = progressiveEvalFlag;
  progressiveEvalFlag = f;
  return oldf;
}

CORE_INLINE long setDefInitialProgressivePrec(long n) {
  long oldn = defInitialProgressivePrec;
  defInitialProgressivePrec = n;
  return oldn;
}

//////////////////////////////////////////////////////////////
// RATIONAL REDUCTION FLAGS
//////////////////////////////////////////////////////////////

CORE_INLINE bool setRationalReduceFlag(bool f) {
  bool oldf = rationalReduceFlag;
  rationalReduceFlag = f;
  return oldf;
}

////////////////////////////////////////////////////////////
// CORE_init(..) is the CORE initialization function.
//	We recommend calling it before anything else.
//	Originally motivated by need to get around gnu's compiler bug
//	in which the variable "defAbsPrec" was not properly initialized
//	But it has other uses,
//	e.g., overriding the default std::cout precision to our own
////////////////////////////////////////////////////////////

CORE_INLINE void CORE_init(long defOutputDig) {
  defAbsPrec = CORE_posInfty;
  defOutputDigits = defOutputDig;
  		// : most systems initializes this value to 6
  std::cout << std::setprecision(defOutputDigits);	
}

// no associated global parameter:
CORE_INLINE void setScientificFormat(std::ostream & o) {
  o.setf(std::ios::scientific, std::ios::floatfield);
}

CORE_INLINE void setPositionalFormat(std::ostream & o) {
  o.setf(std::ios::fixed, std::ios::floatfield);
}

#endif

