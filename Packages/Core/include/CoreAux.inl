/******************************************************************
 * Core Library, Version 1.5, August 2002
 * Copyright (c) 1995-2002 Exact Computation Project
 * 
 * File: CoreAux.inl
 * $Id$ 
******************************************************************/
#ifdef CORE_INLINE

////////////////////////////////////////////////////////////
//  More useful functions to implement:
//
//  To convert digits into bits:
//      given X, compute X * log_2(10)
//  To convert bits into digits:
//      given X, compute X * log_10(2) 
//
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// flrLg(x)
//      returns floor log base 2 of abs(x)
// CONVENTION lg(0) = -1	(Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
CORE_INLINE int flrLg(long x) {
  if (x == LONG_MIN) {
    // special treatment as -LONG_MIN would be not representable as "long"
    return LONG_BIT - 1;
  } else {
    //  1 <= |x| <= LONG_MAX
    if (x < 0) x = - x;
    
    int lg = -1;
    while (x > 0) {lg++; x >>= 1;}
    return lg;
  }
}

////////////////////////////////////////////////////////////
// floor log base 2 of unsigned long x
// CONVENTION lg(0) = -1	(Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
CORE_INLINE int flrLg(unsigned long x) {
  int lg = -1;
  while (x > 0) {lg++; x >>= 1;}
  return lg;
}

////////////////////////////////////////////////////////////
// ceiling log base 2 of abs(x)
// CONVENTION lg(0) = -1	(Slight improvement, Zilin/Chee 8/5/02)
////////////////////////////////////////////////////////////
CORE_INLINE int clLg(long x) {
  if (x == LONG_MIN) return LONG_BIT - 1;
  if (x < 0) x = -x; 		// use absolute value
  if (x > (LONG_MAX >> 1)) 	// the leading data bit is 1
      return (LONG_BIT - 1);	// exclude the sign bit
  if (x >= 2) return flrLg((unsigned long)((x << 1) - 1));
    				// SINCE ceilLog_2(x) = floorLog_2(2x-1) for x>=2
  if (x == 1) return 0;
  return -1;			// x must be 0 here
}

////////////////////////////////////////////////////////////
// ceiling log base 2 of unsigned long x
// CONVENTION lg(0) = -1
////////////////////////////////////////////////////////////
CORE_INLINE int clLg(unsigned long x) {
  if (x > (ULONG_MAX >> 1))	// the leading bit is 1
    return LONG_BIT;
  if (x >= 2) return flrLg((x << 1) - 1); 
  				// SINCE ceilLog_2(x) = floorLog_2(2x-1) for x>=2.
  if (x == 1) return 0;
  return -1;	// x must be equal to 0
}

CORE_INLINE std::ostream& operator<< (std::ostream& o, const std::string& s ) {
  o << s.c_str();
  return o;
}

CORE_INLINE void core_error(std::string msg, std::string filename, int lineno, bool ex) {
  std::cerr << "CORE " << (ex? "ERROR" : "WARNING") << " (at " << filename 
	    << ":" << lineno << "): " << msg << std::endl;
  if (ex) exit(1);
}

#endif
