// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.3-I-75 $
// release_date  : $CGAL_Date: 2001/06/21 $
//
// file          : include/CGAL/LEDA/basic.h
// package       : cgal_window (1.0.3)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.0.3
// revision_date : 25 June 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================



#ifndef CGAL_WINDOW_BASIC_H
#define CGAL_WINDOW_BASIC_H


// include system config file
#if defined(CGAL_USE_CGAL_HEADERS)
#include <CGAL/basic.h>
#else

#if !defined(CGAL_CLIB_STD) 
#if defined(_MSC_VER)
#define CGAL_CLIB_STD
#else
#define CGAL_CLIB_STD std
#endif
#endif

#endif


#include <CGAL/LEDA/system.h>

// include std header files

#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cstddef>
#include <cstdlib>
#include <cmath>


// include basic LEDA headers

#include <CGAL/LEDA/global.h>

#include <string>

namespace CGAL {

extern __exportF void leda_wait(double sec);  
/*{\Mfunc  suspends execution for $sec$ seconds.}*/

// maximal and minimal values for some numerical types

inline int    Max_Value(int& x)     { return x =  MAXINT;   }
inline int    Min_Value(int& x)     { return x = -MAXINT;   }
inline double Max_Value(double& x)  { return x =  MAXDOUBLE;}
inline double Min_Value(double& x)  { return x = -MAXDOUBLE;}

extern __exportF double truncate(double x, int k = 10);
/*{\Mfunc returns a double whose mantissa is truncated after $k-1$ bits after the binary point, i.e, if
$x \not= 0$ then the binary representation of the mantissa of the 
result has the form d.dddddddd, where the number of d's is equal to $k$. 
There is a corresponding function for
|integers|; it has no effect.}*/

}

#endif
