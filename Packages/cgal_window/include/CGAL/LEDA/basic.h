// Copyright (c) 1999  Martin-Luther-University Halle-Wittenberg (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Matthias Baesken, Algorithmic Solutions



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
