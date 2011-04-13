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
// file          : src/CGALWin/_basic.C
// package       : cgal_window (1.0)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9.7
// revision_date : 23 May 2001
// author(s)     : Matthias Baesken, Algorithmic Solutions
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>) 
// ======================================================================



#include <CGAL/LEDA/basic.h>


#if defined(__unix__)

#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>


#if defined(_AIX)
#include<strings.h>
#include<sys/select.h>
#endif


namespace CGAL {

void leda_wait(double sec) 
{ int usec = int(1000000*sec);
  timeval delay;
  delay.tv_sec  = usec / 1000000;
  delay.tv_usec = usec % 1000000;
  select(0, NULL, NULL, NULL, &delay); 
}

}

#else

#include <windows.h>

using namespace std;

namespace CGAL {

void leda_wait(double sec) 
{ int msec = int(1000*sec);
  Sleep(msec);   
}

}

#endif

namespace CGAL {

double truncate(double x, int k)
{ if ( x == 0 ) return 0;
  if ( k >= 52 ) return x;
  int exp;
  double result = frexp(x,&exp);
  if ( k <= 30 )
  { result = result * (1 << k);
    result = (int) result;
    result = result / (1 << k);  
    return ldexp(result,exp);
  }
  result = result * (1 << 30);
  int result1 = (int ) result;
  result = result - result1;   // bit 31 to ...
  result = result1 + truncate(result,k - 30);
  result = result / (1 << 30);
  return ldexp(result,exp);
}

}

