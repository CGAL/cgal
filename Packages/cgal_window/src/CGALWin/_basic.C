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

