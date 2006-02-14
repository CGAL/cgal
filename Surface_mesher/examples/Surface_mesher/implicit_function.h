// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve OUDOT


#ifndef _IMPLICIT_FUNCTION_H
#define _IMPLICIT_FUNCTION_H


int implicit_function (double, double, double);


template <class FT> class Function {
public:  

  enum SURFACE_LOCATION {IN = -1, ON = 0, OUT = 1};
  
  SURFACE_LOCATION operator()(FT x, FT y, FT z) const { 
    int res = implicit_function (CGAL::to_double (x), 
				 CGAL::to_double (y), 
				 CGAL::to_double (z));
    
    switch (res) {
    case -1: return IN;
    case 0: return ON;
    case 1: return OUT;
    default: assert (res == -1 || res == 0 || res == 1);
    }
    return OUT;  // never used
  }
};
  

#endif
