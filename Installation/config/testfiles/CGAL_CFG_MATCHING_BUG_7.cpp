// Copyright (c) 2014  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     :  Andreas Fabri, Laurent Rineau

//| This flag is set if the compiler bugs with some "using Base::Member;" in
//| a derived class.  The workaround is to write a forwarder or not use
//| using.

enum Type { POINT, WPOINT };

#include <cassert>

int I;

struct Point {
};


struct WPoint
{
  WPoint(const Point&) {} // Point is implicitly convertible to WPoint.
};


struct Conv
{
  Conv() {}

  template <typename T>
  Type operator()(const T&){ return POINT; }
};


struct WConv : Conv
{
  WConv() {}

  Type operator()(const WPoint&){ return WPOINT; }

#if MATCHING_BUG_8 // only defined by CGAL_CFG_MATCHING_BUG_8.cpp
  Type operator()(const Point&){ return POINT; }
#endif

  using Conv::operator(); // Import the operator() of Conv, that matches
                          // for every type.

  // Workaround: comment the previous 'using' line, and write:
  //   template <typename T>
  //   Type operator()(const T& t){ return Conv::operator()(t); }
};


int main()
{
  I = 0;
  Point p;
  WPoint wp(p);
  WConv wconv;
  if(wconv(p) != POINT) return 1;
  if(wconv(wp) != WPOINT) return 1;
  return 0;
}


// VC++11 and VC++12 both say that the call to `wconv(p)` is ambiguous:
//
// config\testfiles\CGAL_CFG_MATCHING_BUG_7.cpp(61) : error C2666: 'WConv::operator ()' : 2 overloads have similar conversions
// 
//         config\testfiles\CGAL_CFG_MATCHING_BUG_7.cpp(51): could be 'void WConv::operator ()(const WPoint &)'
// 
//         unable to recover from previous error(s); stopping compilation
