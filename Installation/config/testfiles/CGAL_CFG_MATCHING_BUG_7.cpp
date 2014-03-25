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
//| a derived class.  The workaround is to write a forwarder or not use using.

#include <cassert>

int I;

struct Point {
};


struct WP
{
  WP(const Point& p)
  {}

};


struct Conv
{
  Conv() {}

  template <typename T>
  void operator()(const T& t){ I = 1;}
};


struct WConv : Conv
{
  using Conv::operator();

  WConv() {}

  void operator()(const WP& w){I=2;}

};


int main()
{
  I = 0;
  Point p;

  WConv wconv;

  wconv(p);

  assert(I==1);
  return 0;
}
