// Copyright (c) 2010   INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQUARE_ROOT_2_FIELD_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQUARE_ROOT_2_FIELD_H

#include <complex>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <assert.h>
#include <string>
#include <fstream>

using namespace std;

template<class Int>
class Square_root_2_field 
{

private:
  typedef Int Integer;

public:
  Int l;
  Int r;

  Square_root_2_field(Int l_ = 0, Int r_ = 0) : l(l_), r(r_) {}
  
  Square_root_2_field operator + (const Square_root_2_field& rh) const
  {
    Square_root_2_field temp = *this;
    temp += rh;

    return temp;
  } 

  Square_root_2_field& operator += (const Square_root_2_field& rh)
  {
    l += rh.l;
    r += rh.r;
    return *this;
  }
  
  Square_root_2_field& operator -= (const Square_root_2_field& rh)
  {
    l -= rh.l;
    r -= rh.r;
    return *this;
  }

  Square_root_2_field operator - (const Square_root_2_field& rh) const
  {
    return Square_root_2_field( l - rh.l, r - rh.r );
  }

  Square_root_2_field operator - () const
  {
    return Square_root_2_field(-l, -r);
  }

  Square_root_2_field operator * (const Square_root_2_field& rh) const
  {
    Square_root_2_field temp = *this;
    temp *= rh;

    return temp;
  }

  Square_root_2_field& operator *= (const Square_root_2_field& rh)
  {
    Int l1 = l * rh.l + 2*r*rh.r;
    Int r1 = r*rh.l + l*rh.r;
    l = l1;
    r = r1;

    return *this; 
  }

  Square_root_2_field abs() const
  {
    if (l >= 0 && r >= 0) {
      return *this;
    }

    if (l <= 0 && r <= 0) {
      return Square_root_2_field(-l, -r); 
    }

    if(l*l*l - 2*r*r*l >= 0) {
      return *this;
    }
    return Square_root_2_field(-l, -r);
  }
};

template<class Int>
bool operator < (const Square_root_2_field<Int>& lh, const Square_root_2_field<Int>& rh)
{
  if( (lh.l - rh.l)*(lh.l - rh.l) > 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) < 0) {
      return true;
    }
  }
  
  if( (lh.l - rh.l)*(lh.l - rh.l) < 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((rh.r - lh.r) >= 0) {
      return true;
    }
  }
  return false;
}


template<class Int>
bool operator >= (const Square_root_2_field<Int>& lh, const Square_root_2_field<Int>& rh)
{
  if( (lh.l - rh.l)*(lh.l - rh.l) >= 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) >= 0) {
      return true;
    }
  }
  
  if( (lh.l - rh.l)*(lh.l - rh.l) < 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((rh.r - lh.r) < 0) {
      return true;
    }
  }
  return false;
}

template<class Int>
bool operator == (const Square_root_2_field<Int>& lh, const Square_root_2_field<Int>& rh)
{
  return (lh.l == rh.l && lh.r == rh.r);
}

template<class Int>
bool operator != (const Square_root_2_field<Int>& lh, const Square_root_2_field<Int>& rh)
{
  return !(lh == rh);
}

template<class Int>
Square_root_2_field<Int> operator * (const Int& val, const Square_root_2_field<Int>& rh)
{
  Square_root_2_field<Int> temp = rh;
  temp.l *= val;
  temp.r *= val;

  return temp;
} 

template<class Int>
ostream& operator<<(ostream& os, const Square_root_2_field<Int>& nb)
{
  os << nb.l << " " << nb.r;
  return os;
}

// begin MT - added to please libc++
template<class Int>
bool isnan(const Square_root_2_field<Int>& x)
{
  return false; // MT TBD
}

template<class Int>
bool isinf(const Square_root_2_field<Int>& x)
{
  return false; // MT TBD
}

template<class Int>
Square_root_2_field<Int> copysign(const Square_root_2_field<Int>& x, const Square_root_2_field<Int>& y)
{
  if (y.abs()==y)
    { return x;}
  return -x;
}
// end MT - added to please libc++

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQUARE_ROOT_2_FIELD_H

