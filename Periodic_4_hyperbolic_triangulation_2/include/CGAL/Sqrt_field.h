// to compile with boost headers
// g++ main.cpp -I /opt/local/include -o main

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQRT_FIELD_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQRT_FIELD_H

#include <complex>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <assert.h>
//#include <unordered_set>
#include <string>
#include <fstream>

//#include <boost/optional.hpp>

using namespace std;

template<class Int>
class Sqrt_field 
{
public:
  typedef Int Integer;

  Int l;
  Int r;

public:
  Sqrt_field(Int l_ = 0, Int r_ = 0) : l(l_), r(r_) {}
  
  Sqrt_field operator + (const Sqrt_field& rh) const
  {
    Sqrt_field temp = *this;
    temp += rh;

    return temp;
  } 

  Sqrt_field& operator += (const Sqrt_field& rh)
  {
    l += rh.l;
    r += rh.r;
    return *this;
  }
  
  Sqrt_field& operator -= (const Sqrt_field& rh)
  {
    l -= rh.l;
    r -= rh.r;
    return *this;
  }

  Sqrt_field operator - (const Sqrt_field& rh) const
  {
    return Sqrt_field( l - rh.l, r - rh.r );
  }

  Sqrt_field operator - () const
  {
    return Sqrt_field(-l, -r);
  }

  Sqrt_field operator * (const Sqrt_field& rh) const
  {
    Sqrt_field temp = *this;
    temp *= rh;

    return temp;
  }

  Sqrt_field& operator *= (const Sqrt_field& rh)
  {
    Int l1 = l * rh.l + 2*r*rh.r;
    Int r1 = r*rh.l + l*rh.r;
    l = l1;
    r = r1;

    return *this; 
  }

  Sqrt_field abs() const
  {
    if (l >= 0 && r >= 0) {
      return *this;
    }

    if (l <= 0 && r <= 0) {
      return Sqrt_field(-l, -r); 
    }
  
    if(l*l*l - 2*r*r*l >= 0) {
      return *this;
    }
    return Sqrt_field(-l, -r);
  }
};

template<class Int>
bool operator < (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
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
  
  /*
  if (lh.l < rh.l) {
    return true;
  }
  
  if (lh.l == rh.l) {
    if (lh.r < rh.r) {
      return true;
    }
  }*/
  
  return false;
}
/*Seems correct!
template<class Int>
bool operator > (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  if( (lh.l - rh.l)*(lh.l - rh.l) > 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) >= 0 && (rh.r - lh.r) >= 0) {
      return true;
    }
  }
  
  if( (lh.l - rh.l)*(lh.l - rh.l) < 2*(rh.r - lh.r)*(rh.r - lh.r) ) {
    if ((lh.l - rh.l) <= 0 && (rh.r - lh.r) <= 0) {
      return true;
    }
  }
  return false;
}*/

template<class Int>
bool operator >= (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
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
bool operator == (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  return (lh.l == rh.l && lh.r == rh.r);
}

template<class Int>
bool operator != (const Sqrt_field<Int>& lh, const Sqrt_field<Int>& rh)
{
  return !(lh == rh);
}

template<class Int>
Sqrt_field<Int> operator * (const Int& val, const Sqrt_field<Int>& rh)
{
  Sqrt_field<Int> temp = rh;
  temp.l *= val;
  temp.r *= val;

  return temp;
} 

template<class Int>
ostream& operator<<(ostream& os, const Sqrt_field<Int>& nb)
{
  os << nb.l << " " << nb.r;
  return os;
}

// begin MT - added to please libc++
template<class Int>
bool isnan(const Sqrt_field<Int>& x)
{
  return false; // MT TBD
}

template<class Int>
bool isinf(const Sqrt_field<Int>& x)
{
  return false; // MT TBD
}

template<class Int>
Sqrt_field<Int> copysign(const Sqrt_field<Int>& x, const Sqrt_field<Int>& y)
{
  if (y.abs()==y)
    { return x;}
  return -x;
}
// end MT - added to please libc++

#endif // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQRT_FIELD_H

