// Copyright (c) 2016-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     :   Iordan Iordanov
//

#ifndef CGAL_EXACT_COMPLEX_H
#define CGAL_EXACT_COMPLEX_H

#include <CGAL/license/Hyperbolic_triangulation_2.h>

#include <CGAL/assertions.h>
#include <CGAL/number_utils.h>

#include <iostream>
#include <fstream>

// Complex number in the form a + bi, where a and b are of type NT.
// NT must be an exact number type, model of:
//  + FieldWithRootOf
//  + RealEmbeddable
//  + FromDoubleConstructible
//

namespace CGAL {

template <class NumberType>
class Exact_complex
{
  typedef Exact_complex<NumberType>       Self;

private:
  NumberType _a, _b;

public:
  typedef NumberType NT;

  Exact_complex() : _a(0), _b(0) {}
  Exact_complex(NT a, NT b) : _a(a), _b(b) {}

  NT real() const { return _a; }
  void set_real(NT val) { _a = val; }

  NT imag() const { return _b; }
  void set_imag(NT val) { _b = val; }

  Self conj() const { return Self(_a, -_b); }
  NT square_modulus() const { return (_a*_a + _b*_b); }
  NT modulus() const { return CGAL::sqrt(this->square_modulus()); }

  Self reciprocal()
  {
    NT denom = _a*_a + _b*_b;
    if(denom == NT(0))
      return Self(0,0);
    else
      return Self(_a/denom, -_b/denom);
  }

  Self invert_in_unit_circle() { return this->reciprocal().conj(); }

  template <class Circle_2>
  Self invert_in_circle(Circle_2 c)
  {
    NT r2 = c.squared_radius();
    NT xc = c.center().x();
    NT yc = c.center().y();
    NT denom = (_a - xc)*(_a - xc) + (_b - yc)*(_b - yc);
    return Self(xc + r2*(_a-xc)/denom, yc + r2*(_b-yc)/denom);
  }

  Self operator+(const Self& other) const
  {
    return Self(this->_a + other._a, this->_b + other._b);
  }

  Self operator-(const Self& other) const
  {
    return Self(this->_a - other._a, this->_b - other._b);
  }

  Self operator*(const Self& other) const
  {
    NT rp = _a*other._a - _b*other._b;
    NT ip = _b*other._a + _a*other._b;
    return Self(rp, ip);
  }

  Self operator/(const Self& other) const
  {
    NT denom = other._a*other._a + other._b*other._b;
    CGAL_assertion(denom != 0);
    NT rp = _a*other._a + _b*other._b;
    NT ip = _b*other._a - _a*other._b;
    return Self(rp/denom, ip/denom);
  }
};

template <class NT>
std::ostream& operator<<(std::ostream& s, const Exact_complex<NT>& c)
{
  s << c.real() << (c.imag() >= 0 ? " + " : " - ") << abs(c.imag()) << "i";
  return s;
}

// just to give an ordering
template<class NT>
bool operator==(const Exact_complex<NT>& lh, const Exact_complex<NT>& rh)
{
  return (lh.real() == rh.real() && lh.imag() == rh.imag());
}

template<class NT>
bool operator!=(const Exact_complex<NT>& lh, const Exact_complex<NT>& rh)
{
  return !operator==(lh, rh);
}

// just to give an ordering
template<class NT>
bool operator<(const Exact_complex<NT>& lh, const Exact_complex<NT>& rh)
{
  return lh.square_modulus() < rh.square_modulus();
}

} // namespace CGAL

#endif // CGAL_EXACT_COMPLEX_H
