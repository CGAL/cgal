// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallée (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Marc Pouget, Monique Teillaud

#ifndef CGAL_COMPLEX_NUMBER_H
#define CGAL_COMPLEX_NUMBER_H

#include <fstream>
#include <sstream>
#include <utility>

namespace CGAL {

/*
Templated by a field FT. Represents a complex number over FT.
*/
template <class FT>
class Complex_number
{
  typedef Complex_number<FT> _Self;

  FT _real, _imag;

public:
  Complex_number(const FT& real_part)
    : _real(real_part)
  {}

  Complex_number(const FT& real_part, const FT& imaginary_part)
    : _real(real_part),
      _imag(imaginary_part)
  {}

  Complex_number()
    : Complex_number(FT(0), FT(0))
  {}

  template<class U,class V>
  Complex_number(U&& real_part, V&& imaginary_part)
    : _real(std::forward<U>(real_part)),
      _imag(std::forward<V>(imaginary_part))
  {}

  void real(const FT& real_part) {
    _real = real_part;
  }

  void imag(const FT& imaginary_part) {
    _imag = imaginary_part;
  }

  FT real() const {
    return _real;
  }

  FT imag() const {
    return _imag;
  }

  _Self& operator+=(const _Self& other);
  _Self& operator-=(const _Self& other);
  _Self& operator*=(const _Self& other);
  _Self& operator/=(const _Self& other);

  // These member versions are not working ?
  /* _Self operator+(const _Self& z) const; */
  /* _Self operator-(const _Self& z) const; */

  // Hidden friends
  friend _Self operator+(const _Self& z) {
    return z;
  }

  friend _Self operator-(const _Self& z) {
    return _Self(-z._real,-z._imag);
  }

  friend bool operator==(const _Self& z1, const _Self& z2) {
    return (z1._real==z2._real && z1._imag==z2._imag);
  }

  friend bool operator!=(const _Self& z1, const _Self& z2) {
    return !operator==(z1, z2);
  }

  friend _Self operator+(const _Self& z1, const _Self& z2) {
    return _Self(z1._real+z2._real, z1._imag+z2._imag);
  }

  friend _Self operator-(const _Self& z1, const _Self& z2) {
    return _Self(z1._real-z2._real, z1._imag-z2._imag);
  }

  friend _Self operator*(const _Self& z1, const _Self& z2) {
    return _Self(z1._real*z2._real-z1._imag*z2._imag, z1._real*z2._imag+z1._imag*z2._real);
  }

  friend _Self operator/(const _Self& z1, const _Self& z2) {
    FT m2 = norm(z2);
    return _Self(z1._real/m2, z1._imag/m2)*conj(z2);
  }

  friend std::ostream& operator<<(std::ostream& s, const _Self& z) {
    s << z._real << std::endl << z._imag << std::endl;
    return s;
  }

  friend void operator>>(std::istream& s, _Self& z) {
    FT ft;
    s >> ft;
    z.real(ft);
    s >> ft;
    z.imag(ft);
  }
};

////////////////////////////////////////////////////////////////////////////////
template<class FT>
Complex_number<FT>& Complex_number<FT>::operator+=(const Complex_number<FT>& other)
{
  _real += other.real();
  _imag += other.imag();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator-=(const Complex_number<FT>& other)
{
  _real -= other.real();
  _imag -= other.imag();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator*=(const Complex_number<FT>& other)
{
  _real = _real*other.real() - _imag*other.imag();
  _imag = _real*other.imag() + _imag*other.real();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator/=(const Complex_number<FT>& other)
{
  FT m2 = norm(other);
  _real /= m2;
  _imag /= m2;
  this *= conj(other);
  return *this;
}

template<class FT>
FT norm(const Complex_number<FT>& z)
{
  return z.real()*z.real() + z.imag()*z.imag();
}

template<class FT>
Complex_number<FT> conj(const Complex_number<FT>& z)
{
  return Complex_number<FT>(z.real(), -z.imag());
}

} // namespace CGAL

#endif // CGAL_COMPLEX_NUMBER_H
