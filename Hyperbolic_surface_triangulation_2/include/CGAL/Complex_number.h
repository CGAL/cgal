// Copyright (c) 2024
// INRIA Nancy (France), and Université Gustave Eiffel Marne-la-Vallée (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Loïc Dubois, Monique Teillaud

// This file contains the declaration and the implementation of the class Complex_number

#ifndef CGAL_COMPLEX_NUMBER
#define CGAL_COMPLEX_NUMBER

#include <fstream>
#include <sstream>

namespace CGAL {
/*
Templated by a field FT. Represents a complex number over FT.
*/
template <class FT>
class Complex_number {
private:
  typedef Complex_number<FT>            _Self;
  FT _real, _imag;

public:
  Complex_number();
  Complex_number(const FT& real_part);
  Complex_number(const FT& real_part, const FT& imaginary_part);
  template<class U,class V>
    Complex_number(U&& real_part, V&& imaginary_part): _real(std::forward<U>(real_part)), _imag(std::forward<V>(imaginary_part)) {}

  void real(const FT& real_part);
  void imag(const FT& imaginary_part);
  FT real() const;
  FT imag() const;
  FT squared_modulus() const;
  _Self conjugate() const;

  _Self& operator+=(const _Self& other);
  _Self& operator-=(const _Self& other);
  _Self& operator*=(const _Self& other);
  _Self& operator/=(const _Self& other);
  _Self& operator=(const _Self& other);

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

  friend _Self operator+(const _Self& z1, const _Self& z2){
    return _Self(z1._real+z2._real, z1._imag+z2._imag);
  }

  friend _Self operator-(const _Self& z1, const _Self& z2) {
  return _Self(z1._real-z2._real, z1._imag-z2._imag);
  }

  friend _Self operator*(const _Self& z1, const _Self& z2) {
    return _Self(z1._real*z2._real-z1._imag*z2._imag, z1._real*z2._imag+z1._imag*z2._real);
  }

  friend _Self operator/(const _Self& z1, const _Self& z2){
    FT m2 = z2.squared_modulus();
    return _Self(z1._real/m2, z1._imag/m2)*z2.conjugate();
  }

  friend std::ostream& operator<<(std::ostream& s, const _Self& z){
    s << z._real << std::endl << z._imag << std::endl;
    return s;
  }

  friend void operator>>(std::istream& s, _Self& z){
    std::string line;
    s >> line;
    z.real(FT(line));
    s >> line;
    z.imag(FT(line));
  }

};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class FT>
Complex_number<FT>::Complex_number(){
  _real = FT(0);
  _imag = FT(0);
}

template<class FT>
Complex_number<FT>::Complex_number(const FT& real_part){
  _real = real_part;
  _imag = FT(0);
}

template<class FT>
Complex_number<FT>::Complex_number(const FT& real_part, const FT& imaginary_part){
  _real = real_part;
  _imag = imaginary_part;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
void Complex_number<FT>::real(const FT& real_part){
  _real = real_part;
}

template<class FT>
void Complex_number<FT>::imag(const FT& imaginary_part){
  _imag = imaginary_part;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
FT Complex_number<FT>::real() const{
  return _real;
}

template<class FT>
FT Complex_number<FT>::imag() const{
  return _imag;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
FT Complex_number<FT>::squared_modulus() const{
  return _real*_real + _imag*_imag;
}

template<class FT>
Complex_number<FT> Complex_number<FT>::conjugate() const{
  return Complex_number<FT>(_real, -_imag);
}

////////////////////////////////////////////////////////////////////////////////
template<class FT>
Complex_number<FT>& Complex_number<FT>::operator+=(const Complex_number<FT>& other) {
  _real += other.real();
  _imag += other.imag();
  return *this;
}

template<class FT>
Complex_number<FT>& Complex_number<FT>::operator-=(const Complex_number<FT>& other) {
  _real -= other.real();
  _imag -= other.imag();
  return *this;
}

 template<class FT>
Complex_number<FT>& Complex_number<FT>::operator*=(const Complex_number<FT>& other) {
  _real = _real*other.real() - _imag*other.imag();
  _imag = _real*other.imag() + _imag*other.real();
  return *this;
}

 template<class FT>
Complex_number<FT>& Complex_number<FT>::operator/=(const Complex_number<FT>& other) {
   FT m2 = other.squared_modulus();
   _real /= m2;
   _imag /= m2;
   this *= other.conjugate();
   return *this;
}

 template<class FT>
Complex_number<FT>& Complex_number<FT>::operator=(const Complex_number<FT>& other) {
  _real = other.real();
  _imag = other.imag();
  return *this;
}

} // namespace CGAL

#endif // CGAL_COMPLEX_NUMBER
