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

// This file contains the declaration and the implementation of the class Complex_without_sqrt

#ifndef CGAL_COMPLEX_WITHOUT_SQRT
#define CGAL_COMPLEX_WITHOUT_SQRT

#include <fstream>
#include <sstream>

namespace CGAL {
/*
Templated by a field FT. Represents a complex number over FT.
*/
template <class FT>
class Complex_without_sqrt {
private:
  typedef Complex_without_sqrt<FT>            _Self;
  FT _real, _imag;

public:
  Complex_without_sqrt();
  Complex_without_sqrt(const FT& real_part);
  Complex_without_sqrt(const FT& real_part, const FT& imaginary_part);

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

  template<class FT_>
  friend Complex_without_sqrt<FT_> operator+(const Complex_without_sqrt<FT_>& z);
  template<class FT_>
  friend Complex_without_sqrt<FT_> operator-(const Complex_without_sqrt<FT_>& z);
  template<class FT_>
  friend bool operator==(const Complex_without_sqrt<FT_>& z1, const Complex_without_sqrt<FT_>& z2);
  template<class FT_>
  friend bool operator!=(const Complex_without_sqrt<FT_>& z1, const Complex_without_sqrt<FT_>& z2);
  template<class FT_>
  friend Complex_without_sqrt<FT_> operator+(const Complex_without_sqrt<FT_>& z1, const Complex_without_sqrt<FT_>& z2);
  template<class FT_>
  friend Complex_without_sqrt<FT_> operator-(const Complex_without_sqrt<FT_>& z1, const Complex_without_sqrt<FT_>& z2);
  template<class FT_>
  friend Complex_without_sqrt<FT_> operator*(const Complex_without_sqrt<FT_>& z1, const Complex_without_sqrt<FT_>& z2);
  template<class FT_>
  friend Complex_without_sqrt<FT_> operator/(const Complex_without_sqrt<FT_>& z1, const Complex_without_sqrt<FT_>& z2);

  template<class FT_>
  friend std::ostream& operator<<(std::ostream& s, const Complex_without_sqrt<FT_>& z);
  template<class FT_>
  friend void operator>>(std::istream& s, Complex_without_sqrt<FT_>& z);
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template<class FT>
Complex_without_sqrt<FT>::Complex_without_sqrt(){
  _real = FT(0);
  _imag = FT(0);
}

template<class FT>
Complex_without_sqrt<FT>::Complex_without_sqrt(const FT& real_part){
  _real = real_part;
  _imag = FT(0);
}

template<class FT>
Complex_without_sqrt<FT>::Complex_without_sqrt(const FT& real_part, const FT& imaginary_part){
  _real = real_part;
  _imag = imaginary_part;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
void Complex_without_sqrt<FT>::real(const FT& real_part){
  _real = real_part;
}

template<class FT>
void Complex_without_sqrt<FT>::imag(const FT& imaginary_part){
  _imag = imaginary_part;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
FT Complex_without_sqrt<FT>::real() const{
  return _real;
}

template<class FT>
FT Complex_without_sqrt<FT>::imag() const{
  return _imag;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
FT Complex_without_sqrt<FT>::squared_modulus() const{
  return _real*_real + _imag*_imag;
}

template<class FT>
Complex_without_sqrt<FT> Complex_without_sqrt<FT>::conjugate() const{
  return Complex_without_sqrt<FT>(_real, -_imag);
}

////////////////////////////////////////////////////////////////////////////////
template<class FT>
Complex_without_sqrt<FT>& Complex_without_sqrt<FT>::operator+=(const Complex_without_sqrt<FT>& other) {
  _real += other.real();
  _imag += other.imag();
  return *this;
}

template<class FT>
Complex_without_sqrt<FT>& Complex_without_sqrt<FT>::operator-=(const Complex_without_sqrt<FT>& other) {
  _real -= other.real();
  _imag -= other.imag();
  return *this;
}

 template<class FT>
Complex_without_sqrt<FT>& Complex_without_sqrt<FT>::operator*=(const Complex_without_sqrt<FT>& other) {
  _real = _real*other.real() - _imag*other.imag();
  _imag = _real*other.imag() + _imag*other.real();
  return *this;
}

 template<class FT>
Complex_without_sqrt<FT>& Complex_without_sqrt<FT>::operator/=(const Complex_without_sqrt<FT>& other) {
   FT m2 = other.squared_modulus();
   _real /= m2;
   _imag /= m2;
   this *= other.conjugate();
   return *this;
}

 template<class FT>
Complex_without_sqrt<FT>& Complex_without_sqrt<FT>::operator=(const Complex_without_sqrt<FT>& other) {
  _real = other.real();
  _imag = other.imag();
  return *this;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
Complex_without_sqrt<FT> operator+(const Complex_without_sqrt<FT>& z) {
  return z;
}

template<class FT>
Complex_without_sqrt<FT> operator-(const Complex_without_sqrt<FT>& z) {
  return Complex_without_sqrt<FT>(-z._real,-z._imag);
}

template<class FT>
bool operator==(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2){
  return (z1._real==z2._real && z1._imag==z2._imag);
}

template<class FT>
bool operator!=(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2){
  return !operator==(z1, z2);
}

template<class FT>
Complex_without_sqrt<FT> operator+(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2) {
  return Complex_without_sqrt<FT>(z1._real+z2._real, z1._imag+z2._imag);
}

template<class FT>
Complex_without_sqrt<FT> operator-(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2) {
  return Complex_without_sqrt<FT>(z1._real-z2._real, z1._imag-z2._imag);
}

template<class FT>
Complex_without_sqrt<FT> operator*(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2) {
  return Complex_without_sqrt<FT>(z1._real*z2._real-z1._imag*z2._imag, z1._real*z2._imag+z1._imag*z2._real);
}

template<class FT>
Complex_without_sqrt<FT> operator/(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2) {
  FT m2 = z2.squared_modulus();
  return Complex_without_sqrt<FT>(z1._real/m2, z2._imag/m2)*z2.conjugate();
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
std::ostream& operator<<(std::ostream& s, const Complex_without_sqrt<FT>& z){
  s << z._real << std::endl << z._imag << std::endl;
  return s;
}

template<class FT>
void operator>>(std::istream& s, Complex_without_sqrt<FT>& z){
  std::string line;
  s >> line;
  z.real(FT(line));
  s >> line;
  z.imag(FT(line));
}

} // namespace CGAL

#endif // CGAL_COMPLEX_WITHOUT_SQRT
