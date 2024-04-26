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

  void set_real_part(const FT& real_part);
  void set_imaginary_part(const FT& imaginary_part);

  FT real_part() const;
  FT imaginary_part() const;

  FT squared_modulus() const;
  _Self conjugate() const;
  _Self operator+(const _Self& other) const;
  _Self operator-(const _Self& other) const;
  _Self operator-() const;
  _Self operator*(const _Self& other) const;
  _Self operator/(const _Self& other) const;
};

template<class FT> bool operator==(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2);
template<class FT> bool operator!=(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2);

template<class FT>std::ostream& operator<<(std::ostream& s, const Complex_without_sqrt<FT>& z);
template<class FT>void operator>>(std::istream& s, Complex_without_sqrt<FT>& z);

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
void Complex_without_sqrt<FT>::set_real_part(const FT& real_part){
  _real = real_part;
}

template<class FT>
void Complex_without_sqrt<FT>::set_imaginary_part(const FT& imaginary_part){
  _imag = imaginary_part;
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
FT Complex_without_sqrt<FT>::real_part() const{
  return _real;
}

template<class FT>
FT Complex_without_sqrt<FT>::imaginary_part() const{
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

template<class FT>
Complex_without_sqrt<FT> Complex_without_sqrt<FT>::operator+(const Complex_without_sqrt<FT>& other) const{
  return Complex_without_sqrt<FT>(_real+other.real_part(), _imag+other.imaginary_part());
}

template<class FT>
Complex_without_sqrt<FT> Complex_without_sqrt<FT>::operator-(const Complex_without_sqrt<FT>& other) const{
  return Complex_without_sqrt<FT>(_real-other.real_part(), _imag-other.imaginary_part());
}

template<class FT>
Complex_without_sqrt<FT> Complex_without_sqrt<FT>::operator-() const{
  return Complex_without_sqrt<FT>(-_real, -_imag);
}

template<class FT>
Complex_without_sqrt<FT> Complex_without_sqrt<FT>::operator*(const Complex_without_sqrt<FT>& other) const{
  return Complex_without_sqrt<FT>(_real*other.real_part()-_imag*other.imaginary_part(), _real*other.imaginary_part()+_imag*other.real_part());
}

template<class FT>
Complex_without_sqrt<FT> Complex_without_sqrt<FT>::operator/(const Complex_without_sqrt<FT>& other) const{
  FT m2 = other.squared_modulus();
  return Complex_without_sqrt<FT>(_real/m2, _imag/m2)*other.conjugate();
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
bool operator==(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2){
  return (z1.real_part()==z2.real_part() && z1.imaginary_part()==z2.imaginary_part());
}

template<class FT>
bool operator!=(const Complex_without_sqrt<FT>& z1, const Complex_without_sqrt<FT>& z2){
  return !operator==(z1, z2);
}

////////////////////////////////////////////////////////////////////////////////

template<class FT>
std::ostream& operator<<(std::ostream& s, const Complex_without_sqrt<FT>& z){
  s << z.real_part() << std::endl << z.imaginary_part() << std::endl;
  return s;
}

template<class FT>
void operator>>(std::istream& s, Complex_without_sqrt<FT>& z){
  std::string line;
  std::getline(s, line);
  z.set_real_part(FT(line));
  std::getline(s, line);
  z.set_imaginary_part(FT(line));
}

} // namespace CGAL

#endif // CGAL_COMPLEX_WITHOUT_SQRT
